import queue
import os
import argparse
import logging
from multiprocessing import Queue, Process
from itertools import count

import pysam
from spliceai.utils import Annotator, get_delta_scores
from collections import namedtuple


try:
    from sys.stdin import buffer as std_in
    from sys.stdout import buffer as std_out
except ImportError:
    from sys import stdin as std_in
    from sys import stdout as std_out


def get_options():

    parser = argparse.ArgumentParser(description='Version: 1.3')
    parser.add_argument('-I', metavar='input', nargs='?', default=std_in,
                        help='path to the input VCF file, defaults to standard in')
    parser.add_argument('-O', metavar='output', nargs='?', default=std_out,
                        help='path to the output VCF file, defaults to standard out')
    parser.add_argument('-R', metavar='reference', required=True,
                        help='path to the reference genome fasta file')
    parser.add_argument('-A', metavar='annotation', required=True,
                        help='"grch37" (GENCODE V24lift37 canonical annotation file in '
                             'package), "grch38" (GENCODE V24 canonical annotation file in '
                             'package), or path to a similar custom gene annotation file')
    parser.add_argument('-D', metavar='distance', nargs='?', default=50,
                        type=int, choices=range(0, 5000),
                        help='maximum distance between the variant and gained/lost splice '
                             'site, defaults to 50')
    parser.add_argument('-M', metavar='mask', nargs='?', default=0,
                        type=int, choices=[0, 1],
                        help='mask scores representing annotated acceptor/donor gain and '
                             'unannotated acceptor/donor loss, defaults to 0')
    parser.add_argument('-P', metavar='processes', default=1, type=int)
    args = parser.parse_args()

    return args


def run_serial(args):
    """
    串行运行
    """
    try:
        vcf = pysam.VariantFile(args.I)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()
    header = vcf.header
    header.add_line('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 variant '
                    'annotation. These include delta scores (DS) and delta positions (DP) for '
                    'acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). '
                    'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')
    try:
        output = pysam.VariantFile(args.O, mode='w', header=header)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()
    ann = Annotator(args.R, args.A)
    for record in vcf:
        scores = get_delta_scores(record, ann, args.D, args.M)
        if len(scores) > 0:
            record.info['SpliceAI'] = scores
        output.write(record)
    vcf.close()
    output.close()


ParallelVariantRecord = namedtuple('ParallelVariantRecord', ['id', 'chrom', 'pos', 'ref', 'alts'])


def process_record(records, results, ref_fasta, annotations, dist_var, mask):
    # 创建一个注释助手类
    ann = Annotator(ref_fasta, annotations)
    # 监听队列
    while True:
        # 尝试从队列获得一个待打分的变异
        try:
            record = records.get(True, 0.05)
        except queue.Empty:
            continue
        # 判断队列是否结束
        if record != 'END':
            # 对变异进行打分并把结果放入队列
            scores = get_delta_scores(record, ann, dist_var, mask)
            results.put((record.id, scores))
        else:
            # 队列结束，重新把结束标志放入队列，以终止其他进程
            records.put('END')
            break


def run_parallel(args):
    """
    并行运行
    """
    # 尝试打开文件
    try:
        vcf = pysam.VariantFile(args.I)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()
    header = vcf.header
    header.add_line('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3 variant '
                    'annotation. These include delta scores (DS) and delta positions (DP) for '
                    'acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). '
                    'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')
    try:
        output = pysam.VariantFile(args.O, mode='w', header=header)
    except (IOError, ValueError) as e:
        logging.error('{}'.format(e))
        exit()

    # 计算可用的核心数
    processes_num = min(os.cpu_count(), args.P)
    
    # 创建队列，队列长度10倍于进程数，既避免队列过大占用内存，又避免占用读取IO阻塞
    records, results = Queue(10 * processes_num), Queue()
    
    # 记录已经进入队列，等待运行运行结果的变异
    waiting_records = dict()
    # 生成记录ID
    record_ids = count()
    # 运行标记
    input_finished = False
    output_finished = False
    
    # 启动新进程，监听队列
    for i in range(processes_num):
        Process(target=process_record, args=(records, results, args.R, args.A, args.D, args.M)).start()

    while True:
        # 向队列中加入变异
        while not input_finished and not records.full():
            try:
                # 取得变异
                record_id, record = next(record_ids), next(vcf)
                # 将变异放入队列
                # pysam自带的VariantRecord对象无法在多个进程中共享，使用自己模仿的对象替换
                records.put(ParallelVariantRecord(record_id, record.chrom, record.pos, record.ref, record.alts))
                waiting_records[record_id] = record
            except StopIteration:
                # 读取结束，进行记录并放入队列结束标记
                input_finished = True
                records.put('END')
                break
        
        # 处理获得的结果
        while waiting_records:
            # 尝试从队列获得结果
            try:
                result = results.get(True, 0.05)
            except queue.Empty:
                break
            # 打分结果和变异对应
            record = waiting_records.pop(result[0])
            # 输出打分后的变异
            if len(result[1]) > 0:
                record.info['SpliceAI'] = result[1]
            output.write(record)
        else:
            # 标记队列结束
            output_finished = True
        
        # 全部处理完成
        if output_finished:
            break

    vcf.close()
    output.close()


def main():

    args = get_options()

    if None in [args.I, args.O, args.D, args.M]:
        logging.error('Usage: spliceai [-h] [-I [input]] [-O [output]] -R reference -A annotation '
                      '[-D [distance]] [-M [mask]]')
        exit()
    
    # 根据参数的核心数来选择串行或者并行运行
    run_serial(args) if args.P == 1 else run_parallel(args)


if __name__ == '__main__':
    main()
