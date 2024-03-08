# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 10:05:18 2021

@author: sunxj
"""
import sys
from collections import defaultdict


def generate_ref_dict(fasta):
    rf_dict = {}
    with open(fasta, 'r')as f_in:
        data_list = f_in.readlines()
    for line in data_list:
        if line.startswith('>'):
            seq_id = line.strip('>|\n')
            rf_dict[seq_id] = ''
        else:
            rf_dict[seq_id] += line.strip('\n')
    return rf_dict


def deal_with_wig(perbase_file, rf_dict, out_file):
    sta_list = [['seqName', 'pos', 'ref', 'A', 'C', 'G', 'T', 'N', 'Cpg', 'NonCpg']]
    with open(perbase_file, 'r') as f_in:
        next(f_in)
        for line in f_in:
            __, CHROM, POS, REF_BASE, DEPTH, A, C, G, T, N = line.strip().split("\t")[:10]
            POS, DEPTH, A, C, G, T, N = [int(i) for i in [POS, DEPTH, A, C, G, T, N]]
            ref_seq = rf_dict[CHROM]
            Cpg = 0
            NoCpg = 0
            if ref_seq[POS] == 'G' and ref_seq[POS-1] == 'C':
                if C + T != 0:
                    Cpg = C / (C + T)
            elif (POS+2) > len(ref_seq):
                if C + T != 0:
                    NoCpg = C / (C + T)
            else:
                if C + T != 0:
                    NoCpg = C / (C + T)
            write_line = [CHROM, POS, REF_BASE, A, C, G, T, N, Cpg, NoCpg]
            sta_list.append(write_line)
        write_list_out(sta_list, out_file)
    return sta_list


def final_result(sta_list, out_file):
    title_list = ['seqName', 'Total C to T conversions in CpG context', 'Total methylated C in CpG context', 'C methylated in CpG context']
    dat_dict = defaultdict(dict)
    for line in sta_list[1:]:
        CHROM, POS, REF_BASE, A, C, G, T, N, Cpg, NoCpg = line
        dat_dict[CHROM][POS] = {}
        dat_dict[CHROM][POS]['ref'] = REF_BASE
        dat_dict[CHROM][POS]['methy'] = C
        dat_dict[CHROM][POS]['depth'] = T
    for seqname in dat_dict:
        seqname_methy_reads_on_cpg = 0
        seqname_all_reads_on_cpg = 0
        for pos in dat_dict[seqname]:
            if pos-1 not in dat_dict[seqname]:
                continue
            if dat_dict[seqname][pos]['ref'] == 'G':
                if dat_dict[seqname][pos-1]['ref'] == 'C':
                    print(pos)
                    seqname_all_reads_on_cpg += dat_dict[seqname][pos-1]['depth']
                    seqname_methy_reads_on_cpg += dat_dict[seqname][pos-1]['methy']
    if len(sta_list) < 2:
        write_line = ["NA", "NA", "NA", "NA"]
    else:
        write_line = [seqname, seqname_all_reads_on_cpg, seqname_methy_reads_on_cpg,
                      float(seqname_methy_reads_on_cpg / (seqname_all_reads_on_cpg + seqname_methy_reads_on_cpg))]
    write_list_out([title_list, write_line], out_file)


def write_list_out(sta_list, out_file):
    with open(out_file, 'w') as f_out:
        for line_list in sta_list:
            out_line = '\t'.join([str(i) for i in line_list])
            f_out.write(out_line + '\n')


if __name__ == '__main__':
    fasta_file = sys.argv[1]  # 'mgmt.fa'
    perbase_file = sys.argv[2]  # 'gzy210310186PCR_g_M.result.xls.wig'
    out_file_1 = sys.argv[3]  # 'result.xls'
    out_file_2 = sys.argv[4]  # 'result.methy.xls'

    rf_dict = generate_ref_dict(fasta_file)
    sta_list = deal_with_wig(perbase_file, rf_dict, out_file_1)
    final_result(sta_list, out_file_2)

