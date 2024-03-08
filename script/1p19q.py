# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 14:04:42 2021
@author: xijun
"""
import sys
import numpy as np


def bed2pos(rs_pos_bed: str):
    index_dict = {}
    with open(rs_pos_bed) as pos_bed_in:
        for line in pos_bed_in:
            chrom, start, pos, ref = line.strip().split("\t")
            index_dict[chrom + ":" + pos] = [chrom, pos, '0', ref, ref, '0', '0', 'NA', '2']
    return index_dict


def generate_writeline(perbase_filename: str, rs_pos_dict: dict):
    out_list = [['CHR', 'POS', 'DEPTH', 'REF', 'ALT', 'REF_DEPTH', 'ALT_DEPTH', 'FREQ', 'GENOTYPE']]
    with open(perbase_filename) as depth_in:
        for line in depth_in:
            __, CHROM, POS, REF_BASE, DEPTH, A, C, G, T, N = line.strip().split("\t")[:10]
            current_idx = CHROM + ":" + POS
            while current_idx in rs_pos_dict:
                ordered_idx = next(iter(rs_pos_dict))
                if ordered_idx != current_idx:
                    out_list.append(rs_pos_dict.pop(ordered_idx))
                else:
                    rs_pos_dict.pop(ordered_idx)
                    current_idx = ":"   # set to null key for break while loop!
                    POS, DEPTH, A, C, G, T, N = [int(i) for i in [POS, DEPTH, A, C, G, T, N]]
                    base_depth_list = [A, C, G, T, N]
                    base_name_list = ['A', 'C', 'G', 'T', 'N']
                    ref_depth = base_depth_list[base_name_list.index(REF_BASE)]

                    max_depth = max(base_depth_list)
                    max_depth_base = base_name_list[np.argmax(base_depth_list)]

                    max2_depth = sorted(base_depth_list)[-2]
                    max2_depth_base = base_name_list[base_depth_list.index(max2_depth)]

                    # find alt and alt_depth
                    if REF_BASE == max_depth_base:
                        alt = max2_depth_base
                        alt_depth = max2_depth
                        if max2_depth == max_depth:
                            alt = 'N'
                            if max_depth_base == REF_BASE:
                                alt = max2_depth_base
                            elif max2_depth_base == REF_BASE:
                                alt = max_depth_base
                    else:  # germline
                        alt = max_depth_base
                        alt_depth = max_depth

                    # calculate freq
                    freq = (alt_depth/DEPTH)
                    if DEPTH == 0:
                        freq = '0'
                    # judge genotype
                    if freq < 0.05 or freq > 0.95:
                        genotype = '2'
                    elif 0.40 < freq < 0.60:
                        genotype = '0'
                    else:
                        genotype = '1'
                    out_list.append([CHROM, POS, DEPTH, REF_BASE, alt, ref_depth, alt_depth, freq, genotype])

    while len(rs_pos_dict):
        ordered_idx = next(iter(rs_pos_dict))
        out_list.append(rs_pos_dict.pop(ordered_idx))
    return out_list


def write_list_out(sta_list, out_file: str):
    with open(out_file, 'w') as f_out:
        for line_list in sta_list:
            out_line = '\t'.join([str(i) for i in line_list])
            f_out.write(out_line + '\n')


if __name__ == '__main__':
    rs_pos_bed = sys.argv[1]  # '1p19q.rs.bed'
    perbase_file = sys.argv[2]  # 'out'
    out_file = sys.argv[3]  # 'N_1p19q.result.xls'
    rs_pos_dict = bed2pos(rs_pos_bed)
    out_list = generate_writeline(perbase_file, rs_pos_dict)
    write_list_out(out_list, out_file)

