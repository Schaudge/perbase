#!/usr/bin/env python3
# coding: utf8
"""
complement the long indels for hotspot drug sites from the perbase result
Written By Schaudge King.
Initiate at 2022-06-20.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals
from os import path


def convert_to_vcf(sorted_indel_list, vcf_output="check.vcf", sample_id="XXX",
                   vcf_header_path="/yunying/codes/product/module/software/calling/perbase/script/vcf.header"):
    """
    convert the parsed perbase output to standard vcf
    :param sorted_indel_list:
    :param vcf_output:
    :param vcf_header_path:
    :return:
    """
    with open(vcf_header_path) as head_in, open(vcf_output, "w") as vcf_out:
        for line in head_in:
            vcf_out.write(line.format(sample_id))
        for record in sorted_indel_list:
            chrom, position, ref, alt, del_len, ref_count, alt_count, total_depth, __ = record
            if ref_count < 1:  # for insert case
                mean_ref_count = total_depth - alt_count
                mean_alt_count = alt_count
                mean_depth = total_depth
                mean_vaf = alt_count/total_depth
            else:  # for complex deletion case
                variation_range = len(ref) if len(ref) > len(alt) else len(alt)
                mean_ref_count = round(int(ref_count)/variation_range)
                mean_alt_count = round(int(alt_count)/int(del_len))
                mean_depth = round(int(total_depth)/variation_range)
                mean_vaf = int(mean_alt_count)/mean_depth
            if ref[-1] == alt[-1]:
                pentail_trim_len = min(len(ref), len(alt))
                for reverse in range(-1, -1 * pentail_trim_len, -1):
                    if ref[reverse] != alt[reverse]:
                        break
                ref = ref[0:reverse+1]
                alt = alt[0:reverse+1]
            vcf_record_line = chrom + "\t" + str(position) + "\t.\t" + ref + "\t" + alt + "\t.\t.\tDP={}".format(mean_depth) + "\tAD:AF:DP\t" + \
                    str(mean_ref_count) + "," + str(mean_alt_count) + ":" + str(mean_vaf) + ":" + str(mean_depth)
            vcf_out.write(vcf_record_line + "\n")


def fetch_pential_long_indels(perbase_result_path, min_indels_counts=4, min_indels_length=5, prefix_context_length=10,
                              partition_coefficient=1.2) -> list:
    """
    create standard review information from perbase output file
    only for long deletion (currently)!
    :param perbase_result_path: perbase output file
    :param min_indels_counts: the minimum indel counts for positive mutation
    :param min_indels_length: the minimum length for successive indels
    :param prefix_context_length: site up/down stream seq length for indel context
    :param partition_coefficient: used for separate two adjacent indels events
    :return: list []
    """
    complement_indels_list = []
    successive_segment_location = [0, 0, "", ""]     # [chromosome, position, ref, alt]
    successive_segment_stats = [0, 0, 0, 0, 0]       # [del_length, ref_counts, alt_counts, total_depth, master_count]
    insert_segment_location = [0, 0, "", ""]
    insert_segment_stats = [0, 0, 0, 0, 0]
    insert_successive_seq = ""
    preceding_sequence = ""
    preceding_del_list = []
    if path.exists(perbase_result_path):
        opened_successive_seq, opened_insert_seq = False, False
        with open(perbase_result_path) as pinput:
            next(pinput)
            for line in pinput:
                __, chrom, pos, ref, depth, a, c, g, t, __, __, ins_seq, ins_count, dels, context_seq, master_count, __, __, near_max = line.strip().split("\t")
                position = int(pos)
                del_master_count = int(master_count)
                ref_count = int(a) if ref == "A" else (int(c) if ref == "C" else (int(g) if ref == "G" else int(t)))
                # 1. parse the insert sequence case:
                if chrom != successive_segment_location[0] or position > successive_segment_location[1] + 120:
                    preceding_sequence = ""
                    preceding_del_list = []
                    if successive_segment_stats[0] >= min_indels_length:
                        if len(successive_segment_location[2]) <= successive_segment_stats[0] + len(successive_segment_location[3]):
                            alt_trim_length = len(successive_segment_location[2]) - successive_segment_stats[0]
                            successive_segment_location[3] = successive_segment_location[3][0:alt_trim_length]
                        complement_indels_list.append(successive_segment_location + successive_segment_stats)
                    if opened_insert_seq:
                        complement_indels_list.append(insert_segment_location + insert_segment_stats)
                        opened_insert_seq = False
                    successive_segment_location = [chrom, position, "", ""]
                    successive_segment_stats = [0, 0, 0, 0, 0]
                preceding_sequence += ref
                if del_master_count > min_indels_counts:
                    if del_master_count > successive_segment_stats[4] * partition_coefficient:
                        preceding_seq_pair = (ref, ref)
                        preceding_del_size = 0
                        deletion_back_size = 0
                        for del_pos, del_count, del_depth in preceding_del_list:
                            if del_master_count * 0.98 <= del_count < del_master_count * 1.5 and \
                                    position <= del_pos + 10:
                                preceding_del_size += 1
                                if deletion_back_size < position - del_pos:
                                    deletion_back_size = position - del_pos + 1
                        under_trimmed_seq = preceding_sequence
                        context_seq_length = prefix_context_length
                        if preceding_del_size > 0:
                            under_trimmed_seq = preceding_sequence[0:len(preceding_sequence) - deletion_back_size]
                            context_seq_length = 10 - deletion_back_size + preceding_del_size
                            preceding_seq_pair = (preceding_sequence[-1 * deletion_back_size:], context_seq[context_seq_length:10])
                            position -= (deletion_back_size - 1)
                        preceding_seq_len = len(under_trimmed_seq)
                        min_compare_size = min(preceding_seq_len, context_seq_length)
                        if context_seq[context_seq_length - min_compare_size:context_seq_length] != \
                                under_trimmed_seq[preceding_seq_len - min_compare_size:]:
                            for px in range(min_compare_size, 0, -1):
                                if context_seq[context_seq_length - px] != under_trimmed_seq[preceding_seq_len - px]:
                                    preceding_seq_pair = (preceding_sequence[len(preceding_sequence)-px:], context_seq[10-px:10])
                                    position = position - px + 1
                                    break
                        successive_segment_location = [chrom, position, preceding_seq_pair[0], preceding_seq_pair[1] + context_seq[10:]]
                        successive_segment_stats = [preceding_del_size, ref_count, 0, int(depth), del_master_count]
                        opened_successive_seq = True
                    elif opened_successive_seq and len(successive_segment_location[2]) < successive_segment_stats[0] + len(successive_segment_location[3]):
                        successive_segment_location[2] += ref
                        this_del = 1 if int(dels) > successive_segment_stats[4] * 0.768 else 0
                        successive_segment_stats = [raw + add for raw, add
                                                    in zip(successive_segment_stats, [this_del, ref_count, int(dels), int(depth), 0])]
                elif opened_successive_seq and int(dels) > min_indels_counts and \
                        len(successive_segment_location[2]) < successive_segment_stats[0] + len(successive_segment_location[3]):
                    estimated_del_depth = int(dels) * (len(successive_segment_location[2]) - 1)
                    same_del_increase = 1 if estimated_del_depth >= successive_segment_stats[2] * 0.9 or (
                            int(dels) > 200 and estimated_del_depth >= successive_segment_stats[2] * 0.7) or (
                            int(dels) < 13 and estimated_del_depth >= successive_segment_stats[2] - 2) else 0
                    successive_segment_location[2] += ref
                    successive_segment_stats = [raw + add for raw, add
                                                in zip(successive_segment_stats, [same_del_increase, ref_count, int(dels), int(depth), 0])]
                elif not opened_successive_seq:
                    if 0 < successive_segment_stats[0] < min_indels_length:
                        successive_segment_location = [chrom, position, "", ""]
                        successive_segment_stats = [0, 0, 0, 0, 0]
                    if int(dels) > min_indels_counts:
                        preceding_del_list.append((position, int(dels), int(depth)))
                elif opened_successive_seq and successive_segment_stats[0] >= min_indels_length:
                    if len(successive_segment_location[2]) < successive_segment_stats[0] + len(successive_segment_location[3]):
                        successive_segment_location[2] += ref
                        successive_segment_stats = [raw + add for raw, add
                                                    in zip(successive_segment_stats, [0, ref_count, 0, int(depth), 0])]
                        # rare case for small insertion after long deletion, need padding the insertion seq
                        # TODO: too long insertion (beyond to context length) will cause wrong result
                        if successive_segment_stats[4] * 0.96 <= int(ins_count) < successive_segment_stats[4] * 1.2:
                            successive_segment_stats[0] -= len(ins_seq)
                    else:
                        opened_successive_seq = False

                # 2. parse the insert sequence case:
                if not opened_successive_seq and not opened_insert_seq and int(ins_count) > min_indels_counts and len(ins_seq) > 2:
                    insert_segment_location = [chrom, pos, ref, ref + ins_seq]
                    insert_segment_stats = [len(ins_seq), 0, int(ins_count), int(depth), int(ins_count)]
                    insert_successive_seq = ""
                    opened_insert_seq = True
                elif not opened_successive_seq and opened_insert_seq and len(insert_segment_stats[2]) < 6:
                    insert_successive_seq += ref
                    if int(ins_count) > min_indels_counts and len(ins_seq) > 2 and int(ins_count) > 1.2 * insert_segment_stats[-1]:
                        insert_segment_location = [chrom, pos, ref, ref + ins_seq]
                        insert_segment_stats = [len(ins_seq), 0, int(ins_count), int(depth), int(ins_count)]
                        insert_successive_seq = ""
                    elif int(ins_count) > min_indels_counts and len(ins_seq) > 2 and int(ins_count) > 0.9 * insert_segment_stats[-1]:
                        insert_segment_location[2] += insert_successive_seq
                        insert_segment_location[3] += (insert_successive_seq + ins_seq)
                        insert_successive_seq = ""
                elif not opened_successive_seq and opened_insert_seq and len(insert_segment_stats[2]) > 6:
                    complement_indels_list.append(insert_segment_location + insert_segment_stats)
                    opened_insert_seq = False

            if successive_segment_stats[0] >= min_indels_length:
                if len(successive_segment_location[2]) <= successive_segment_stats[0] + len(successive_segment_location[3]):
                    alt_trim_length = len(successive_segment_location[2]) - successive_segment_stats[0]
                    successive_segment_location[3] = successive_segment_location[3][0:alt_trim_length]
                complement_indels_list.append(successive_segment_location + successive_segment_stats)

    return complement_indels_list


if __name__ == '__main__':
    from sys import argv
    vcf_header_path = path.dirname(path.realpath(argv[0])) + "/vcf.header"
    if len(argv) < 2:
        print("\nBasic Usage:\n python3 convert_indels_vcf.py perbase_out_file output_vcf_file [SAMPLE_ID]\n")
        exit(1)
    elif len(argv) < 3:
        parsed_indels_list = fetch_pential_long_indels(argv[1])
        convert_to_vcf(parsed_indels_list, vcf_header_path=vcf_header_path)
    else:
        sample_id = argv[3] if len(argv) > 3 else argv[1].split("_")[0]
        parsed_indels_list = fetch_pential_long_indels(argv[1])
        convert_to_vcf(parsed_indels_list, argv[2], sample_id, vcf_header_path)

