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
            chrom, position, ref, alt, del_len, ref_count, alt_count, total_depth = record
            print(record)
            variation_range = len(ref) if len(ref) > len(alt) else len(alt)
            mean_ref_count = round(int(ref_count)/variation_range)
            mean_alt_count = round(int(alt_count)/int(del_len))
            mean_depth = round(int(total_depth)/variation_range)
            mean_vaf = int(mean_alt_count)/mean_depth
            pentail_trim_len = min(len(ref), len(alt))
            for reverse in range(-1, -1 * pentail_trim_len, -1):
                if ref[reverse] != alt[reverse]:
                    break
            trimed_ref = ref[0:reverse+1]
            trimed_alt = alt[0:reverse+1]
            vcf_record_line = chrom + "\t" + str(position) + "\t.\t" + trimed_ref + "\t" + trimed_alt + "\t.\t.\tDP={}".format(mean_depth) + "\tAD:AF:DP\t" + \
                    str(mean_ref_count) + "," + str(mean_alt_count) + ":" + str(mean_vaf) + ":" + str(mean_depth)
            vcf_out.write(vcf_record_line + "\n")


def fetch_pential_long_indels(perbase_result_path, min_indels_counts=5, min_indels_length=9) -> list:
    """
    create standard review information from perbase output file
    only for long deletion (currently)!
    :param perbase_result_path: perbase output file
    :param min_indels_counts: the minimum indel counts for postive mutation
    :param min_indels_length: the minimum length for successive indels
    :return: list []
    """
    complement_indels_list = []
    successive_segment_location = [0, 0, "", ""]     # [chromosome, start, ref, alt]
    successive_segment_stats = [0, 0, 0, 0]          # [del_length, ref_counts, alt_counts, total_depth]
    if path.exists(perbase_result_path):
        with open(perbase_result_path) as pinput:
            next(pinput)
            for line in pinput:
                __, chrom, start, ref, depth, a, c, g, t, __, ins, __, __, dels, successive_seq, master_count, __, __, near_max = line.strip().split("\t")
                ref_count = int(a) if ref == "A" else (int(c) if ref == "C" else (int(g) if ref == "G" else int(t)))
                if chrom != successive_segment_location[0] or int(start) > successive_segment_location[1] + 500:
                    if successive_segment_stats[0] >= min_indels_length:
                        complement_indels_list.append(successive_segment_location + successive_segment_stats)
                    successive_segment_location = [0, 0, "", ""]
                    successive_segment_stats = [0, 0, 0, 0]
                if int(master_count) > min_indels_counts:
                    if successive_segment_stats[2] < 1 or int(master_count) > successive_segment_stats[2]:
                        successive_segment_location = [chrom, int(start), ref, ref + successive_seq]
                        successive_segment_stats = [0, ref_count, 0, int(depth)]
                    else:
                        successive_segment_location[2] += ref
                        successive_segment_stats = [raw + add for raw, add
                                                    in zip(successive_segment_stats, [1, ref_count, int(dels), int(depth)])]
                elif int(dels) > min_indels_counts:
                    successive_segment_location[2] += ref
                    successive_segment_stats = [raw + add for raw, add
                                                in zip(successive_segment_stats, [1, ref_count, int(dels), int(depth)])]
                elif 0 < successive_segment_stats[0] < min_indels_length:
                    successive_segment_location = [0, 0, "", ""]
                    successive_segment_stats = [0, 0, 0, 0]
                elif successive_segment_stats[0] >= min_indels_counts:
                    if len(successive_segment_location[2]) < successive_segment_stats[0] + \
                            len(successive_segment_location[3]):
                        successive_segment_location[2] += ref
                        successive_segment_stats = [raw + add for raw, add
                                                    in zip(successive_segment_stats, [0, ref_count, 0, int(depth)])]
            if successive_segment_stats[0] >= min_indels_length:
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

