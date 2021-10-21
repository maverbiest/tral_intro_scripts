#!/usr/bin/env python3
"""
Collect repeatlists from a directory and filter the associated repeats based on command line arguments.
The filtered repeatlists will then be serialized to a new output location.
"""

import argparse
import os

from tral.repeat.repeat import Repeat
from tral.repeat_list.repeat_list import RepeatList, pvalue

from tral_score import load_repeatlists

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--input", "-i", type=str, required=True, help="Directory containing pickled RepeatLists to score"
    )
    parser.add_argument(
        "--output", "-o", type=str, required=True, help="Directory where scored RepeatLists will be stored"
    )
    parser.add_argument(
        "--model", "-m", type=str, required=True, help="Model to use for filtering repeats. Options: phylo, phylo_gap01, phylo_gap001"
    )
    parser.add_argument(
        "--pvalue", "-p", type=float, required=True, help="P value threshold (upper bound)"
    )
    parser.add_argument(
        "--divergence", "-d", type=float, required=True, help="Divergence threshold (upper bound)"
    )
    parser.add_argument(
        "--units", "-u", type=float, required=False, help="(Optional) Number of repeat units required (lower bound, inclusive)"
    )
    parser.add_argument(
        "--unit_len", "-l", type=int, required=False, help="(Optional) Minimum unit length number (lower bound, inclusive)"
    )

    return parser.parse_args()

def filter_and_correct_tails(repeat_list, model, pval, div):
    filtered_list = []
    for repeat in repeat_list.repeats:
        if reject_repeat(repeat, model, pval, div):
            # repeat is rejected, attempt to salvage by trimming gappy leading and/or trailing unit
            leading_unit = repeat.msa[0].replace("-", "")
            trailing_unit = repeat.msa[-1].replace("-", "")     

            if len(trailing_unit) <= 0.5 * repeat.l_effective and len(repeat.msa) >= 3:
                # trailing unit is gappy, attempt trim                
                repeat = trim_rescore(repeat, model, to_trim="last")
                if not reject_repeat(repeat, model, pval, div):
                    # repeat is now significant, add and continue to next iteration to prevent further trimming
                    filtered_list.append(repeat)                
                    continue 

            if len(leading_unit) <= 0.5 * repeat.l_effective and len(repeat.msa) >= 3:
                # leading unit is gappy, attempt trim
                repeat = trim_rescore(repeat, model, to_trim="first")
                if not reject_repeat(repeat, model, pval, div):
                    # repeat is now significant, add to filtered list
                    filtered_list.append(repeat) 
            continue # continue to next iteration, otherwise rejected repeats will be appended below        
        # repeat is significant without need for trimming, add to filtered list
        filtered_list.append(repeat)
    
    return filtered_list

def reject_repeat(repeat, model, pvalue, divergence):
    """return True if repeat divergence and pvalue are below thresholds, False otherwise"""
    if repeat.d_pvalue[model] >= pvalue or repeat.d_divergence[model] >= divergence:
        return True
    return False

def trim_rescore(repeat, model, to_trim):
    if not to_trim in {"first", "last"}:
        raise ValueError("Can only accepted 'first' or 'last' position to trim")

    if to_trim == "first":
        # new_msa is old msa without fist unit, update begin value to account for this
        new_msa = repeat.msa[1:]
        new_begin = repeat.begin + len(repeat.msa[0].replace("-", ""))
    elif to_trim == "last":
        # new_msa is old msa without last unit, begin value stays the same
        new_msa = repeat.msa[0:-1]
        new_begin = repeat.begin

    # construct new, trimmed Repeat
    trimmed_repeat = Repeat(
        new_msa,
        begin=new_begin,
        sequence_type=repeat.sequence_type,
        scoreslist=[model], 
        calc_score=True, 
        calc_pvalue=True
        )
    trimmed_repeat.TRD = repeat.TRD
    trimmed_repeat.msa_original = trimmed_repeat.msa    

    return trimmed_repeat


def main():
    args = cla_parser()

    input_dir = args.input 
    output_dir = args.output
    model = args.model
    pvalue = args.pvalue
    divergence = args.divergence    
        
    for file_name, repeat_list in load_repeatlists(input_dir):        
    
        # optional: filtering for number of repeat units
        if args.units:            
            repeat_list = repeat_list.filter(
                "attribute",
                "n_effective",
                "min",
                args.units)

        # optional: filtering for unit length
        if args.unit_len:
            repeat_list = repeat_list.filter(
                "attribute",
                "l_effective",
                "min",
                args.unit_len)            

        repeat_list_filt = RepeatList(filter_and_correct_tails(repeat_list, model=model, pval=pvalue, div=divergence))

        # Clustering
        # De novo repeats are clustered for overlap (common ancestry). In case of overlap, best repeat
        # (i.e. lowest p-value and lowest divergence) is retained.
        criterion_list = [("pvalue", model), ("divergence", model)]
        repeat_list_clust = repeat_list_filt.filter("none_overlapping", ["common_ancestry"], criterion_list)

        new_file_name = file_name.replace(".pickle", "_filt.pickle")
        output_path = os.path.join(output_dir, new_file_name)
        repeat_list_clust.write(output_format="pickle", file=output_path)        

if __name__ == "__main__":
    main()
