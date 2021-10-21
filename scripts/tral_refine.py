#!/usr/bin/env python3
import sys
import os

import argparse
import time
# import shutil

from Bio import SeqIO
import numpy as np
from tral.sequence.sequence import Sequence
from tral.repeat.repeat import Repeat
import tral.repeat_list.repeat_list as repeat_list
import tral.repeat.repeat_align as repeat_realign
from tral.hmm.hmm import HMM
from tral.hmm import hmm_viterbi

from tral_score import load_repeatlists
from tral_filter import filter_and_correct_tails

def cla_parser():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--repeat_dir", "-r", type=str, required=True, help="Directory where pickled RepeatLists are stored"
    )    
    parser.add_argument(
        "--output_dir", "-o", type=str, required=False, help="Refined TRs will be deposited here"
    )
    parser.add_argument(
        "--fasta_file", "-f", type=str, required=True, help="Fasta file containing the sequence that the TRs originate from. Fasta file is expected to contain multiple sequences, out of which the one matching the repeatlist will be selected"
    )
    parser.add_argument(
        "-s", "--seq_type", type=str, required=True, help="Type of sequence that is being analyzed, currently supports 'AA' and 'DNA'"
    )
    parser.add_argument(
        "--score_type", "-c", type=str, required=True, help="Type of score to use options: phylo, phylo_gap01, phylo_gap001"
    )
    parser.add_argument(
        "--pvalue", "-p", type=float, default=0.05, help="P value threshold for accepting HMM refined repeats (upper bound)"
    )
    parser.add_argument(
        "--divergence", "-d", type=float, default=0.1, help="Divergence threshold for accepting HMM refined repeats (upper bound)"
    )
    parser.add_argument(
        "--max_seqs", "-m", type=int, default=False, help="Maximum number of sequences from the fasta file that will be analyzed"
    )    
    

    return parser.parse_args()

def refine_repeatlist(seq, tag, score_type, seq_type, pval, div):
    """ Take a sequence with one or more repeat lists associated to it, and a tag specifying repeat list of interest
    (one sequence can have multiple repeat lists associated to it, each identifiable with a tag). Then, create a cpHMM
    for each tandem repeat, detect TRs in sequence again and see if this improves (refines) the de novo TR
    seq (tral.sequence.Sequence):
                    A sequence with one or more RepeatLists associated to it    
    tag (str):      Which RepeatList associated to seq should be used?
    score_type (str):    
                    Score used to evaluate TRs
    seq_type (str): Are we analyzing DNA or protein sequence?
    Returns
    refined_list (list):
                    A list containing HMM refined tandem repeats
    """
    if not score_type in {"phylo", "phylo_gap01", "phylo_gap001"}:
        raise ValueError("Options for score_type: phylo, phylo_gap01, phylo_gap001")
    if not seq_type in {"DNA", "AA"}:
        raise ValueError("Options for seq_type: DNA, AA")

    refined_list = repeat_list.RepeatList(repeats=[])
    skipped_count, realigned_count, use_refined_count = 0, 0, 0

    print(f"Length of sequence: {len(seq.seq)}")
    print(f"Total number of repeats in RepeatList: {len(seq.get_repeatlist(tag).repeats)}")
    for counter, repeat in enumerate(seq.get_repeatlist(tag).repeats):  
        denovo_hmm = HMM.create(input_format='repeat', repeat=repeat)    

        repeat_begin, unaligned_msa = get_unaligned_msa(repeat, denovo_hmm, seq)
        if not unaligned_msa:
            refined_list.repeats.append(repeat)
            skipped_count += 1
            continue        
        # continue

        if new_units_same_as_old(repeat, unaligned_msa):               
            refined_list.repeats.append(repeat)
            skipped_count += 1
            continue

        skip_realign, filled_msa = skip_realign_new_check(unaligned_msa)        
        if not skip_realign:
            realigned_msa = repeat_realign.realign_repeat(filled_msa,
                                                        "mafft",
                                                        seq_type)
            hmm_repeat = Repeat(begin=repeat_begin, msa=realigned_msa)
            realigned_count += 1            
        else:            
            hmm_repeat = Repeat(begin=repeat_begin, msa=filled_msa)
            skipped_count += 1  
        
        if hasattr(repeat, "TRD"):
            hmm_repeat.TRD = repeat.TRD
        else:
            hmm_repeat.TRD = None
        hmm_repeat.model = "cpHMM"
        hmm_repeat.calculate_scores(scoreslist=[score_type])
        hmm_repeat.calculate_pvalues(scoreslist=[score_type])
        seq.repeat_in_sequence(hmm_repeat)
        
        try:
            hmm_repeat_filtered = filter_and_correct_tails(repeat_list.RepeatList([hmm_repeat]), model=score_type, pval=pval, div=div)[0]
            if repeat_list.two_repeats_overlap(
                    "shared_char",
                    repeat,
                    hmm_repeat_filtered):
                refined_list.repeats.append(hmm_repeat_filtered)
                use_refined_count += 1
                # continue to next iteration otherwise old repeat will be added to refined list below
                continue 
        except IndexError:
            # HMM repeat was rejected based on pvalue or divergence threshold  
            pass
        refined_list.repeats.append(repeat)
                  
    print(f"skipped: {skipped_count}, realigned: {realigned_count}, used refined: {use_refined_count}")

    # Clustering
    # HMM repeats are clustered for overlap (common ancestry). In case of overlap, best repeat
    # (i.e. lowest p-value and lowest divergence) is retained.
    criterion_list = [("pvalue", score_type), ("divergence", score_type)]
    refined_list = refined_list.filter("none_overlapping", ["common_ancestry"], criterion_list)

    # print(f"Number of repeats after refinement and reclustering: {len(refined_list.repeats)}", end="\r")

    return refined_list

def get_unaligned_msa(repeat, denovo_hmm, seq, units_before=10, units_after=10):
    update_window = True
    while update_window:                     
        slice_index, seq_slice = get_subsequence(sequence=seq, repeat=repeat, before=units_before, after=units_after)                    
        most_likely_path = denovo_hmm.viterbi(seq_slice.seq)        
        if not most_likely_path:
            return None, None
        slice_match_index = path_match_indices(most_likely_path)
        
        if len(seq_slice.seq) == len(seq.seq):
            break

        update_window = False   
        if slice_match_index[0] == 0 and not slice_index[0] == 0:
            units_before += 10
            update_window = True
        if slice_match_index[1] == len(seq_slice.seq) - 1 and not slice_index[1] == len(seq.seq) - 1:
            units_after += 10
            update_window = True
      
    repeat_begin = slice_index[0] + slice_match_index[0] + 1                 

    unaligned_msa = hmm_viterbi.hmm_path_to_non_aligned_tandem_repeat_units(
                seq_slice.seq,
                most_likely_path,
                denovo_hmm.l_effective)
    return repeat_begin, unaligned_msa

def skip_realign_current_check(msa):
    length_tuple = (len(unit) for unit in msa)
    msa = fill_out_msa_units(msa, max(length_tuple))
    
    if len(msa) <= 2 and len(msa[0]) <= 10:
        return True, msa
    return False, msa

def skip_realign_new_check(msa):
    length_tuple = tuple((len(unit) for unit in msa))
    msa = fill_out_msa_units(msa, max(length_tuple)) 

    # are all units of length 1?
    if len(msa[0]) == 1:
        return True, msa

    # existing check in TRAL: is repeat made up of 2 units and are they length 10 or shorter?
    if len(msa) <= 2 and len(msa[0]) <= 10:
        return True, msa    

    # is the naive alignment already optimal?
    if check_if_homogeneous(msa):
        return True, msa

    return False, msa

def fill_out_msa_units(repeat_units, max_unit_length):
    pre_processed_units = []
    for unit in repeat_units:
        filled_unit = f"{unit}{'-' * (max_unit_length - len(unit))}"
        pre_processed_units.append(filled_unit)
    return pre_processed_units

def new_units_same_as_old(unrefined_repeat, hmm_msa):
    # in order for the new and old msa to match, they have to contain an equal number of units
    if not len(unrefined_repeat.msa) == len(hmm_msa):
        return False
    
    # do the nucleotides/amino acids in each unit (excluding gaps) of the old msa match the new msa? 
    for i, unit in enumerate(unrefined_repeat.msa):
        if not unit.replace("-", "") == hmm_msa[i]:
            return False

    # all units match, therefore it is safe to use the old repeat and skip further refinement/realignment
    return True

def check_if_homogeneous(msa):
    # For each column in the transposed MSA, transform the list of characters into a set of unique
    # characters (excluding gaps) in that column in the following steps: 
    # ['A', 'A', 'A', '-'] -> {'A', '-'} -> {'A'}
    msa_columns = np.array([list(unit) for unit in msa], dtype=object).transpose()
    
    unique_chars_per_col = tuple(set(col) for col in msa_columns)
    for i in unique_chars_per_col:
        i.discard("-")

    # If all sets contain only one character it means all columns in the MSA are homogeneous -> return True and skip realignment
    return all(len(col) == 1 for col in unique_chars_per_col)

def get_subsequence(sequence, repeat, before, after):
    begin_index = max(
        repeat.begin - 1 - (before * repeat.l_effective),
        0)
    end_index = min(
        repeat.begin - 1 + repeat.repeat_region_length + (after * repeat.l_effective),
        len(sequence.seq) - 1)

    return (begin_index, end_index), Sequence(sequence.seq[begin_index : end_index + 1], sequence_type=sequence.sequence_type)

def path_match_indices(path):
    begin = None
    for i, char in enumerate(path):
        if begin is None and char != "N":
            begin = i
        if char == "C":
            return (begin, i)
    return (begin, i)


def main():
    args = cla_parser()

    seq_dict = dict()
    seq_parser = SeqIO.parse(args.fasta_file, "fasta")
    for record in seq_parser:       
        try: 
            seq_dict[record.id] = record
        except KeyError:
            raise Exception(f"During making of sequence for: {record.id}")

    finished_sequences = {i.replace("_refined.pickle", "") for i in os.listdir(args.output_dir) if i.endswith("_refined.pickle")}

    seq_counter = 0
    for file_name, repeatlist in load_repeatlists(args.repeat_dir):                        
        seq_of_origin = file_name.replace(".pickle", "")
        seq_of_origin = seq_of_origin.replace("_filt", "") # may or may not be present in file name
        if seq_of_origin in finished_sequences:
            continue

        print(f"Starting HMM refinement of RepeatList originating from {seq_of_origin}")
        seq_counter += 1

        try:
            current_seq = seq_dict[seq_of_origin]
            current_seq = Sequence(
                seq=str(current_seq.seq),
                name=current_seq.id,
                sequence_type=args.seq_type
            )
        except KeyError:
            print(f"WARNING!!! No sequence was found in fasta file for RepeatList originating from {seq_of_origin}, skipping analysis!\n")
            continue

        current_seq.set_repeatlist(repeatlist, "unrefined")
        start = time.time()
        refined_list = refine_repeatlist(
            seq=current_seq, 
            tag="unrefined",
            score_type=args.score_type, 
            seq_type=args.seq_type,
            pval=args.pvalue,
            div=args.divergence
            )  
        print(f"Refining Repeats took {time.time() - start} seconds\n")     

        output_file_name = os.path.join(args.output_dir, f"{current_seq.name}_refined.pickle")
        refined_list.write(output_format="pickle", file=output_file_name)

        # local_scratch = os.environ['LSFM_CLUSTER_LOCAL_SCRATCH_JOB_PATH']
        # for item in os.listdir(local_scratch):
        #     if item.startswith("tmp"):
        #         shutil.rmtree(os.path.join(local_scratch, item))

        if args.max_seqs:
            if seq_counter == args.max_seqs:
                break

if __name__ == "__main__":
    main()
