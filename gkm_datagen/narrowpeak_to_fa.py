import argparse
import os
import numpy as np
import twobitreader
from collections import defaultdict

sizes = {
    'chr1': 249250621,
    'chr2': 243199373,
    'chr3': 198022430,
    'chr4': 191154276,
    'chr5': 180915260,
    'chr6': 171115067,
    'chr7': 159138663,
    'chrX': 155270560,
    'chr8': 146364022,
    'chr9': 141213431,
    'chr10': 135534747,
    'chr11': 135006516,
    'chr12': 133851895,
    'chr13': 115169878,
    'chr14': 107349540,
    'chr15': 102531392,
    'chr16': 90354753,
    'chr17': 81195210,
    'chr18': 78077248,
    'chr20': 63025520,
    'chrY': 59373566,
    'chr19': 59128983,
    'chr22': 51304566,
    'chr21': 48129895
}


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


def convert_coord_to_seq_letters(narrow_filename, genome_dict):
    line_count = get_line_count(filename=narrow_filename)
    coord_dict = defaultdict(list)
    seq_tuple_list = []
    with open(narrow_filename, 'r') as narrow_file:
        for line_index, line in enumerate(narrow_file):
            tokens = line.strip().split()
            chrom, start, stop = tokens[0], int(tokens[1]), int(tokens[2])
            
            dna_sequence = genome_dict[chrom].get_slice(start, stop)
            
            if len(dna_sequence) != stop - start:
                print(f'The DNA sequence is not {stop - start} bp in length. Will skip.'
                      f'chrom: {chrom} start: {start} stop: {stop}')
                continue
            
            seq_tuple_list.append((dna_sequence, chrom, start, stop))
            coord_dict[chrom].append((start, stop))
            
            if line_index % 1000 == 0:
                print(f'=> {line_index}/{line_count} = {line_index / line_count:.2%}', end='\r')
        
        print(f'\n=> Processed {line_count} lines from {narrow_filename}')
    
    return seq_tuple_list, coord_dict


def take_random_split(datapoints, second_segment_ratio):
    num_datapoints = len(datapoints)
    second_segment_indices = np.random.choice(num_datapoints, int(num_datapoints * second_segment_ratio), replace=False)
    first_segment = [datapoints[i] for i in range(num_datapoints) if i not in second_segment_indices]
    second_segment = [datapoints[i] for i in second_segment_indices]
    return first_segment, second_segment


def narrowpeak_to_fa(narrowpeak_filename, output_prefix):
    """
    Given a narrowpeak file containing the positive sequence coordinates,
    generate both positive and negative training sequences and store them into .fa.ir files
    for further processing by species_letters_from_coord.py
    :param narrowpeak_filename: narrowpeak file containing the positive sequence coordinates
    :param output_prefix: a dirctory named output_prefix will be created with everything generated saved under
    """
    
    test_ratio = 0.3
    genome_dict = twobitreader.TwoBitFile('hg19.2bit')
    
    print(f'=> Creating directory {output_prefix}')
    os.makedirs(output_prefix, exist_ok=False)
    
    positive_seq_tuple_list, positive_coord_dict = convert_coord_to_seq_letters(narrowpeak_filename, genome_dict)
    
    original_dir = os.getcwd()
    print(f'=> Changing the working directory to {output_prefix}')
    os.chdir(output_prefix)
    
    with open(f'{output_prefix}.pos.coord', 'w') as pos_coord_file:
        # coordinates are left inclusive right exclusive
        for _, chrom, start, stop in positive_seq_tuple_list:
            pos_coord_file.write(f'{chrom} {start} {stop}\n')
    
    positive_train_tuple_list, positive_test_tuple_list = take_random_split(positive_seq_tuple_list, test_ratio)
    
    print(f'positive_train_list length: {len(positive_train_tuple_list)}\n'
          f'positive_test_list length: {len(positive_test_tuple_list)}')
    
    train_positive_coord_filename = f'{output_prefix}.train.pos.coord'
    test_positive_coord_filename = f'{output_prefix}.test.pos.coord'
    # The suffix ir stands for intermediate representation
    # because the line preceding each sequence letters is of the form >chrom start stop
    # whereas that in the final fasta file will be of the form >chrom:start-stop
    # The intermediate representation is for the ease of data generation
    write_seq_intermediate_rep_and_coord_file(positive_train_tuple_list, f'{output_prefix}.train.pos.fa.ir',
                                              train_positive_coord_filename)
    write_seq_intermediate_rep_and_coord_file(positive_test_tuple_list, f'{output_prefix}.test.pos.fa.ir',
                                              test_positive_coord_filename)
    
    print('\n=> Generating negative sequences')
    negative_coord_dict = generate_negative_sequence_coord(positive_coord_dict, sizes, len(positive_seq_tuple_list),
                                                           genome_dict)
    
    negative_coord_filename = f'{output_prefix}.neg.coord'
    with open(negative_coord_filename, 'w') as neg_coord_file:
        for chrom, start_length_tuple_list in negative_coord_dict.items():
            for start_coord, length in start_length_tuple_list:
                # Writing chromosome name, start coordiante, stop coordinate
                neg_coord_file.write(f'{chrom} {start_coord} {start_coord + length}\n')
    
    neg_seq_tuple_list, _ = convert_coord_to_seq_letters(negative_coord_filename, genome_dict)
    
    neg_train_tuple_list, neg_test_tuple_list = take_random_split(neg_seq_tuple_list, test_ratio)
    print(f'negative_train_list length: {len(neg_train_tuple_list)}\n'
          f'negative_test_list length: {len(neg_test_tuple_list)}')
    
    train_neg_coord_filename = f'{output_prefix}.train.neg.coord'
    test_neg_coord_filename = f'{output_prefix}.test.neg.coord'
    
    write_seq_intermediate_rep_and_coord_file(neg_train_tuple_list, f'{output_prefix}.train.neg.fa.ir',
                                              train_neg_coord_filename)
    write_seq_intermediate_rep_and_coord_file(neg_test_tuple_list, f'{output_prefix}.test.neg.fa.ir',
                                              test_neg_coord_filename)
    
    print(f'=> Changing back to the original directory {original_dir}')
    os.chdir(original_dir)


def write_seq_intermediate_rep_and_coord_file(seq_tuple_list, seq_ir_filename, coord_filename):
    """
    Writing to two files, saving the intermediate representation (IR) of the sequences and the coordinates of those sequences.
    The IR differs from the regular fasta file only in that the description line of each sequence
    IR will be >chrom start stop
    while the regular fasta file will be >chrom:start-stop
    :param seq_tuple_list: a list of tuples of the form (seq, chrom, start, stop)
    :param seq_ir_filename: filename for the intermediate representation fa.ir file e.g. output_prefix.train.neg.fa.ir
    :param coord_filename: filename of the file to save the coordinates in
    """
    with open(seq_ir_filename, 'w') as seq_ir_file, open(coord_filename, 'w') as coord_file:
        print(f'=> Writing to {seq_ir_file.name} and {coord_file.name}')
        for test_seq, chrom, start, stop in seq_tuple_list:
            seq_ir_file.write(f'>{chrom} {start} {stop}\n{test_seq}\n')
            coord_file.write(f'{chrom} {start} {stop}\n')


def generate_negative_sequence_coord(positive_coord_dict, chrom_sizes, num_samples_required, genome_dict):
    """
    :param positive_coord_dict: mapping chromosome name to a list of tuple of (start, end) of positive sequences, so that
    we do not sample from those coordiantes.
    :param chrom_sizes: mapping chromosome name to their sizes
    :param num_samples_required: total number of negative samples required.
    :param genome_dict: use twobitreader on hg19.2bit
    """
    total_num_coordinates = sum(chrom_sizes.values())
    print(f'-> Total number of coordinates: {total_num_coordinates}')
    sample_coord = {}
    
    print(positive_coord_dict)
    seq_lengths = [stop - start for tuple_list in positive_coord_dict.values() for start, stop in tuple_list]
    max_len = max(seq_lengths)
    print(f'-> Longest positive sequence length encountered: {max_len}')
    
    # Subtract the maximum possible sequence length max_len from each of the start coordinates
    # The negative samples will have length equal to at most max_len so that adding the actual
    # length of the negative samples to the start coordinate of the sample will never overlap
    # any of the positive sequences.
    for chrom, positive_seq_start_stop_list in positive_coord_dict.items():
        forbidden_coord_list = []
        for start, stop in positive_seq_start_stop_list:
            forbidden_coord_list.append((start - max_len, stop))
        forbidden_coord_list.sort(key=lambda start_stop_tuple: start_stop_tuple[0])
        
        upper_bound = chrom_sizes[chrom]
        print(f'-> Chromosome {chrom} has max coordinate: {upper_bound}')
        
        for start, stop in forbidden_coord_list:
            assert start < stop
            upper_bound -= stop - start
        print(f'-> Upper bound changed to: {upper_bound} in the mapping range')
        
        num_samples = int(num_samples_required * chrom_sizes[chrom] / total_num_coordinates)
        # sample twice the amount and get rid of indices that are too close together
        # Then resample one more time to make the length become num_samples
        sampled_indices = list(np.random.choice(upper_bound, num_samples * 2, replace=False))
        length_sampling_indices = np.random.choice(len(seq_lengths), len(sampled_indices), replace=True)
        start_length_tuple_list = list(zip(sampled_indices, [seq_lengths[i] for i in length_sampling_indices]))
        
        print('=> Filtering coordinates that are too close')
        filter_close_coordinates(start_length_tuple_list, max_len)
        
        print('=> Mapping sampled numbers back to the real coordinates')
        for index in range(len(start_length_tuple_list)):
            start, length = start_length_tuple_list[index]
            for forbidden_start, forbidden_stop in forbidden_coord_list:
                if start >= forbidden_start:
                    start += forbidden_stop - forbidden_start
                    start_length_tuple_list[index] = (start, length)
                else:
                    break
        
        # Remove samples containing N
        current_index = len(start_length_tuple_list) - 1
        chrom_reader = genome_dict[chrom]
        num_uncertain_samples = 0
        while current_index >= 0:
            start, length = start_length_tuple_list[current_index]
            if 'N' in chrom_reader.get_slice(start, start + length):
                start_length_tuple_list.pop(current_index)
                num_uncertain_samples += 1
            current_index -= 1
        
        print(f'-> Removed {num_uncertain_samples} potential sequences containing N')
        
        if len(start_length_tuple_list) < num_samples:
            raise ValueError('not enough number of indices after filtering out indices that are close together')
        
        start_length_tuple_list = downsample(start_length_tuple_list, num_samples)
        sample_coord[chrom] = start_length_tuple_list
        
        print('=> Mapping complete\n')
    
    return sample_coord


def filter_close_coordinates(start_length_tuple_list, min_distance):
    """
    Filters the sample_indices in place. The first pivot will be the last element of the list. At each iteration,
    if the value of the element to the left of the pivot is less than min_distance from the pivot value,
    it will be removed. Otherwise, the pivot moves to the left by one.
    The returned list will be sorted in increasing order and contain elements that are at least
    min_distance away from each other.
    :param sample_indices: the list of coordinates to be filtered in place.
    :param min_distance: the minimum distance between any two elements in the sample_indices.
    :return: a list containing a subset of the sample_indices, sorted in increasing order and containing
    elements which are at least min_distance away from each other.
    """
    start_length_tuple_list.sort(key=lambda pair: pair[0])
    current_value, current_length = start_length_tuple_list[-1]
    current_index = len(start_length_tuple_list) - 2
    
    while current_index >= 0:
        if current_value - start_length_tuple_list[current_index][0] < min_distance:
            start_length_tuple_list.pop(current_index)
            current_index -= 1
        else:
            current_value, current_length = start_length_tuple_list[current_index]
            current_index -= 1


def downsample(samples, target_count):
    if len(samples) < target_count:
        raise ValueError(
            f'The number of elements in the samples {len(samples)} is less than the target_count {target_count}')
    target_indices = np.random.choice(len(samples), target_count, replace=False)
    return [samples[i] for i in target_indices]


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('filename')
    arg_parser.add_argument('output_prefix')
    args = arg_parser.parse_args()
    # example output_prefix uw_gm12878_ctcf
    narrowpeak_to_fa(args.filename, args.output_prefix)
