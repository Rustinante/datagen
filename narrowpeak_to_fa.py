import argparse
import numpy as np
import twobitreader
from collections import defaultdict

dirname = 'uw_gm12878_ctcf'

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
    'chr19': 59128983,
    'chr22': 51304566,
    'chr21': 48129895
}

padding_each_side = 0
target_seq_length = 250


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


def convert_coord_to_seq_letters(narrow_filename, genome_dict):
    line_count = get_line_count(filename=narrow_filename)
    coord_dict = defaultdict(list)
    seq_list = []
    with open(narrow_filename, 'r') as narrow_file:
        for line_index, line in enumerate(narrow_file):
            tokens = line.strip().split()
            chrom, start, stop = tokens[0], int(tokens[1]), int(tokens[2])
            seq_length = stop - start
            
            if seq_length == target_seq_length:
                dna_sequence = genome_dict[chrom].get_slice(start, stop)
                
                if len(dna_sequence) != target_seq_length:
                    print(f'The DNA sequence is not {target_seq_length} bp in length. Will skip.'
                          f'chrom: {chrom} start: {start} stop: {stop}')
                    continue
                
                seq_list.append((dna_sequence, chrom, start, stop))
                coord_dict[chrom].append((start, stop))
            
            if line_index % 1000 == 0:
                print('=> {}/{} = {:.2f}%'.format(line_index, line_count, line_index / line_count * 100), end='\r')
        
        print(f'\n=> Processed {line_count} lines from {narrow_filename}')
    
    return seq_list, coord_dict


def take_random_split(datapoints, second_segment_ratio):
    num_datapoints = len(datapoints)
    second_segment_indices = np.random.choice(num_datapoints, int(num_datapoints * second_segment_ratio), replace=False)
    first_segment = [datapoints[i] for i in range(num_datapoints) if i not in second_segment_indices]
    second_segment = [datapoints[i] for i in second_segment_indices]
    return first_segment, second_segment


def narrowpeak_to_fa(filename):
    test_ratio = 0.2
    genome_dict = twobitreader.TwoBitFile('hg19.2bit')
    
    seq_list, positive_coord_dict = convert_coord_to_seq_letters(filename, genome_dict)
    
    with open('uw_gm12878_ctcf.pos.coord', 'w') as pos_coord_file:
        # coordinates are left inclusive right exclusive
        for _, chrom, start, stop in seq_list:
            pos_coord_file.write('{} {} {}\n'.format(chrom, start, stop))
    
    positive_train_list, positive_test_list = take_random_split(seq_list, test_ratio)
    
    print(f'positive_train_list length: {len(positive_train_list)}\n'
          f'positive_test_list length: {len(positive_test_list)}')
    
    with open('uw_gm12878_ctcf.train.pos.fa', 'w') as train_posfile, \
            open('uw_gm12878_ctcf.test.pos.fa', 'w') as test_posfile:
        
        print(f'=> Writing to uw_gm12878_ctcf.train.pos.fa')
        for train_seq, chrom, start, stop in positive_train_list:
            train_posfile.write(f'>{chrom} {start}\n{train_seq}\n')
        
        print('=> Writing to uw_gm12878_ctcf.test.pos.fa')
        for test_seq, chrom, start, stop in positive_test_list:
            test_posfile.write(f'>{chrom} {start}\n{test_seq}\n')
    
    print('\n=> Generating negative sequences')
    negative_coord_dict = generate_negative_sequence_coord(positive_coord_dict, sizes, len(seq_list))

    negative_coord_filename = 'uw_gm12878_ctcf.neg.coord'
    with open(negative_coord_filename, 'w') as neg_coord_file:
        for chrom, start_coord_list in negative_coord_dict.items():
            for start_coord in start_coord_list:
                # Writing chromosome name, start coordiante, stop coordinate
                neg_coord_file.write(f'{chrom}, {start_coord}, {start_coord + target_seq_length}\n')
                
    neg_seq_list, _ = convert_coord_to_seq_letters(negative_coord_filename, genome_dict)
    
    neg_train_list, neg_test_list = take_random_split(neg_seq_list, test_ratio)
    print(f'negative_train_list length: {len(neg_train_list)}\n'
          f'negative_test_list length: {len(neg_test_list)}')
    
    neg_train_filename = 'uw_gm12878_ctcf.train.neg.fa'
    neg_test_filename = 'uw_gm12878_ctcf.test.neg.fa'
    with open(neg_train_filename, 'w') as train_negfile, \
            open(neg_test_filename, 'w') as test_negfile:
        
        print(f'=> Writing to {neg_train_filename}')
        for train_seq, chrom, start, stop in neg_train_list:
            train_negfile.write(f'>{chrom} {start}\n{train_seq}\n')
            
        print(f'=> Writing to {neg_test_filename}')
        for test_seq, chrom, start, stop in neg_test_list:
            test_negfile.write(f'>{chrom} {start}\n{test_seq}\n')
        

def generate_negative_sequence_coord(positive_coord_dict, chrom_sizes, num_samples_required):
    """
    :param positive_coord_dict: mapping chromosome name to a list of tuple of (start, end) of positive sequences, so that
    we do not samplee efrom those coordiantes.
    :param chrom_sizes: mapping chromosome name to their sizes
    :param num_samples_required: total number of negative samples required.
    """
    total_num_coordinates = sum(chrom_sizes.values())
    print(f'-> Total number of coordinates: {total_num_coordinates}')
    sample_coord = {}
    
    for chrom, positive_seq_start_stop_list in positive_coord_dict.items():
        forbidden_coord_list = []
        for start, stop in positive_seq_start_stop_list:
            forbidden_coord_list.append((start - target_seq_length, stop))
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
        
        print('=> Filtering coordinates that are too close')
        filter_close_coordinates(sampled_indices, target_seq_length)
        sampled_indices = downsample(sampled_indices, num_samples)
        
        if len(sampled_indices) < num_samples:
            raise ValueError('not enough number of indices after filtering out indices that are close together')
        
        print('=> Mapping sampled numbers back to the real coordinates')
        for index, start in enumerate(sampled_indices):
            for forbidden_start, forbidden_stop in forbidden_coord_list:
                if start >= forbidden_start:
                    increment = forbidden_stop - forbidden_start
                    sampled_indices[index] += increment
                    start += increment
                else:
                    break
        
        sample_coord[chrom] = sampled_indices
    
        print('=> Mapping complete\n')
    
    return sample_coord


def filter_close_coordinates(sample_indices, min_distance):
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
    sample_indices.sort()
    current_value = sample_indices[-1]
    current_index = len(sample_indices) - 2
    
    while current_index >= 0:
        if current_value - sample_indices[current_index] < min_distance:
            sample_indices.pop(current_index)
            current_index -= 1
        else:
            current_value = sample_indices[current_index]
            current_index -= 1


def downsample(samples, target_count):
    target_indices = np.random.choice(len(samples), target_count, replace=False)
    return [samples[i] for i in target_indices]


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('filename')
    args = arg_parser.parse_args()
    narrowpeak_to_fa(args.filename)
