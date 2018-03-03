import argparse
import numpy as np
import twobitreader
from collections import defaultdict
import os
from matplotlib import pyplot as plt


dirname = 'uw_gm12878_ctcf'

def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


def narrowpead_to_fa(filename):
    genome_dict = twobitreader.TwoBitFile('hg19.2bit')
    padding_each_side = 50
    
    line_count = get_line_count(filename=filename)
    
    # if not os.path.isdir(dirname):
    #     print('\n=> Creatingn directory {}'.format(dirname))
    #     os.makedirs(dirname)
    
    coord_dict = defaultdict(list)
    
    with open(filename, 'r') as narrow_file, open('uw_gm12878_ctcf.pos.coord', 'w') as pos_coord_file, \
            open('uw_gm12878_ctcf.train.pos.fa', 'w') as train_posfile, open('uw_gm12878_ctcf.test.pos.fa', 'w') as test_posfile:
        seq_list = []
        header_line_list = []
        for line_index, line in enumerate(narrow_file):
            tokens = line.strip().split()
            chrom, start, stop = tokens[0], int(tokens[1]), int(tokens[2])
            seq_length = stop - start
            
            if seq_length == 250:
                dna_sequence = genome_dict[chrom].get_slice(start - padding_each_side, stop + padding_each_side)
                if len(dna_sequence) != 350:
                    print('The DNA sequence is not 350 bp in length starting Will skip')
                    print('chrom: {} start: {} stop: {}'.format(chrom, start, stop))
                    continue
                seq_list.append(dna_sequence)
                header_line_list.append(start-padding_each_side)
                # posfile.write('>{}\n{}\n'.format(start - padding_each_side, dna_sequence))
                
                # coordinates are left inclusive right exclusive
                start_coord, end_coord = start - padding_each_side, stop + padding_each_side
                pos_coord_file.write('{} {} {}\n'.format(chrom, start_coord, end_coord))
                
                coord_dict[chrom].append(start_coord)
                
            if line_index % 1000 == 0:
                print('=> {}/{} = {:.2f}%'.format(line_index, line_count, line_index/line_count * 100), end='\r')

        print('=> Processed {} lines from {}'.format(line_count, filename))
        
        num_sequences = len(seq_list)
        test_indices = np.random.choice(num_sequences, int(num_sequences * 0.2), replace=False)

        test_list = [seq_list[i] for i in test_indices]
        test_header_list =[header_line_list[i] for i in test_indices]
        train_list = [seq_list[i] for i in range(num_sequences) if i not in test_list]
        train_header_list = [header_line_list[i] for i in range(num_sequences) if i not in test_list]

        print('train_list length: {}\n'
              'test_list length: {}'
              .format(len(train_list), len(test_list)))
        
        for start, train_seq in zip(train_header_list, train_list):
            train_posfile.write('>{}\n{}\n'.format(start - padding_each_side, train_seq))
        
        for start, test_seq in zip(test_header_list, test_list):
            test_posfile.write('>{}\n{}\n'.format(start - padding_each_side, test_seq))
            
        
    # for key in coord_dict.keys():
    #     coord_dict[key].sort()
        
    return coord_dict
        
    
    
if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('filename')
    args = arg_parser.parse_args()
    d = narrowpead_to_fa(args.filename)