import argparse
import h5py
import numpy as np
import os


def split_pos_neg_for_gkm(chr, purpose, label_index):
    def sample_index(index):
        return index // 2
    
    label_index = int(label_index)
    gkm_data_dirname = 'gkm_fasta'
    pure_label_filename = os.path.join('pure_labels', '{}_{}.pure_label.hdf5'.format(chr, purpose))
    feature_filename = os.path.join(gkm_data_dirname, '{}_{}'.format(chr, purpose))
    pos_filename = 'gkm_{}_{}.pos.fa'.format(chr, purpose)
    neg_filename = 'gkm_{}_{}.neg.fa'.format(chr, purpose)
    pos_filepath = os.path.join(gkm_data_dirname, pos_filename)
    neg_filepath = os.path.join(gkm_data_dirname, neg_filename)
    
    print('=> labels are from {} at index {}\n'
          '=> features are from {}\n'
          '=> pos_filename: {}\n'
          '=> neg_filename: {}\n'
          '=> label_index: {}'
          .format(pure_label_filename, label_index, feature_filename, pos_filename, neg_filename, label_index))
    
    with h5py.File(pure_label_filename, 'r') as file, open(feature_filename, 'r') as feature_file,\
            open(pos_filepath, 'w') as pos_file, open(neg_filepath, 'w') as neg_file:
        label = file['label/data']
        single_label_column = label[:, label_index]
        positive_sample_indices = np.nonzero(single_label_column)[0]
        
        print('=> There are {} positive sequences out of {}'.format(len(positive_sample_indices), len(single_label_column)))
        
        processed_line_count = 0
        for line_index, line in enumerate(feature_file):
            # coordinate = int(line[1:].strip())
            # print(coordinate)
            actual_index = sample_index(line_index)
            if actual_index in positive_sample_indices:
                pos_file.write(line)
            else:
                neg_file.write(line)
            processed_line_count += 1
            
    print('=> processed {} for a total of {} lines'.format(feature_filename, processed_line_count))


if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('chr')
    arg_parser.add_argument('purpose')
    arg_parser.add_argument('label_index')
    args = arg_parser.parse_args()
    split_pos_neg_for_gkm(args.chr, args.purpose, args.label_index)