import argparse
import h5py
import numpy as np
import os


def split_pos_neg_for_gkm(chr, purpose, label_index):
    """
    This relies on 2_coord_to_fasta to prepare the corresponding fasta files that are not
    split into the positives and the negatives.
    The output of 2_coord_to_fasta should be in the directory feature_dirname.
    :param chr: the chromosome name, e.g. chr1
    :param purpose: train, valid or test
    :param label_index: The index of the chromatin features out of the 919 from the DeepSEA data
    """
    def sample_index(index):
        return index // 2
    
    label_index = int(label_index)
    feature_dirname = 'hg19_fasta_aggregate'
    pure_label_filename = os.path.join('pure_labels', '{}_{}.pure_label.hdf5'.format(chr, purpose))
    feature_filename = os.path.join(feature_dirname, '{}_{}'.format(chr, purpose))
    
    pos_filename = 'gkm_{}_{}.pos.fa'.format(chr, purpose)
    neg_filename = 'gkm_{}_{}.neg.fa'.format(chr, purpose)
    
    output_dirname = 'hg19_chromatin_feature_{}_fasta'.format(label_index)
    pos_filepath = os.path.join(output_dirname, pos_filename)
    neg_filepath = os.path.join(output_dirname, neg_filename)
    
    if not os.path.isdir(output_dirname):
        print('\n=> Creating directory {} as it does not exist.'.format(output_dirname))
        os.makedirs(output_dirname, exist_ok=False)
    
    print('\n'
          '=> labels are from {} at index {}\n'
          '=> features are from {}\n'
          '=> pos_filename: {}\n'
          '=> neg_filename: {}\n'
          '=> label_index: {}'
          .format(pure_label_filename, label_index, feature_filename, pos_filename, neg_filename, label_index))
    
    with h5py.File(pure_label_filename, 'r') as file, open(feature_filename, 'r') as feature_file,\
            open(pos_filepath, 'w') as pos_file, open(neg_filepath, 'w') as neg_file:
        label = file['label/data']
        single_label_column = label[:, label_index]
        pos_sample_indices = np.nonzero(single_label_column)[0]
        
        num_pos_samples = len(pos_sample_indices)
        total_num_samples = len(single_label_column)
        num_neg_samples = total_num_samples - num_pos_samples
        
        neg_sample_indices = np.array([i for i in range(total_num_samples) if i not in pos_sample_indices])
        
        check_neg_indices = False
        
        multiplier = 2
        if num_pos_samples * multiplier < num_neg_samples:
            check_neg_indices = True
            print('=> # negative samples more than {} times the # positive samples.\n'
                  '=> Will randomly sample the negative samples.'.format(multiplier))
            random_neg_indices = np.random.choice(num_neg_samples, size=num_pos_samples, replace=False)
            neg_sample_indices = neg_sample_indices[random_neg_indices]
            num_neg_samples = len(neg_sample_indices)
        
        print('=> {} positives {} negatives {} total\n'
              '=> Got rid of {} negatives'
              .format(num_pos_samples,
                      num_neg_samples,
                      total_num_samples,
                      total_num_samples - num_pos_samples - num_neg_samples))
        
        processed_line_count = 0
        
        if check_neg_indices:
            for line_index, line in enumerate(feature_file):
                actual_index = sample_index(line_index)
                if actual_index in pos_sample_indices:
                    pos_file.write(line)
                elif actual_index in neg_sample_indices:
                    neg_file.write(line)
                processed_line_count += 1
        else:
            for line_index, line in enumerate(feature_file):
                actual_index = sample_index(line_index)
                if actual_index in pos_sample_indices:
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