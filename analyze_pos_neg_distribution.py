import argparse
import h5py
import numpy as np
import os


def analyze_pos_neg_distribution():
    purpose = 'train'
    chrom_list = [
        'chr1',
        'chr2',
        'chr3',
        'chr4',
        'chr5',
        'chr6',
        'chr7',
        'chr10',
        'chr11',
        'chr12',
        'chr13',
        'chr14',
        'chr15',
        'chr16',
        'chr17',
        'chr18',
        'chr19',
        'chr20',
        'chr21',
        'chr22',
    ]
    for chrom in chrom_list:
        pure_label_filename = os.path.join('pure_labels', '{}_{}.pure_label.hdf5'.format(chrom, purpose))
        
        print('\n'
              '=> {}\n'
              '=> {}\n'
              '=> features are from {}\n'
              .format(chrom, purpose, pure_label_filename))
        
        stats = {}
        
        with h5py.File(pure_label_filename, 'r') as file:
            label = file['label/data']
            num_chromatin_features = label.shape[1]
            
            print('=> {} chromatin features in total'.format(num_chromatin_features))
            
            for chrom_feature_index in range(num_chromatin_features):
                single_label_column = label[:, chrom_feature_index]
                positive_sample_indices = np.nonzero(single_label_column)[0]
            
                num_pos_samples = len(positive_sample_indices)
                total_num_samples = len(single_label_column)
                postive_ratio = num_pos_samples / total_num_samples
                
                stats[(chrom, chrom_feature_index)] = postive_ratio
                
                print('=> [{}] positives/total = {}/{} = {:2f}'
                      .format(chrom_feature_index,
                              num_pos_samples,
                              total_num_samples,
                              postive_ratio))
            

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser()
    args = arg_parser.parse_args()
    analyze_pos_neg_distribution()