import numpy as np
import h5py
import argparse
import time


def generate_five_channel_reverse_complement(chr, purpose):
    complement_map = {
        'a': 't',
        't': 'a',
        'c': 'g',
        'g': 'c',
        'A': 'T',
        'T': 'A',
        'C': 'G',
        'G': 'C',
        'N': 'N',
        'n': 'n',
        'X': 'X',
        'x': 'x'
    }
    
    a_array = np.array([1, 0, 0, 0, 0], dtype='uint8')
    g_array = np.array([0, 1, 0, 0, 0], dtype='uint8')
    c_array = np.array([0, 0, 1, 0, 0], dtype='uint8')
    t_array = np.array([0, 0, 0, 1, 0], dtype='uint8')
    x_array = np.array([0, 0, 0, 0, 1], dtype='uint8')
    zero_array = np.array([0, 0, 0, 0, 0], dtype='uint8')
    
    hdf5_filename = '{}_{}.five_channel.mix.hdf5'.format(chr, purpose)
    feature_key = 'feature/data'
    label_key = 'label/data'
    
    number_of_speices = 5
    number_of_channels = 5
    
    rev_comp_filename = '{}_{}.rev_comp.five_channel.mix.hdf5'.format(chr, purpose)
    
    print('opening {} and creating {}'.format(hdf5_filename, rev_comp_filename))
    
    with h5py.File(hdf5_filename, 'r') as hdf5_file, h5py.File(rev_comp_filename) as revcomp_file:
        original_feature = hdf5_file[feature_key]
        print('feature has shape: ', original_feature.shape)
        num_samples = original_feature.shape[0]
        
        feature_group = revcomp_file.create_group('feature')
        # new_feature is the reverse complement data
        new_feature = feature_group.create_dataset('data', (num_samples, number_of_speices, 1000, number_of_channels), dtype='uint8')
        
        print('=> copying the features')
        stamp = time.time()
        
        for sample_index in range(num_samples):
            sample = original_feature[sample_index]
            # sample should be of size number_of_speices x 1000 x number_of_channels
            revcomp_matrix = np.zeros_like(sample, dtype='uint8')
            
            for species_index in range(number_of_speices):
                species_row = revcomp_matrix[species_index]
                
                for nucleotide_index in range(1000):
                    original_encoding = sample[species_index][nucleotide_index]
                    
                    if np.array_equal(original_encoding, a_array):
                        species_row[999 - nucleotide_index] = t_array
                    
                    elif np.array_equal(original_encoding, g_array):
                        species_row[999 - nucleotide_index] = c_array
                    
                    elif np.array_equal(original_encoding, c_array):
                        species_row[999 - nucleotide_index] = g_array
                    
                    elif np.array_equal(original_encoding, t_array):
                        species_row[999 - nucleotide_index] = a_array
                        
                    elif np.array_equal(original_encoding, x_array):
                        species_row[999 - nucleotide_index] = x_array
                    
                    else:
                        # all zeros
                        species_row[999 - nucleotide_index] = zero_array
            
            new_feature[sample_index] = revcomp_matrix
            
            if sample_index % 1000 == 1:
                print('{}/{} {:5f} done in {:2f}s'
                      .format(sample_index, num_samples, sample_index / num_samples, time.time() - stamp))
        
        print('=> copying the labels')
        label_group = revcomp_file.create_group('label')
        label_array = hdf5_file[label_key]
        label_group.create_dataset('data', data=label_array, dtype='uint8')
        print('=> copied the labels')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('chr')
    parser.add_argument('purpose')
    args = parser.parse_args()
    generate_five_channel_reverse_complement(args.chr, args.purpose)