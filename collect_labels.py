from collections import defaultdict
import numpy as np
import h5py
import time


def collect_labels():
    print('=> collecting labels')
    
    labels_dict = defaultdict(list)
    stamp = time.time()
    
    with open('train_coord', 'r') as train_coord_file, h5py.File('train.mat', 'r') as train_file:
        training_labels = train_file['traindata']
        print('training labels has shape: {}'.format(training_labels.shape))
        
        for line_index, line in enumerate(train_coord_file):
            tokens = line.split()
            chr, start, stop = tokens[0], int(tokens[1]), int(tokens[2])
            labels_dict[chr].append(training_labels[:, line_index])
            
            if line_index % 10000 == 9999:
                elapsed = time.time() - stamp
                processed_line_count = line_index + 1
                print('processed {} lines in {:5f}s, averaging {:5f}s per line'
                      .format(processed_line_count, elapsed, elapsed / processed_line_count))
            
    return labels_dict


def add_labels_to_dataset(labels_dict):
    for chr, label_list in labels_dict.items():
        num_labels = len(label_list)
        label_array = np.array(label_list)
        print('=> {} has {} labels with a corresponding array of shape {}'.format(chr, num_labels, label_array.shape))
        hdf5_filename = '{}_train.align.hdf5'.format(chr)
        print('=> opening {}'.format(hdf5_filename))
        with h5py.File(hdf5_filename, 'r+') as hdf5_file:
            group = hdf5_file.create_group('label')
            dataset = group.create_dataset('data', data=label_array, dtype='uint8')
        print('=> added the labels to {}'.format(hdf5_filename))
            

if __name__ == '__main__':
    collected_labels = collect_labels()
