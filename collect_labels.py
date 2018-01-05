from collections import defaultdict
import numpy
import h5py
import time


def collect_labels():
    labels_dict = defaultdict(list)
    stamp = time.time()
    print('=> collecting labels')
    with open('train_coord', 'r') as train_coord_file, h5py.File('train.mat', 'r') as train_file:
        labels = train_file['traindata']
        print('labels has shape: {}'.format(labels.shape))
        for line_index, line in enumerate(train_coord_file):
            tokens = line.split()
            chr, start, stop = tokens[0], int(tokens[1]), int(tokens[2])
            labels_dict[chr].append(labels[:, line_index])
            
            if line_index % 1000 == 1:
                elapsed = time.time() - stamp
                print('processed {} lines in {}s, averagin {}s per line'
                      .format(line_index + 1, elapsed, elapsed / (line_index + 1)))
            
    return labels_dict


# def add_labels_to_dataset(labels_dict):
#     opened_files = {}
#     data_groups = {}
#     for chr, value in labels_dict.items():
#         if chr not in opened_files:
#             print('opening {}'.format())
            # opened_files[chr] = h5py.File('{}_train.align.hdf5'.format(chr), 'r+')
            # data_groups[chr] = opened_files[chr].create_group('label')
            # data_groups[chr].create_dataset('data', (len(array_list), 100, 1000, 4), dtype='uint8', compression='gzip')
    #
    # for chr_index, file in opened_files.items():
    #     print('closing {}'.format(get_destination_filename(chr_index, data_purpose)))
    #     file.close()
        

if __name__ == '__main__':
    collected_labels = collect_labels()
