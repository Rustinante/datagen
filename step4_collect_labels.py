from collections import defaultdict
import numpy as np
import h5py
from scipy import io
import time


class Timer:
    def __init__(self):
        self.stamp = time.time()

    def reset(self):
        self.stamp = time.time()

    def get_elapsed_time(self):
        return time.time() - self.stamp


class ProgressTracker(Timer):
    def __init__(self, frequency):
        super().__init__()
        self.frequency = frequency

    def print_progress(self, line_index):
        if line_index % self.frequency == 0:
            elapsed = self.get_elapsed_time()
            processed_line_count = line_index + 1
            print(f'processed {processed_line_count} lines in {elapsed:5f}s, '
                  f'averaging {elapsed / processed_line_count:5f}s per line')


def collect_labels():
    print('=> collecting training labels')

    labels_dict = defaultdict(list)
    tracker = ProgressTracker(10000)

    with open('train_coord', 'r') as train_coord_file, h5py.File('train.mat', 'r') as train_file:
        training_labels = train_file['traindata']
        print(f'training labels has shape: {training_labels.shape}')

        for line_index, line in enumerate(train_coord_file):
            tokens = line.split()
            chrom, start, stop = tokens[0], int(tokens[1]), int(tokens[2])
            labels_dict[chrom].append(training_labels[:, line_index])
            tracker.print_progress(line_index)

    return labels_dict


def add_labels_to_dataset(labels_dict):
    for chrom, label_list in labels_dict.items():
        label_array = np.tile(label_list, (2, 1))
        print(f'=> {chrom} has {len(label_list)} labels with a corresponding array of shape {label_array.shape}')

        hdf5_filename = f'{chrom}_train.hundred.hdf5'
        print(f'=> opening {hdf5_filename}')

        with h5py.File(hdf5_filename, 'r+') as hdf5_file:
            group = hdf5_file.create_group('label')
            group.create_dataset('data', data=label_array, dtype='uint8')
        print(f'=> added the labels to {hdf5_filename}')


def collect_validation_labels():
    print('=> collecting validation labels')

    labels_dict = defaultdict(list)
    tracker = ProgressTracker(10000)

    validation_file = io.loadmat('valid.mat')

    with open('valid_coord', 'r') as valid_coord_file:
        validation_labels = validation_file['validdata']
        print(f'validation labels has shape: {validation_labels.shape}')

        for line_index, line in enumerate(valid_coord_file):
            tokens = line.split()
            chrom = tokens[0]
            labels_dict[chrom].append(validation_labels[line_index, :])
            tracker.print_progress(line_index)

    return labels_dict


def add_labels_to_validation_dataset(labels_dict):
    for chrom, label_list in labels_dict.items():
        label_array = np.tile(label_list, (2, 1))
        print(f'=> {chrom} has {len(label_list)} labels with a corresponding array of shape {label_array.shape}')

        hdf5_filename = f'{chrom}_valid.hundred.hdf5'
        print(f'=> opening {hdf5_filename}')

        with h5py.File(hdf5_filename, 'r+') as hdf5_file:
            group = hdf5_file.create_group('label')
            group.create_dataset('data', data=label_array, dtype='uint8')
        print(f'=> added the labels to {hdf5_filename}')


def collect_test_labels():
    print('=> collecting test labels')
    labels_dict = defaultdict(list)
    tracker = ProgressTracker(10000)

    test_file = io.loadmat('test.mat')
    test_labels = test_file['testdata']
    print(f'testing labels has shape {test_labels.shape}')

    with open('test_coord', 'r') as test_coord_file:
        for line_index, line in enumerate(test_coord_file):
            tokens = line.split()
            chrom = tokens[0]
            labels_dict[chrom].append(test_labels[line_index, :])
            tracker.print_progress(line_index)

    return labels_dict


def add_labels_to_test_dataset(labels_dict):
    for chrom, label_list in labels_dict.items():
        label_array = np.tile(label_list, (2, 1))
        print(f'=> {chrom} has {len(label_list)} labels with a corresponding array of shape {label_array.shape}')

        hdf5_filename = f'{chrom}_test.hundred.hdf5'
        print(f'=> opening {hdf5_filename}')

        with h5py.File(hdf5_filename, 'r+') as hdf5_file:
            group = hdf5_file.create_group('label')
            group.create_dataset('data', data=label_array, dtype='uint8')
        print(f'=> added the labels to {hdf5_filename}')


if __name__ == '__main__':
    print('=> adding labels to the training dataset')
    collected_train_labels = collect_labels()
    add_labels_to_dataset(collected_train_labels)

    collected_validation_labels = collect_validation_labels()
    print('=> adding labels to the validation dataset')
    add_labels_to_validation_dataset(collected_validation_labels)

    collected_test_labels = collect_test_labels()
    print('=> adding labels to the testing dataset')
    add_labels_to_test_dataset(collected_test_labels)
