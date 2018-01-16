import h5py
import argparse


def extract_primates(chr, purpose):
    primate_indices = [0, 10, 18, 40, 49, 61, 71, 74, 75, 78, 83, 85]
    
    hdf5_filename = '{}_{}.align.hdf5'.format(chr, purpose)
    feature_key = 'feature/data'
    label_key = 'label/data'
    
    primate_filename = '{}_{}.primate.hdf5'.format(chr, purpose)
    print('opening {}'.format(primate_filename))
    
    with h5py.File(hdf5_filename, 'r') as hdf5_file, h5py.File(primate_filename) as primate_file:
        feature = hdf5_file[feature_key]
        print('feature has shape: ', feature.shape)
        num_samples = feature.shape[0]
        feature_group = primate_file.create_group('feature')
        primate_feature_data = feature_group.create_dataset('data', (num_samples, len(primate_indices), 1000, 4), dtype='uint8')
        
        print('=> copying the features')
        for i in range(num_samples):
            primate_feature_data[i] = feature[i][primate_indices]
            if i % 1000 == 1:
                print('{}/{} {:5f}'.format(i, num_samples, i / num_samples))
        
        label_group = primate_file.create_group('label')
        label_array = hdf5_file[label_key]
        label_group.create_dataset('data', data=label_array, dtype='uint8')
        print('=> copied the labels')
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('chr')
    parser.add_argument('purpose')
    args = parser.parse_args()
    extract_primates(args.chr, args.purpose)