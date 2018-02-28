import numpy as np
import argparse


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


def generate_test_data(filename, test_dataset_ratio=0.2):
    train_filename = 'roci.train.' + filename
    test_filename = 'roci.test.' + filename

    total_line_count = get_line_count(filename)
    
    print('train_filename: {}\n'
          'test_filename: {}'
          .format(train_filename, test_filename))
    
    test_data_indices = np.random.choice(total_line_count, int(total_line_count * test_dataset_ratio), replace=False)
    # multiply the index by 2 because one data point in the fasta file has two lines.
    test_data_indices *= 2
    
    with open(filename, 'r') as original_file, open(test_filename, 'w') as test_file, open(train_filename, 'w') as train_file:
        should_write = False
        for index, line in enumerate(original_file):
            if index in test_data_indices:
                test_file.write(line)
                should_write = True
            elif should_write:
                test_file.write(line)
                should_write = False
            else:
                train_file.write(line)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pos_filename')
    parser.add_argument('neg_filename')
    args = parser.parse_args()

    generate_test_data(args.pos_filename)
    generate_test_data(args.neg_filename)
