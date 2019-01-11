import argparse
import os
import time

import biotool.file_binary_search as file_binary_search
import h5py
import numpy as np

from binary_search import scan_through_line_for_number
from line_cache import LineCache
from nucleotide_mapping import mapping, complement_mapping, map_counts_to_vec, map_counts_to_revcomp_vec
from util.file import get_line_count


def extend_dataset(chrom, purpose):
    checkpoint_time_str = time.strftime('%a %b %d %Y %H:%M:%S UTC%z', time.localtime(time.time()))
    print('Current time: {}'.format(checkpoint_time_str))

    cache = LineCache()

    matrix_list = []
    revcomp_matrix_list = []

    coordinate_filename = os.path.join('data', '{}_{}'.format(chrom, purpose))
    alignment_filename = '{}_maf_sequence.csv'.format(chrom)
    hdf5_filename = '{}_{}.counting.hdf5'.format(chrom, purpose)
    num_rows = 2
    num_non_humans = 99

    print('=> coordinate_filename: {}'.format(coordinate_filename))
    print('=> alignment_filename: {}'.format(alignment_filename))
    print('=> target hdf5_filename: {}'.format(hdf5_filename))
    file_byte_size = os.stat(alignment_filename).st_size

    line_count = get_line_count(coordinate_filename)

    flanking_number = 400
    seq_len = 200 + 2 * flanking_number
    last_bp_index = seq_len - 1
    feature_dim = 5

    serializing_index = 0

    with open(coordinate_filename, 'r') as coord_file, open(alignment_filename, 'r') as alignment_file, \
            h5py.File(hdf5_filename, 'w') as hdf5_file:

        feature_group = hdf5_file.create_group('feature')
        feature_data = feature_group.create_dataset('data', (line_count * 2, seq_len, num_rows, feature_dim),
                                                    dtype='uint8')

        header = alignment_file.readline().strip().split(',')
        human_index = header.index('hg19')

        processed_line_count = 0
        start_time = time.time()

        for line in coord_file:
            processed_line_count += 1
            (start_coordinate, sequence) = line.strip().split(',')
            start_coordinate = int(start_coordinate)

            # 1000 x 2 x 5
            seq_matrix = np.zeros((seq_len, num_rows, feature_dim), dtype='uint8')
            revcomp_seq_matrix = np.zeros((seq_len, num_rows, feature_dim), dtype='uint8')

            start_line_hint = None
            for bp_index, (hg_letter, coordinate) in enumerate(
                    zip(sequence, range(start_coordinate - flanking_number, start_coordinate + 200 + flanking_number))):

                if coordinate in cache:
                    # the entry for revcomp_alignment_matrix in the cache corresponds to the
                    # coordinate last_seq_index - bp_index on the reverse complement strand
                    (seq_matrix[bp_index, :, :],
                     revcomp_seq_matrix[last_bp_index - bp_index, :, :],
                     start_line_hint) = cache[coordinate]
                    continue

                # insert the human letter first
                seq_matrix[bp_index, 0, :] = mapping[hg_letter]
                revcomp_seq_matrix[last_bp_index - bp_index, 0, :] = complement_mapping[hg_letter]

                if not start_line_hint:
                    # TODO: test they are equal
                    # result = search(alignment_file, coordinate, file_byte_size)
                    result = file_binary_search.search(alignment_filename, coordinate, num_header_lines=1)
                else:
                    result = scan_through_line_for_number(alignment_file=alignment_file,
                                                          start_line_hint=start_line_hint, number=coordinate)

                if result:
                    tokens, start_line_hint = result[0].strip().split(','), result[1]
                    # Important: pop the human_index first before removing the start index,
                    # so that the human_index will be the correct index.
                    del tokens[human_index]
                    del tokens[0]

                    a = g = c = t = x = 0
                    for remaining_token_index, letter in enumerate(tokens):
                        a += letter.upper() == 'A'
                        g += letter.upper() == 'G'
                        c += letter.upper() == 'C'
                        t += letter.upper() == 'T'
                        x += letter.upper() == 'X'

                    assert a + g + c + t + x == num_non_humans
                    seq_matrix[bp_index, 1, :] = map_counts_to_vec(a=a, g=g, c=c, t=t, x=x)
                    revcomp_seq_matrix[last_bp_index - bp_index, 1, :] = map_counts_to_revcomp_vec(a=a, g=g, c=c, t=t,
                                                                                                   x=x)

                    cache[coordinate] = (seq_matrix[bp_index],
                                         revcomp_seq_matrix[last_bp_index - bp_index],
                                         start_line_hint)

                else:
                    # broadcasting along the species dimension
                    seq_matrix[bp_index, 1, :] = map_counts_to_vec(a=0, g=0, c=0, t=0, x=num_non_humans)
                    revcomp_seq_matrix[last_bp_index - bp_index, 1, :] = map_counts_to_revcomp_vec(a=0, g=0, c=0, t=0,
                                                                                                   x=num_non_humans)

            matrix_list.append(seq_matrix)
            revcomp_matrix_list.append(revcomp_seq_matrix)

            if processed_line_count % 100 == 1:
                for matrix, revcomp_matrix in zip(matrix_list, revcomp_matrix_list):
                    feature_data[serializing_index] = matrix
                    # For the reverse complement strand
                    feature_data[serializing_index + line_count] = revcomp_matrix
                    serializing_index += 1

                matrix_list = []
                revcomp_matrix_list = []

                elapsed_time = time.time() - start_time
                time_per_line = elapsed_time / processed_line_count
                print(
                    f'{processed_line_count}/{line_count} = {processed_line_count / line_count:.2%} in {elapsed_time:.4f}s '
                    f'averaging {time_per_line:2f}s per line '
                    f'current serializing index: {serializing_index}',
                    end='\r')

        # Serializing the remaining data
        if matrix_list:
            assert len(matrix_list) == len(revcomp_matrix_list)
            for matrix, revcomp_matrix, in zip(matrix_list, revcomp_matrix_list):
                feature_data[serializing_index] = matrix
                # For the reverse complement strand
                feature_data[serializing_index + line_count] = revcomp_matrix
                serializing_index += 1
            elapsed_time = time.time() - start_time
            time_per_line = elapsed_time / processed_line_count
            print(
                f'{processed_line_count}/{line_count} = {processed_line_count / line_count:.2%} in {elapsed_time:.4f}s '
                f'averaging {time_per_line:2f}s per line '
                f'current serializing index: {serializing_index}')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('chr')
    parser.add_argument('purpose')
    args = parser.parse_args()
    extend_dataset(args.chr, args.purpose)
