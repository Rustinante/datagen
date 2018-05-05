import argparse
from collections import defaultdict
import os
import time

from binary_search import search, scan_through_line_for_number


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


def distribute_tokens_to_species_sequence(tokens, species_sequence):
    for species_index, letter in enumerate(tokens):
        species_sequence[species_index] += letter


def get_alignment_filename(chrom):
    return f'{chrom}_maf_sequence.csv'


def open_alignment_files():
    def convert_number_to_chrom_str(number):
        return 'chr' + str(number)
    
    filenames = [(get_alignment_filename(chrom), chrom) for chrom in
                 list(map(convert_number_to_chrom_str, range(1, 23))) + ['chrX', 'chrY']]
    print(f'-> alignment filenames: {filenames!r}')
    file_dict = {}
    species_header = None
    for filename, chrom in filenames:
        file = open(filename, 'r')
        header = file.readline().strip().split(',')
        del header[0]
        assert len(header) == 100
        file_dict[chrom] = file
        
        species_header = header
    
    return file_dict, species_header


def close_alignment_files(file_dict):
    for file in file_dict.values():
        file.close()


def is_informative_sequence(seq):
    for letter in seq:
        if letter != 'X' and letter != 'N':
            return True
    return False


def generate_informative_species_filename(species_index, species_code):
    return f'{species_index}_{species_code}.informative'


def get_species_letters_from_coord(coord_filename, target_dirname, ignore_noninformative):
    checkpoint_time_str = time.strftime('%a %b %d %Y %H:%M:%S UTC%z', time.localtime(time.time()))
    print('Current time: {}'.format(checkpoint_time_str))
    print(f'-> ignore noninformative: {ignore_noninformative}')
    
    species_file_dict = {}
    informative_seq_file_dict = {}
    
    total_line_count = get_line_count(coord_filename)
    print(f'=> coordinate_filename: {coord_filename}\n'
          f'-> {coord_filename} has {total_line_count} lines')
    
    alignment_file_dict, header = open_alignment_files()
    
    with open(coord_filename, 'r') as coord_file:
        dir_name = os.path.join(target_dirname, f'{os.path.basename(os.path.normpath(coord_filename))}.mult_species')
        os.makedirs(dir_name, exist_ok=False)
        
        for index, species_code in enumerate(header):
            species_filename = f'{index}_{species_code}.fa.ir'
            print('=> Creating {} under {}'.format(species_filename, dir_name))
            species_file_dict[index] = open(os.path.join(dir_name, species_filename), 'w')
            
            informative_seq_filename = generate_informative_species_filename(index, species_code)
            informative_seq_file_dict[index] = open(os.path.join(dir_name, informative_seq_filename), 'w')
        
        processed_line_count = 0
        start_time = time.time()
        number_of_n_substituted = 0
        
        for line_index, line in enumerate(coord_file):
            processed_line_count += 1
            tokens = line.strip().split()
            chrom, start_coord, stop_coord = tokens[0], int(tokens[1]), int(tokens[2])
            alignment_file = alignment_file_dict[chrom]
            
            species_sequence = defaultdict(str)
            
            start_line_hint = None
            for letter_index, coordinate in enumerate(range(start_coord, stop_coord)):
                if not start_line_hint:
                    result = search(alignment_file, coordinate, get_alignment_filename(chrom))
                else:
                    result = scan_through_line_for_number(alignment_file=alignment_file,
                                                          start_line_hint=start_line_hint, number=coordinate)
                
                if result:
                    start_line_hint = result[1]
                    tokens = result[0].strip().split(',')
                    del tokens[0]
                    assert len(tokens) == 100
                    
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
                
                else:
                    number_of_n_substituted += 1
                    tokens = ['N'] * (stop_coord - start_coord)
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
            
            if ignore_noninformative:
                for species_index in range(len(header)):
                    seq = species_sequence[species_index]
                    if is_informative_sequence(seq):
                        species_file_dict[species_index].write(f'>{chrom} {start_coord} {stop_coord}\n'
                                                               f'{seq}\n')
                        informative_seq_file_dict[species_index].write(f'{line_index}\n')
            else:
                for species_index in range(len(header)):
                    seq = species_sequence[species_index]
                    species_file_dict[species_index].write(f'>{chrom} {start_coord} {stop_coord}\n'
                                                           f'{seq}\n')
                    if is_informative_sequence(seq):
                        informative_seq_file_dict[species_index].write(f'{line_index}\n')
            
            if processed_line_count % 1000 == 0:
                elapsed_time = time.time() - start_time
                time_per_line = elapsed_time / processed_line_count
                print(f'Processed [{processed_line_count}/{total_line_count}] lines in {elapsed_time:5f}s, '
                      f'averaging: {time_per_line:5f}s per line, {processed_line_count / total_line_count:3f} done.')
    
    for filename, species_file in species_file_dict.items():
        print(f'=> Closing {filename}')
        species_file.close()
    
    for filename, file in informative_seq_file_dict.items():
        print(f'=> Closing {filename}')
        file.close()
    
    close_alignment_files(alignment_file_dict)
    
    print(f'-> #N substituted: {number_of_n_substituted}\n'
          f'=> Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('coord_file')
    parser.add_argument('target_dirname')
    parser.add_argument('--ignore-noninformative', action='store_true')
    args = parser.parse_args()
    get_species_letters_from_coord(args.coord_file, args.target_dirname, args.ignore_noninformative)
