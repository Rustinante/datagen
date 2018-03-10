import argparse
from collections import defaultdict
import os
import time

from binary_search import search, scan_through_line_for_number
from LineCache import LineCache


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
    
    filenames = [(get_alignment_filename(chrom), chrom) for chrom in list(map(convert_number_to_chrom_str,range(1, 22))) + ['chrX']]
    print(f'filenames: {filenames!r}')
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


def get_species_letters_from_coord(coord_filename, purpose):
    checkpoint_time_str = time.strftime('%a %b %d %Y %H:%M:%S UTC%z', time.localtime(time.time()))
    print('Current time: {}'.format(checkpoint_time_str))
    
    cache = LineCache()
    
    species_file_dict = {}
    
    print(f'=> coordinate_filename: {coord_filename}')
    
    total_line_count = get_line_count(coord_filename)
    
    dir_name = f'{coord_filename}.{purpose}.mult_species'
    os.makedirs(dir_name, exist_ok=True)
    
    alignment_file_dict, header = open_alignment_files()
    
    with open(coord_filename, 'r') as coord_file:
        for index, species_code in enumerate(header):
            species_filename = f'{index}_{species_code}_{purpose}.fa'
            print('=> Creating {} under {}'.format(species_filename, dir_name))
            species_file_dict[index] = open(os.path.join(dir_name, species_filename), 'w')
        
        processed_line_count = 0
        start_time = time.time()
        
        for line in coord_file:
            processed_line_count += 1
            tokens = line.strip().split()
            chrom, start_coord, stop_coord = tokens[0], int(tokens[1]), int(tokens[2])
            alignment_file = alignment_file_dict[chrom]
            
            species_sequence = defaultdict(str)
            
            start_line_hint = None
            for letter_index, coordinate in enumerate(range(start_coord, stop_coord)):
                if coordinate in cache:
                    cached_result = cache[coordinate]
                    tokens = cached_result[0]
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
                    start_line_hint = cached_result[1]
                    continue
                
                elif not start_line_hint:
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
                    
                    cache[coordinate] = (tokens, start_line_hint)
                
                else:
                    print(f'Failed to find the letters at coordinate {coordinate}. '
                          f'Substituting with N instead')
                    tokens = ['N'] * (stop_coord - start_coord)
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
            
            for species_index in range(len(header)):
                species_file_dict[species_index].write(f'>{chrom} {start_coord} {stop_coord}\n'
                                                       f'{species_sequence[species_index]}\n')
            
            if processed_line_count % 1000 == 0:
                elapsed_time = time.time() - start_time
                time_per_line = elapsed_time / processed_line_count
                print(f'Processed [{processed_line_count}/{total_line_count}] lines in {elapsed_time:5f}s, '
                      f'averaging: {time_per_line:5f}s per line, {processed_line_count / total_line_count:3f} done.')
    
    for filename, species_file in species_file_dict.items():
        print(f'=> Closing {filename}')
        species_file.close()
        
    close_alignment_files(alignment_file_dict)
    
    print('=> Done!')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('coord_file')
    parser.add_argument('purpose')
    args = parser.parse_args()
    get_species_letters_from_coord(args.coord_file, args.purpose)