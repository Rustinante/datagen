import argparse
import h5py
import numpy as np
import os
import time

from collections import defaultdict

from binary_search import search, get_location_from_line


class LineCache:
    def __init__(self):
        self.capacity = 1000
        self.count_per_eviction = 200
        self.cache = {}
        self.key_list = []
    
    def evict(self):
        for key in self.key_list[:self.count_per_eviction]:
            del self.cache[key]
        self.key_list = self.key_list[self.count_per_eviction:]
    
    def __contains__(self, item):
        return item in self.cache
    
    def __getitem__(self, item):
        return self.cache[item]
    
    def __setitem__(self, key, value):
        if len(self.cache) >= self.capacity:
            self.evict()
        self.cache[key] = value
        if key not in self.key_list:
            self.key_list.append(key)
    
    def __len__(self):
        return len(self.cache)
    
    def keys(self):
        return self.cache.keys()
    
    def values(self):
        return self.cache.values()
    
    def items(self):
        return self.cache.items()


def scan_through_line_for_number(alignment_file, start_line_hint, number):
    alignment_file.seek(start_line_hint)
    
    for line in alignment_file:
        location = get_location_from_line(line)
        
        if not location:
            raise ValueError('There still lines in the alignment file but cannot obtain coordinate location')
        
        if location == number:
            return line, start_line_hint
        
        # The lines are sorted so if the current location is already greater than the one
        # we're searching for we know what we search does not exist.
        elif location > number:
            return None
        
        start_line_hint += len(bytes(line, 'ascii'))
    
    return None


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for _ in file:
            count += 1
    return count


def distribute_tokens_to_species_sequence(tokens, species_sequence):
    for species_index, letter in enumerate(tokens):
        species_sequence[species_index] += letter
        
        
def extend_dataset(chr, purpose):
    checkpoint_time_str = time.strftime('%a %b %d %Y %H:%M:%S UTC%z', time.localtime(time.time()))
    print('Current time: {}'.format(checkpoint_time_str))
    
    cache = LineCache()
    
    array_list = []
    
    coordinate_filename = os.path.join('data', '{}_{}'.format(chr, purpose))
    alignment_filename = '{}_maf_sequence.csv'.format(chr)
    
    species_file_dict = {}
    
    print('=> coordinate_filename: {}'.format(coordinate_filename))
    print('=> alignment_filename: {}'.format(alignment_filename))

    total_line_count = get_line_count(coordinate_filename)
    
    dir_name = 'gkm_fasta_chr{}'.format(chr)
    os.makedirs(dir_name, exist_ok=True)
    
    with open(coordinate_filename, 'r') as file, open(alignment_filename, 'r') as alignment_file:
        
        header = alignment_file.readline().strip().split(',')
        del header[0]
        assert len(header) == 100
        human_index = header.index('hg19')
        
        for index, species_code in enumerate(header):
            species_filename = '{}_{}_{}.fasta'.format(species_code, chr, purpose)
            print('=> Creating {} under {}'.format(species_filename, dir_name))
            species_file_dict[index] = open(os.path.join(dir_name, species_filename), 'w')
        
        processed_line_count = 0
        start_time = time.time()
        flanking_number = 400
        
        for line in file:
            t1 = time.time()
            processed_line_count += 1
            (start_coordinate, sequence) = line.strip().split(',')
            start_coordinate = int(start_coordinate)
            
            species_sequence = defaultdict(str)
            
            start_line_hint = None
            
            for letter_index, coordinate in enumerate(range(start_coordinate - flanking_number, start_coordinate + 200 + flanking_number)):
                if coordinate in cache:
                    cached_result = cache[coordinate]
                    tokens = cached_result[0]
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
                    start_line_hint = cached_result[1]
                    continue
                
                elif not start_line_hint:
                    result = search(alignment_file, coordinate, alignment_filename)
                else:
                    result = scan_through_line_for_number(alignment_file=alignment_file,
                                                          start_line_hint=start_line_hint, number=coordinate)
                
                if result:
                    start_line_hint = result[1]
                    tokens = result[0].strip().split(',')
                    del tokens[0]
                    # assert len(tokens) == 100
                    
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
                    
                    cache[coordinate] = (tokens, start_line_hint)
                
                else:
                    print('Failed to find the letters at coordinate {}. '
                          'Substituting with N instead'.format(coordinate))
                    tokens = ['N'] * 100
                    distribute_tokens_to_species_sequence(tokens, species_sequence)
                    
            t2 = time.time()
            for species_index, _ in enumerate(header):
                species_file_dict[species_index].write('>middle 200 bp start coordinate {}\n'
                                                       '{}\n'
                                                       .format(start_coordinate, species_sequence[species_index]))
                
            t3 = time.time()
            print('t2-t1 {} t3-t2 {}'.format(t2-t1, t3-t2))
            
            if processed_line_count % 1000 == 0:
                elapsed_time = time.time() - start_time
                time_per_line = elapsed_time / processed_line_count
                print('Processed [{}/{}] lines in {:5f}s, averaging: {:5f}s per line, {:3f} done.'
                      .format(processed_line_count, total_line_count,
                              elapsed_time, time_per_line, processed_line_count/total_line_count))
        
    for filename, species_filename in species_file_dict.values():
        print('=> Closing {}'.format(filename))
        species_filename.close()
        
    print('=> Done!')
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('chr')
    parser.add_argument('purpose')
    args = parser.parse_args()
    extend_dataset(args.chr, args.purpose)
