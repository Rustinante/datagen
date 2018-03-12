import shutil
import os
import argparse
import numpy as np


def get_line_count(filename):
    count = 0
    with open(filename, 'r') as file:
        for line in file:
            if line:
                count += 1
    return count


def map_line_to_coord_triple(line):
    tokens = line.strip().split()
    # In the fasta files > precedes the label line
    if tokens[0][0] == '>':
        return tokens[0][1:], int(tokens[1]), int(tokens[2])
    else:
        return tokens[0], int(tokens[1]), int(tokens[2])


def downsample(in_filename, out_filename, downsample_coord_filename):
    coord_set = set()
    with open(downsample_coord_filename, 'r') as coord_file:
        for line in coord_file:
            chrom, start, stop = map_line_to_coord_triple(line)
            coord_set.add((chrom, start, stop))
            
    with open(in_filename, 'r') as infile,  open(out_filename, 'w') as outfile:
        for index, line in enumerate(infile):
            if map_line_to_coord_triple(line) in coord_set:
                outfile.write(line)
                outfile.write(infile.readline())
            else:
                # Discard the current line and the next line
                infile.readline()


def create_target_filename(filename, purpose, label):
    index = filename.split('_')[0]
    return f'{index}_{purpose}_{label}.fa'


def get_filename_with_downsample_suffix(filename, downsample_ratio):
    return f'{filename}.at{downsample_ratio:.2f}'
    

def generate_downsample_coord(coord_filename, out_ratio):
    print(f'=> Downdsampling coordinates from {coord_filename}')
    downsample_coord_filename = get_filename_with_downsample_suffix(coord_filename, out_ratio)
    line_count = get_line_count(coord_filename)
    out_indices = np.random.choice(line_count, int(line_count * out_ratio), replace=False)
    with open(coord_filename, 'r') as file, open(downsample_coord_filename, 'w') as outfile:
        for line_index, line in enumerate(file):
            if line_index in out_indices:
                outfile.write(line)
                
    print(f'-> Downsampled coordinates saved to {downsample_coord_filename}')
    
    return downsample_coord_filename
    

def transport_files(target_dirname, downsample_ratio):
    """
    :param downsample_ratio: the ratio of the output to the total number of samples
    :param target_dirname: an abosolute path
    """
    purpose_list = ['train', 'test']
    label_list = ['pos', 'neg']
    should_downsample = downsample_ratio < 1

    os.chdir('uw_gm12878_ctcf')
    for purpose in purpose_list:
        for label in label_list:
            source_sub_dirname = f'uw_gm12878_ctcf.{purpose}.{label}.coord.mult_species'
            
            if should_downsample:
                target_sub_dirname = os.path.join(
                    target_dirname,
                    get_filename_with_downsample_suffix(source_sub_dirname, downsample_ratio))
                
                # We will save the downsampled data to source_downsample_sub_dirname before copying them
                # to the target subdirectory target_sub_dirname
                source_downsample_sub_dirname = get_filename_with_downsample_suffix(source_sub_dirname, downsample_ratio)
                os.makedirs(source_downsample_sub_dirname, exist_ok=False)
                downsampled_coord_filename = generate_downsample_coord(f'uw_gm12878_ctcf.{purpose}.{label}.coord', downsample_ratio)
            else:
                target_sub_dirname = os.path.join(target_dirname, source_sub_dirname)
                
            print(f'\n=> Creating target subdirectory with target_sub_dirname: {target_sub_dirname}')
            os.makedirs(target_sub_dirname, exist_ok=False)
            
            source_subdir_files = os.listdir(source_sub_dirname)
            print(f'=> Source subdirectory: {source_sub_dirname}')
            
            for filename in source_subdir_files:
                if filename.endswith('.fa'):
                    target_filepath = os.path.join(target_sub_dirname, create_target_filename(filename, purpose, label))
                    source_filepath = os.path.join(source_sub_dirname, filename)
                    
                    if should_downsample:
                        source_downsampled_filepath = os.path.join(source_downsample_sub_dirname,
                                                                   get_filename_with_downsample_suffix(filename, downsample_ratio))
                        downsample(source_filepath, source_downsampled_filepath, downsampled_coord_filename)
                        
                        print(f'=> Copying {source_downsampled_filepath} to {target_filepath}')
                        shutil.copyfile(source_downsampled_filepath, target_filepath)
                        
                    else:
                        print(f'=> Copying {source_filepath} to {target_filepath}')
                        shutil.copyfile(source_filepath, target_filepath)
                    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('target_dirname')
    parser.add_argument('-d', '--downsample', type=float, default=1)
    args = parser.parse_args()
    # /Users/aaron/lsgkm/tests
    transport_files(args.target_dirname, downsample_ratio=args.downsample)
