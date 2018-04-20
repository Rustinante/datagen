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
    
    with open(in_filename, 'r') as infile, open(out_filename, 'w') as outfile:
        for index, line in enumerate(infile):
            triple = map_line_to_coord_triple(line)
            if triple in coord_set:
                chrom, start, stop = triple
                outfile.write(f'>{chrom}:{start}-{stop}\n')
                outfile.write(infile.readline())
            else:
                # Discard the current line and the next line
                infile.readline()


def create_target_filename(filename, purpose, label):
    index = filename.split('_')[0]
    return f'{index}_{purpose}_{label}.fa'


def get_filename_with_downsample_suffix(filename, downsample_ratio):
    # return f'{filename}.at{downsample_ratio:.2f}'
    return f'{filename}.at1.00'


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


def determine_downsample_ratio(coord_filename):
    max_samples = 15000
    num_samples = get_line_count(coord_filename)
    if num_samples <= max_samples:
        return 1.
    else:
        return max_samples / num_samples
    

def transport_files(source_dirname, target_dirname):
    """
    :param source_dirname: can be a relative path or an absolute path
    :param target_dirname: an absolute path
    """
    purpose_list = ['train', 'test']
    label_list = ['pos', 'neg']
    
    original_dirname = os.getcwd()
    print(f'=> Changing the current working directory to {source_dirname}')
    os.chdir(source_dirname)
    src_dir_basename = os.path.basename(os.path.normpath(source_dirname))
    
    for purpose in purpose_list:
        for label in label_list:
            source_sub_dirname = f'{src_dir_basename}.{purpose}.{label}.coord.mult_species'
            
            coord_filename = f'{src_dir_basename}.{purpose}.{label}.coord'
            downsample_ratio = determine_downsample_ratio(coord_filename)
            print(f'-> Downsample ratio for {coord_filename}: {downsample_ratio:.2f}')
            
            target_sub_dirname = os.path.join(
                target_dirname,
                get_filename_with_downsample_suffix(source_sub_dirname, downsample_ratio))
            
            # We will save the downsampled data to source_downsample_sub_dirname before copying them
            # to the target subdirectory target_sub_dirname
            source_downsample_sub_dirname = get_filename_with_downsample_suffix(source_sub_dirname, downsample_ratio)
            os.makedirs(source_downsample_sub_dirname, exist_ok=False)
            downsampled_coord_filename = generate_downsample_coord(coord_filename, downsample_ratio)
            
            print(f'\n=> Creating target subdirectory with target_sub_dirname: {target_sub_dirname}')
            os.makedirs(target_sub_dirname, exist_ok=False)
            
            shutil.copyfile(downsampled_coord_filename, os.path.join(target_sub_dirname, downsampled_coord_filename))
            
            source_subdir_files = os.listdir(source_sub_dirname)
            print(f'=> Source subdirectory: {source_sub_dirname}')
            
            for filename in source_subdir_files:
                if filename.endswith('.fa.ir'):
                    suffix_len = len('.ir')
                    target_filepath = os.path.join(target_sub_dirname,
                                                   create_target_filename(filename[:-suffix_len - 1], purpose, label))
                    source_filepath = os.path.join(source_sub_dirname, filename)
                    
                    source_downsampled_filepath = os.path.join(source_downsample_sub_dirname,
                                                               get_filename_with_downsample_suffix(filename, downsample_ratio))
                    downsample(source_filepath, source_downsampled_filepath, downsampled_coord_filename)
                    
                    print(f'=> Copying {source_downsampled_filepath} to {target_filepath}')
                    shutil.copyfile(source_downsampled_filepath, target_filepath)
                    
    print(f'=> Changing the working directory back to {original_dirname}')
    os.chdir(original_dirname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        usage='To be called directly with python3 transport_files.py and not as a module with the -m option.')
    parser.add_argument('source_dirname', help='The top level directory to copy from')
    parser.add_argument('target_dirname', help='The target directory to copy data to.')
    args = parser.parse_args()
    if args.target_dirname[0] != '/':
        print(f'target dirname must be an absolute path')
        exit(1)
    # /Users/aaron/lsgkm/tests
    transport_files(args.source_dirname, args.target_dirname)
