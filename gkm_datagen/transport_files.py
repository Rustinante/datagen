import shutil
import os
import argparse
import numpy as np


def downsample(in_filename, out_filename, out_ratio):
    print(f'\n=> Downsampling...\n'
          f'in_filename: {in_filename}\n'
          f'out_filename: {out_filename}\n'
          f'out_ratio: {out_ratio}')
    with open(in_filename, 'r') as infile, open(out_filename, 'w') as outfile:
        keep = True
        for index, line in enumerate(infile):
            if index % 2 == 0:
                keep = np.random.binomial(1, out_ratio)
            if keep:
                outfile.write(line)


def create_target_filename(filename, purpose, label):
    index = filename.split('_')[0]
    return f'{index}_{purpose}_{label}.fa'


def transport_files(target_dirname, downsample_ratio):
    """
    :param target_dirname: an abosolute path
    """
    purpose_list = ['train', 'test']
    label_list = ['pos', 'neg']
    should_downsample = downsample_ratio < 1
    
    for purpose in purpose_list:
        for label in label_list:
            source_sub_dirname = f'uw_gm12878_ctcf.{purpose}.{label}.coord.mult_species'
            
            if should_downsample:
                target_sub_dirname = os.path.join(target_dirname, source_sub_dirname + f'.at{downsample_ratio:.2f}')
            else:
                target_sub_dirname = os.path.join(target_dirname, source_sub_dirname)
                
            print(f'\n=> Creating target subdirectory with target_sub_dirname: {target_sub_dirname}')
            os.makedirs(target_sub_dirname, exist_ok=False)
            
            source_subdir_files = os.listdir(source_sub_dirname)
            print(f'=> Copying from directory: {source_sub_dirname}')
            
            for filename in source_subdir_files:
                if filename.endswith('.fa'):
                    target_filepath = os.path.join(target_sub_dirname, create_target_filename(filename, purpose, label))
                    source_filepath = os.path.join(source_sub_dirname, filename)
                    
                    if should_downsample:
                        downsampled_filepath = os.path.join(source_sub_dirname, filename + f'.at{downsample_ratio:.2f}')
                        downsample(source_filepath, downsampled_filepath, 0.2)
                        
                        print(f'=> Copying {downsampled_filepath} to {target_filepath}')
                        shutil.copyfile(downsampled_filepath, target_filepath)
                        
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