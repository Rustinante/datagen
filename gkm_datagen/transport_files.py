import shutil
import os
import argparse


def create_target_filename(filename, purpose, label):
    index = filename.split('_')[0]
    return f'{index}_{purpose}_{label}.fa'


def transport_files(target_dirname):
    """
    :param target_dirname: an abosolute path
    """
    purpose_list = ['train', 'test']
    label_list = ['pos', 'neg']
    for purpose in purpose_list:
        for label in label_list:
            source_sub_dirname = f'uw_gm12878_ctcf.{purpose}.{label}.coord.mult_species'
            
            target_sub_dirname = os.path.join(target_dirname, source_sub_dirname)
            print(f'\n=> Creating target subdirectory with target_sub_dirname: {target_sub_dirname}')
            os.makedirs(target_sub_dirname, exist_ok=False)
            
            source_subdir_files = os.listdir(source_sub_dirname)
            print(f'=> Copying from directory: {source_sub_dirname}')
            
            for filename in source_subdir_files:
                if filename.endswith('.fa'):
                    target_filepath = os.path.join(target_sub_dirname, create_target_filename(filename, purpose, label))
                    source_filepath = os.path.join(source_sub_dirname, filename)
                    print(f'=> Copying {source_filepath} to {target_filepath}')
                    shutil.copyfile(source_filepath, target_filepath)
    

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('target_dirname')
    args = parser.parse_args()
    # /Users/aaron/lsgkm/tests
    transport_files(args.target_dirname)