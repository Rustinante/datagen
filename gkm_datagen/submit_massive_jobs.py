import argparse
import os
from subprocess import check_call

parser = argparse.ArgumentParser()
parser.add_argument('narrowfile_list_filename')
parser.add_argument('narrow_dirname')
parser.add_argument('line_index', type=int)
args = parser.parse_args()

filename = args.narrowfile_list_filename
narrow_dirname = args.narrow_dirname
line_index = args.line_index

print(f'=> submit massive jobs for line index {line_index}')


with open(f'{filename}', 'r') as f:
    for index, line in enumerate(f):
        if index != line_index:
            continue
        tokens = line.strip().split()
        filename = tokens[0]
        output_prefix = tokens[1]
        print(f'=> Calling ./run.sh {filename} {output_prefix}')
        with open(f'{filename}.stdout', 'w') as outfile, open(f'{filename}.stderr', 'w') as errfile:
            check_call(['./run.sh', os.path.join(narrow_dirname, filename), output_prefix], stdout=outfile, stderr=errfile)
        print(f'=> Finished running ./run.sh')
