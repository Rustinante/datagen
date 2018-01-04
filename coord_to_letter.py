import os
import twobitreader


def get_destination_filename(chr, data_purpose):
    return os.path.join('data', '{}_{}'.format(chr, data_purpose))


def coord_to_letter(coord_filename, data_purpose, genome_dict):
    print('Handling coordinate file: {}'.format(coord_filename))
    opened_files = {}
    with open(coord_filename, 'r') as coordinate_file:
        for line in coordinate_file:
            tokens = line.split()
            chr, start, stop = tokens[0], int(tokens[1]), int(tokens[2])
            if chr not in opened_files:
                print('opening {}'.format(get_destination_filename(chr, data_purpose)))
                opened_files[chr] = open(get_destination_filename(chr, data_purpose), 'w')

            dna_sequence = genome_dict[chr].get_slice(start - 400, stop + 400)
            if len(dna_sequence) != 1000:
                print('The DNA sequence is not 200 bp in length!')
                print('chr: {} start: {} stop: {}'.format(chr, start, stop))
                continue
            opened_files[chr].write('{},{}'.format(start, dna_sequence) + '\n')

    for chr_index, file in opened_files.items():
        print('closing {}'.format(get_destination_filename(chr_index, data_purpose)))
        file.close()


def run():
    genome_dict = twobitreader.TwoBitFile('hg19.2bit')
    os.makedirs('data', exist_ok=False)
    coord_to_letter('train_coord', 'train', genome_dict=genome_dict)
    coord_to_letter('valid_coord', 'valid', genome_dict)
    coord_to_letter('test_coord', 'test', genome_dict)


if __name__ == '__main__':
    run()
