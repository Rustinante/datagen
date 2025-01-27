import os
import twobitreader

dirname = 'hg19_fasta_aggregate'


def get_destination_filename(chr, data_purpose):
    return os.path.join(dirname, '{}_{}'.format(chr, data_purpose))


def coord_to_letter(coord_filename, data_purpose, genome_dict):
    print('Handling coordinate file: {}'.format(coord_filename))
    opened_files = {}
    flanking_number = 50
    with open(coord_filename, 'r') as coordinate_file:
        for line in coordinate_file:
            tokens = line.split()
            chr, start, stop = tokens[0], int(tokens[1]), int(tokens[2])
            if chr not in opened_files:
                print('opening {}'.format(get_destination_filename(chr, data_purpose)))
                opened_files[chr] = open(get_destination_filename(chr, data_purpose), 'w')

            target_length = 2 * flanking_number + stop - start
            dna_sequence = genome_dict[chr].get_slice(start - flanking_number, stop + flanking_number)
            if len(dna_sequence) != target_length:
                print('The DNA sequence is not 200 bp in length!')
                print('chr: {} start: {} stop: {}'.format(chr, start, stop))
                
                middle = genome_dict[chr].get_slice(start, stop)
                assert (len(middle) == 200)
                
                left = genome_dict[chr].get_slice(start - flanking_number, start)
                right = genome_dict[chr].get_slice(stop + 1, stop + 401)
                
                if len(left) != flanking_number:
                    print('The left flanking sequence has length {}'.format(len(left)))
                    padding = 'N' * (flanking_number - len(left))
                    left = padding + left
                
                if len(right) != flanking_number:
                    print('The right flanking sequence has length {}'.format(len(right)))
                    padding = 'N' * (flanking_number - len(right))
                    right = right + padding
                
                padded_sequence = left + middle + right
                assert (len(padded_sequence) == target_length)
                print('Padded the sequence to {}'.format(padded_sequence))
                opened_files[chr].write('>{}\n{}'.format(start, padded_sequence) + '\n')
                continue
            
            opened_files[chr].write('>{}\n{}'.format(start, dna_sequence) + '\n')
    
    for chr_index, file in opened_files.items():
        print('closing {}'.format(get_destination_filename(chr_index, data_purpose)))
        file.close()


def run():
    genome_dict = twobitreader.TwoBitFile('hg19.2bit')
    os.makedirs(dirname, exist_ok=False)
    coord_to_letter('train_coord', 'train', genome_dict=genome_dict)
    coord_to_letter('valid_coord', 'valid', genome_dict)
    coord_to_letter('test_coord', 'test', genome_dict)


if __name__ == '__main__':
    run()
