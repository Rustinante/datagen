import os
import twobitreader

flanking_bp = 400
center_bp = 200
total_bp = flanking_bp * 2 + center_bp


def get_destination_filename(chrom, data_purpose):
    return os.path.join('data', f'{chrom}_{data_purpose}')


def coord_to_letter(coord_filename, data_purpose, genome_dict):
    print(f'Handling coordinate file: {coord_filename}')
    opened_files = {}
    with open(coord_filename, 'r') as coord_file:
        for line in coord_file:
            tokens = line.split()
            chrom, start, stop = tokens[0], int(tokens[1]), int(tokens[2])
            if chrom not in opened_files:
                print(f'opening {get_destination_filename(chrom, data_purpose)}')
                opened_files[chrom] = open(get_destination_filename(chrom, data_purpose), 'w')

            dna_sequence = genome_dict[chrom].get_slice(start - flanking_bp, stop + flanking_bp)
            if len(dna_sequence) != total_bp:
                print(f'The DNA sequence is not {total_bp} bp in length!')
                print(f'chr: {chrom} start: {start} stop: {stop}')

                middle = genome_dict[chrom].get_slice(start, stop)
                assert (len(middle) == center_bp)

                left = genome_dict[chrom].get_slice(start - flanking_bp, start)
                right = genome_dict[chrom].get_slice(stop + 1, stop + flanking_bp + 1)

                if len(left) != flanking_bp:
                    print(f'The left flanking sequence has length {len(left)}')
                    padding = 'N' * (flanking_bp - len(left))
                    left = padding + left

                if len(right) != flanking_bp:
                    print(f'The right flanking sequence has length {len(right)}')
                    padding = 'N' * (flanking_bp - len(right))
                    right = right + padding

                padded_sequence = left + middle + right
                assert (len(padded_sequence) == total_bp)
                print(f'Padded the sequence to {padded_sequence}')
                opened_files[chrom].write(f'{start},{padded_sequence}\n')
                continue

            opened_files[chrom].write(f'{start},{dna_sequence}\n')

    for chr_index, file in opened_files.items():
        print(f'closing {get_destination_filename(chr_index, data_purpose)}')
        file.close()


def run():
    genome_dict = twobitreader.TwoBitFile('hg19.2bit')
    os.makedirs('data', exist_ok=False)
    coord_to_letter('train_coord', 'train', genome_dict=genome_dict)
    coord_to_letter('valid_coord', 'valid', genome_dict)
    coord_to_letter('test_coord', 'test', genome_dict)


if __name__ == '__main__':
    run()
