def create_chrom_state_files():
    file_dict = {}
    # chromosome 1 to 22
    for i in range(1, 23):
        filename = f'chr{i}_phast_cons.txt'
        file_dict[f'chr{i}'] = open(filename, 'w')
    file_dict['chrX'] = open('chrX_phast_cons.txt', 'w')
    file_dict['chrY'] = open('chrY_phast_cons.txt', 'w')
    return file_dict


if __name__ == '__main__':
    chrom_file_dict = create_chrom_state_files()
    
    linecount = 0
    with open('phastConsElements100way.txt', 'r') as file:
        for _ in file:
            linecount += 1
    
    error_set = set()
    error_count = 0
    with open('phastConsElements100way.txt', 'r') as file:
        
        for index, line in enumerate(file):
            tokens = line.split()
            chrom = tokens[1]
            
            try:
                chrom_file_dict[chrom].write(line)
            except KeyError as error:
                error_count += 1
                if error not in error_set:
                    error_set.add(error)
                continue
            
            index += 1
            if index % 1000 == 0:
                print(f'[{index}/{linecount}] = {index/linecount:.2%} {chrom} #key_errors: {error_count}', end='\r')
    
    print('\n=> Closing chrom files...')
    for file in chrom_file_dict.values():
        file.close()
    
    print('=> Done')
