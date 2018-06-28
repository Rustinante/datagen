import os
import argparse


def get_location_from_line(line):
    if not line:
        return None
    return int(line.split(',')[0])


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


def binary_search(low, high, number, file):
    if low > high:
        return None
    
    # There is only one potential candidate.
    if low == high:
        # Always make sure to seek to the actual offset that we want to examine.
        file.seek(low)
        line = file.readline()
        location = get_location_from_line(line)
        if not location:
            return None
        return (line, low) if location == number else None
    
    # There are at least two bytes to work with.
    mid = low + (high - low) // 2
    original_mid = mid
    file.seek(mid)
    
    # Move right until the next new line character
    while file.read(1) != '\n':
        mid += 1
        # However, if we hit the upper bound before hitting a newline character,
        # we'd know that the line we're searching for could only be in the left half.
        if mid == high:
            return binary_search(low, original_mid, number, file)
    # When we hit the newline character, the mid would not be incremented in the while loop body,
    # so we have to increment it here to keep it synchronized with the position of the file pointer.
    mid += 1
    
    line = file.readline()
    location = get_location_from_line(line)
    # print(low, original_mid, mid, high, location)
    if not location:
        if mid < high:
            raise ValueError("mid < high but no location can be obtained")
        return binary_search(low, original_mid, number, file)
    elif location == number:
        return line, mid
    elif location > number:
        return binary_search(low, original_mid, number, file)
    else:
        return binary_search(mid, high, number, file)


def search(file, number, file_byte_size):
    """
    Performs binary search on the lines of the file for the line containing the nucleotide coordinate number.

    WARNING: Assumes there is one header line that ends with a newline.
    WARNING: Assumes all other lines start with a coordinate number followed by a comma.

    :param file: a file object that's opened in read mode
    :param number: the nucleotide coordinate number we're searching for
    :param file_byte_size: can be obtained by calling os.stat(filename).st_size
    an example filename is chr2_maf_sequence.csv

    :return: (line, byte-offset) or None
    If there is a line containing the number, the function returns a tuple whose first element is the line we
    are searching for as a string ending with a newline character, and whose second element is the byte offset
    of the first byte of the line from the beginning of the file.

    It returns None if there is no line containing the number.
    """
    file.seek(0)
    
    # high is the byte offset of the last byte in the file.
    high = file_byte_size - 1
    
    # Getting rid of the header line.
    start_index = 0
    while file.read(1) != '\n':
        start_index += 1
        if start_index >= high:
            print('The header line extends onto or beyond the upper bound.')
            return None
    
    # Start at the position after the newline.
    start_index += 1
    
    # print('binary search in range {} to {}'.format(start_index, high))
    return binary_search(start_index, high, number, file)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # an example filename would be chr2_maf_sequence.csv
    parser.add_argument('filename')
    parser.add_argument('number')
    args = parser.parse_args()
    
    with open(args.filename, 'r') as sequence_file:
        print(search(file=sequence_file, number=int(args.number), file_byte_size=os.stat(args.filename).st_size))
