import os
import argparse


def get_start_end_location_from_line(line):
    tokens = line.split()
    return int(tokens[1]), int(tokens[2])


def scan_through_line_for_number(alignment_file, start_line_hint, number):
    alignment_file.seek(start_line_hint)

    for line in alignment_file:
        start, end_exclusive = get_start_end_location_from_line(line)

        if start <= number < end_exclusive:
            return line, start_line_hint

        # The lines are sorted so if the current location is already greater than the one
        # we're searching for we know what we search does not exist.
        elif number < start:
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
        start, end_exclusive = get_start_end_location_from_line(line)
        return (line, low) if start <= number < end_exclusive else None
    
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
    start, end_exclusive = get_start_end_location_from_line(line)
    # print(low, original_mid, mid, high, location)
    if start <= number < end_exclusive:
        return line, mid
    elif start > number:
        return binary_search(low, original_mid, number, file)
    else:
        return binary_search(mid, high, number, file)


def search(file, number, filename):
    """
    Performs binary search on the lines of the chr*_segmentation file for the line containing the nucleotide coordinate number.
    Each line of the file will contain four items separated by either tabs or white spaces:
    chromosome number, start(inclusive), end(exclusive), conservation state
    e.g. chr19 0 60000 U96

    :param file: a file object that's opened in read mode
    :param number: the nucleotide coordinate number we're searching for
    :param filename: the name of the file that's opened that contains the maf sequence e.g. chr19_segmentation.bed

    :return: (line, byte-offset) or None
    
    If there is a line containing the number, the function returns a tuple whose first element is the line we
    are searching for as a string ending with a newline character, and whose second element is the byte offset
    of the first byte of the line from the beginning of the file.
    
    Otherwise returns None.
    """
    start_offset = 0
    file.seek(start_offset)
    
    # high is the byte offset of the last byte in the file.
    high = os.stat(filename).st_size - 1
    
    return binary_search(start_offset, high, number, file)


if __name__ == '__main__':
    # This is statement is used for debugging
    parser = argparse.ArgumentParser()
    # an example filename would be chr19_segmentation.bed
    parser.add_argument('filename')
    parser.add_argument('number', type=int)
    args = parser.parse_args()
    with open(args.filename, 'r') as sequence_file:
        result = search(file=sequence_file, number=args.number, filename=args.filename)
        if result:
            print(result[0].split(), result[1])
        else:
            print('Not found')
