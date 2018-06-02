import os
import argparse


def get_start_end_location_from_line(line):
    # returns start, end_exclusive
    tokens = line.split()
    return int(tokens[2]), int(tokens[3])


def scan_through_line_for_number(constrained_element_file, start_byteoffset_hint, number):
    constrained_element_file.seek(start_byteoffset_hint)
    
    for line in constrained_element_file:
        start, end_exclusive = get_start_end_location_from_line(line)
        
        if start <= number < end_exclusive:
            return True, line, start_byteoffset_hint
        
        # The lines are sorted so if the current location is already greater than the one
        # we're searching for we know what we search does not exist.
        # We also know that the start coordinate in line at the hint is <= the target coordinate,
        # so the moment number < start we know the line is the first line that has a start coordinate
        # greater than the target coordinate.
        elif number < start:
            return False, line, start_byteoffset_hint
        
        start_byteoffset_hint += len(bytes(line, 'ascii'))
    
    return False, None, None


def binary_search(low, high, number, file):
    if low > high:
        raise ValueError(f'The low argument ({low}) is larger than the high argument ({high})')
    
    # There is only one potential candidate.
    if low == high:
        # Always make sure to seek to the actual offset that we want to examine.
        file.seek(low)
        line = file.readline()
        start, end_exclusive = get_start_end_location_from_line(line)
        if start <= number < end_exclusive:
            return True, line, low
        elif number < start:
            # The current line must be the first line that contains a start coordinate higher than the target coordinate,
            # because even if there is a line L preceding the current line, the end coordinate in L must be less than the target
            # coordinate due to the nature of binary search.
            return False, line, low
        else:
            # number >= end_exclusive
            # Check if there is a line after the current line.
            # If there is none, then there is no line containing a start coordinate greater than the target coordiante.
            # If there is a line L after the current line, the start coordinate in L must be larger than the target coordinate.
            # So the line L is the first line that contiains a start coordinate higher than the target number.
            next_line = file.readline()
            if not next_line or next_line == '\n':
                return False, None, None
            else:
                return False, next_line, low + len(bytes(line, 'ascii'))
        
        # return (True, line, low) if start <= number < end_exclusive else None
    
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
        return True, line, mid
    elif start > number:
        return binary_search(low, original_mid, number, file)
    else:
        return binary_search(mid, high, number, file)


def search(file, number, file_byte_size):
    """
    Performs binary search on the lines of the file phastConsElements100way.txt for the line containing the
    nucleotide coordinate number.
    Each line of the file looks like the following:
    585	chr1	11991	11995	lod=12	240

    :param file: a file object that's opened in read mode
    :param number: the nucleotide coordinate number we're searching for
    :param file_byte_size: can be obtained by calling os.stat(filename).st_size

    :return: (line, byte-offset) or None

    If there is a line containing the number, the function returns a tuple whose first element is the line we
    are searching for as a string ending with a newline character, and whose second element is the byte offset
    of the first byte of the line from the beginning of the file.

    Otherwise returns None.
    """
    start_offset = 0
    file.seek(start_offset)
    
    # high is the byte offset of the last byte in the file.
    high = file_byte_size - 1
    
    return binary_search(start_offset, high, number, file)


if __name__ == '__main__':
    # This is used for debugging
    parser = argparse.ArgumentParser()
    # an example filename is chr1_phast_cons.txt
    parser.add_argument('filename')
    parser.add_argument('number', type=int)
    args = parser.parse_args()
    
    byte_size = os.stat(args.filename).st_size
    with open(args.filename, 'r') as sequence_file:
        result = search(file=sequence_file, number=args.number, file_byte_size=byte_size)
        if result:
            print(result[0].split(), result[1])
        else:
            print('Not found')
