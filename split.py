
line_number = 0
with open ('coord', 'r') as f, open('train_coord', 'w') as train_file, open('valid_coord', 'w') as valid_file, open('test_coord', 'w') as test_file:
    for i in range(2200000):
        line = f.readline()
        train_file.write(line)
        line_number += 1
    print("{} lines for train_coord".format(line_number))
    line_number = 0
    for i in range(2200001, 2204000 + 1):
        line = f.readline()
        valid_file.write(line)
        line_number += 1
    print("{} lines for valid_coord".format(line_number))
    line_number = 0
    while True:
        line = f.readline()
        if not line:
            break
        test_file.write(line)
        line_number += 1
    print("{} lines for test_coord".format(line_number))
    
    
