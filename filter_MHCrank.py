#### Arguments to be passed ####
# Argument 1 - input file from netCTLpan
# Argument 2 - output file

import sys

with open(sys.argv[1], 'r') as mhc:
    for line in mhc:
        words = line.split()
        number = float(words[8])
        number1 = float(words[14])
        number2 = float(words[20])
        number3 = float(words[26])
        if number < 1.0:
            print(words[1], words[2], words[8],file=open(sys.argv[2], 'a'))
        else:
            if number1 < 1.0:
                print(words[1], words[2], words[14],file=open(sys.argv[2], 'a'))
            else:
                if number2 < 1.0:
                    print(words[1], words[2], words[20],file=open(sys.argv[2], 'a'))
                else:
                    if number3 < 1.0:
                        print(words[1], words[2], words[26],file=open(sys.argv[2], 'a'))
