#### Arguments to be passed ####
# Argument 1 - input file from netCTLpan
# Argument 2 - output file

import sys
print("Transcript", "Peptide", "MHC_binding_log50(nM)", "TAP", "Cleavage", file=open(sys.argv[2], 'a')) 
with open(sys.argv[1], 'r') as mhc:
    for line in mhc:
        words = line.split()
        #print(len(words))
        if len(words) == 39:
            number = float(words[4])
            number1 = float(words[10])
            number2 = float(words[16])
            number3 = float(words[22])
            number4 = float(words[28])
            number5 = float(words[34])
            if number < 0.58910:
                print(words[1], words[2], words[4], words[5], words[6], file=open(sys.argv[2], 'a'))
            else:
                if number1 < 0.58910:
                    print(words[1], words[2], words[10], words[11], words[12],file=open(sys.argv[2], 'a'))
                else:
                    if number2 < 0.58910:
                        print(words[1], words[2], words[16], words[17], words[18], file=open(sys.argv[2], 'a'))
                    else:
                        if number3 < 0.58910:
                            print(words[1], words[2], words[22], words[23], words[24], file=open(sys.argv[2], 'a'))
                        else:
                            if number4 < 0.58910:
                                print(words[1], words[2], words[28], words[29], words[30], file=open(sys.argv[2], 'a'))
                            else:
                                if number5 < 0.58910:
                                    print(words[1], words[2], words[34], words[35], words[36], file=open(sys.argv[2], 'a'))
        else:
            if len(words) == 33:
                number = float(words[4])
                number1 = float(words[10])
                number2 = float(words[16])
                number3 = float(words[22])
                number4 = float(words[28])
                if number < 0.58910:
                    print(words[1], words[2], words[4], words[5], words[6], file=open(sys.argv[2], 'a'))
                else:
                    if number1 < 0.58910:
                        print(words[1], words[2], words[10], words[11], words[12],file=open(sys.argv[2], 'a'))
                    else:
                        if number2 < 0.58910:
                            print(words[1], words[2], words[16], words[17], words[18], file=open(sys.argv[2], 'a'))
                        else:
                            if number3 < 0.58910:
                                print(words[1], words[2], words[22], words[23], words[24], file=open(sys.argv[2], 'a'))
                            else:
                                if number4 < 0.58910:
                                    print(words[1], words[2], words[22], words[23], words[24], file=open(sys.argv[2], 'a'))
            else:
                if len(words) == 27:
                    number = float(words[4])
                    number1 = float(words[10])
                    number2 = float(words[16])
                    number3 = float(words[22])
                    if number < 0.58910:
                        print(words[1], words[2], words[4], words[5], words[6], file=open(sys.argv[2], 'a'))
                    else:
                        if number1 < 0.58910:
                            print(words[1], words[2], words[10], words[11], words[12],file=open(sys.argv[2], 'a'))
                        else:
                            if number2 < 0.58910:
                                print(words[1], words[2], words[16], words[17], words[18], file=open(sys.argv[2], 'a'))
                            else:
                                if number3 < 0.58910:
                                    print(words[1], words[2], words[22], words[23], words[24], file=open(sys.argv[2], 'a'))
                else:
                    if len(words) == 21:
                        number = float(words[4])
                        number1 = float(words[10])
                        number2 = float(words[16])
                        if number < 0.58910:
                            print(words[1], words[2], words[4], words[5], words[6], file=open(sys.argv[2], 'a'))                
                        else:
                            if number1 < 0.58910:
                                print(words[1], words[2], words[10], words[11], words[12], file=open(sys.argv[2], 'a'))
                            else:
                                if number2 < 0.58910:
                                    print(words[1], words[2], words[16], words[17], words[18], file=open(sys.argv[2], 'a'))
                    else:
                        print("Number of HLA alleles not currently accounted for, update program")
