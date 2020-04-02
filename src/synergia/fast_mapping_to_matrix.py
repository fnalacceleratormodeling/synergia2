#!/usr/bin/env python

import numpy

def fast_mapping_to_matrix(filename):
    f = open(filename, 'r')
    A = numpy.zeros((6, 6), 'd')
    b = numpy.zeros((6,), 'd')
    f.readline()
    f.readline()
    f.readline()
    order = int(f.readline())
    print("fast_mapping of order", order)
    f.readline()
    expected_length = int(f.readline())
    length = 0
    line = f.readline()
    while line != "end_fast_mapping\n":
        index = int(line.split('=')[1].split(',')[0])
        order = int(line.split('=')[2])
        f.readline()
        num = int(f.readline())
        if num > 0:
            f.readline()
        for i in range(0, num):
            pieces = f.readline().split()
            if order == 0:
                b[index] = float(pieces[2])
            elif order == 1:
                j = int(pieces[1])
                A[index, j] = float(pieces[3])
        line = f.readline()
        length += 1
    if length == expected_length:
        print("found correct number of terms")
    return b,A


    f.close()

if __name__ == '__main__':
    import sys
    b,A = fast_mapping_to_matrix(sys.argv[1])
    print(b)
    print(A)

