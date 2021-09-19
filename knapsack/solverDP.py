#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
from subprocess import Popen, PIPE


def solve_it(input_data):
    # Writes the inputData to a temporay file
    tmpFileName = 'tmp.data'
    tmpFile = open(tmpFileName, 'w')
    tmpFile.write(input_data)
    tmpFile.close()


    # proc = Popen(['julia', 'greedy.jl',tmpFileName],stdout=PIPE)
    proc = Popen(['julia', 'DP.jl',tmpFileName],stdout=PIPE)
    #proc = Popen(['julia', 'BnB.jl',tmpFileName],stdout=PIPE)
    

    (out,err) = proc.communicate()
    os.remove(tmpFileName)

    return(out.decode("ascii"))


import sys

if __name__ == '__main__':
    if len(sys.argv) > 1:
        file_location = sys.argv[1].strip()
        with open(file_location, 'r') as input_data_file:
            input_data = input_data_file.read()
        print(solve_it(input_data))
    else:
        print('This test requires an input file.  Please select one from the data directory. (i.e. python solver.py ./data/ks_4_0)')




