# -*- coding: utf-8 -*-
"""
Created on Thu Aug 10 13:59:32 2017

@author: Peter Nachtwey
Delta Computer Systems, Inc.
"""

import csv
import numpy as np

def readCSV(path, ch):
    """ Reads the CSV file designated by the path.
        ch is the deliminator. It could be a comma, tab or spaces 
        The file must have a 1 line header that is ignored.
        The data must be in three columns of time, control and actual """
    lTime = []                              # define empty lists
    lAct = []
    lCtrl = []
    with open(path, newline='') as csvfile:
        readCSV = csv.reader(csvfile, delimiter=ch)
        header = next(readCSV, None)    # the header is ignored
        for row in readCSV:
            # row = [time, control, actual,.......]
            lTime.append(float(row[0]))
            lAct.append(float(row[2]))
            lCtrl.append(float(row[1]))
    csvfile.close()
    aTime = np.array(lTime)             # convert from lists to np.array
    aAct = np.array(lAct)
    aControl = np.array(lCtrl)
    return [aTime, aControl, aAct]      # return a list of three np.arrays


if __name__ == "__main__":
    """ Run this file To test the readCSV program by itself.
        The path must point to the file you want to read.
        The delminantor character must be set correctly.
        Type the arrays to verify the data is read correctly """
    path="c:\\Users\peter\Downloads\System RS Mode 1.csv"
    aTime, aControl, aAct = readCSV(path,',')
        