import argparse
from argparse import RawTextHelpFormatter
import os
import sys
import subprocess
import multiprocessing
from multiprocessing import Manager
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
import pickle
import shutil
import re
import time
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
=================================================
 Classify exons based on output from rsrBulk.py
=================================================
""")
parser.add_argument('--exons', action='store',
                    required=True,
                    help='Input list of exons')
parser.add_argument('--p', action='store',
                    type=int,
                    required=True,
                    help='Number of processors')
parser.add_argument('--cfg', action='store',
                    required=True,
                    help='Config file preset')
parser.add_argument('--gtf', action='store',
                    required=True,
                    help='GTF file')
parser.add_argument('--out', action='store',
                    default='AltExonData.tab',
                    help='Output exon datatable filename')
parser.add_argument('--min', action='store',
                    type=int,
                    default=30,
                    help='Minimum junction count for LPSI/RPSI')
parser.add_argument('--new', action='store',
                    type=int,
                    default=0,
                    help='Make new exon list')
parser.add_argument('--deldir', action='store',
                    type=int,
                    default=1,
                    help='Delete output directory')
args = parser.parse_args()


def worker1(sdir, ddir, LPSI, RPSI):
    from rsrFunctions import classifyExons
    from rsrFunctions import writeNewInput
    from rsrFunctions import writeTimestamp
    exid2loc = pickle.load(open(ddir + '/exid2loc.pickle', 'rb'))
    exid2strand = pickle.load(open(ddir + '/exid2strand.pickle', 'rb'))
    PSImax = pickle.load(open(ddir + '/PSImax.pickle', 'rb'))
    PSImin = pickle.load(open(ddir + '/PSImin.pickle', 'rb'))

    inputCount = 1
    checkpointCount = 1
    currentTime = time.time()
    classifyExons('y_input_v1.csv', 6, sdir, exid2loc, exid2strand, PSImax, PSImin, LPSI, RPSI)
    writeTimestamp(sdir, checkpointCount, currentTime)
    currentTime = time.time()
    checkpointCount = checkpointCount + 1

    dropList = []
    with open(sdir + '/df2_LASSRASS.csv', 'r') as f:
        df_main = pd.read_csv(f, sep=',', index_col=0)
    for ix in range(len(df_main.index)):
        if df_main['MergedCount'][ix] >= 2:
            dropList.append(df_main['ExonID'][ix])
        elif df_main['ExonClass'][ix].split(':')[0] == 'Filtered':
            dropList.append(df_main['ExonID'][ix])
    with open(sdir + '/df6_REFINE.csv', 'r') as f:
        df_main = pd.read_csv(f, sep=',', index_col=0)
    for ix in range(len(df_main.index)):
        LASSinfo = df_main['MetadataLASS'][ix].split(',')
        RASSinfo = df_main['MetadataRASS'][ix].split(',')
        LINKEDinfo = df_main['MetadataLINKED'][ix].split(',')
        if df_main['Cassette'][ix] == 'Yes' or df_main['MetadataMUTEX'][ix] != 'None':
            if len(LASSinfo) > 1:
                for x in LASSinfo:
                    if x in dropList:
                        dropList.remove(x)
            if len(RASSinfo) > 1:
                for x in RASSinfo:
                    if x in dropList:
                        dropList.remove(x)
        if df_main['MetadataLINKED'][ix] != 'None':
            if len(LINKEDinfo) == 0:
                if df_main['ExonID'][ix] in dropList:
                    dropList.remove(df_main['ExonID'][ix])

    writeNewInput(sdir, inputCount, dropList)
    inputCount = inputCount + 1
    writeTimestamp(sdir, checkpointCount, currentTime)
    currentTime = time.time()
    checkpointCount = checkpointCount + 1

    if os.path.getsize(sdir + '/y_input_v' + str(inputCount) + '.csv') > 0:
        classifyInput = 'y_input_v' + str(inputCount) + '.csv'
        classifyExons(classifyInput, 6, sdir, exid2loc, exid2strand, PSImax, PSImin, LPSI, RPSI)
        writeTimestamp(sdir, checkpointCount, currentTime)
        currentTime = time.time()
        checkpointCount = checkpointCount + 1

        with open(sdir + '/df6_REFINE.csv', 'r') as f:
            df_main = pd.read_csv(f, sep=',', index_col=0)
        df_main['Temp'] = ':'
        df_main['ExonLocation'] = df_main['ExonChr'].map(str) + df_main['Temp'] + \
                                  df_main['ExonStart'].astype(str) + df_main['Temp'] + \
                                  df_main['ExonEnd'].astype(str) + df_main['Temp'] + \
                                  df_main['ExonStrand'].astype(str)
        df_main = df_main.drop('Temp', axis=1)

        dropList = []
        LINKEDmatch = []
        for ix in range(len(df_main.index)):
            if df_main['ExonClass'][ix] == 'LINKED':
                if df_main['ExonLocation'][ix] not in LINKEDmatch:
                    LINKEDmatch.append(df_main['ExonLocation'][ix])
        for ix in range(len(df_main.index)):
            if df_main['ExonClass'][ix] != 'LINKED':
                if df_main['ExonLocation'][ix] in LINKEDmatch:
                    dropList.append(df_main['ExonID'][ix])
            if df_main['ExonClass'][ix] == 'Unclassified':
                if df_main['Cassette'][ix] == 'No':
                    if df_main['MetadataLASS'][ix] == 'None':
                        if df_main['MetadataRASS'][ix] == 'None':
                            if df_main['MetadataLINKED'][ix] == 'None':
                                if df_main['MetadataMUTEX'][ix] == 'None':
                                    dropList.append(df_main['ExonID'][ix])
            if df_main['Cassette'][ix] == 'Yes':
                if df_main['LPSImax1'][ix] <= 10:
                    dropList.append(df_main['ExonID'][ix])
                elif df_main['RPSImax1'][ix] <= 10:
                    dropList.append(df_main['ExonID'][ix])
        writeNewInput(sdir, inputCount, dropList)
        inputCount = inputCount + 1
        writeTimestamp(sdir, checkpointCount, currentTime)
        currentTime = time.time()
        checkpointCount = checkpointCount + 1

        if os.path.getsize(sdir + '/y_input_v' + str(inputCount) + '.csv') > 0:
            classifyInput = 'y_input_v' + str(inputCount) + '.csv'
            classifyExons(classifyInput, 6, sdir, exid2loc, exid2strand, PSImax, PSImin, LPSI, RPSI)
            writeTimestamp(sdir, checkpointCount, currentTime)
        else:
            os.remove(sdir + '/df6_REFINE.csv')


def worker2(df_main, fname, metaLASS, metaRASS, LPSI, RPSI):
    def addPSI(exID_list, leftOrRight):
        sumList = []
        if leftOrRight == 'LEFT':
            for exID in exID_list:
                sumList.append(list(map(float, LPSI[exID].split(','))))
        elif leftOrRight == 'RIGHT':
            for exID in exID_list:
                sumList.append(list(map(float, RPSI[exID].split(','))))
        else:
            sys.exit('Error in addPSI(): improper leftOrRight value')
        PSIsum = [sum(x) for x in zip(*sumList)]
        PSIsum = [round(x, 3) for x in PSIsum]
        for col in range(len(PSIsum)):
            neg1count = 0
            for row in range(len(sumList)):
                if sumList[row][col] == -1:
                    neg1count = neg1count + 1
            if neg1count < len(sumList):
                PSIsum[col] = PSIsum[col] + neg1count
            elif neg1count == len(sumList):
                PSIsum[col] = -1
            elif len(sumList) == 0:
                sys.exit('Error in addPSI(): sumList is empty')
        PSIsum = [str(x) for x in PSIsum]
        return ','.join(PSIsum)

    dropList = []
    for ix in df_main.index:
        if df_main['ExonClass'][ix] == 'LINKED':
            linkIDs = df_main['MetadataLINKED'][ix].split(',')
            exonClass = [linkIDs[0].split('_')[0], linkIDs[1].split('_')[0]]
            if exonClass[0] in ['LASS', 'RASS'] and exonClass[1] in ['LASS', 'RASS'] and exonClass[0] != exonClass[1]:
                dropList.append(ix)
            else:
                if exonClass[0] == 'LASS':
                    df_main['LeftID'][ix] = metaLASS[linkIDs[0]]
                elif exonClass[0] == 'RASS':
                    df_main['LeftID'][ix] = metaRASS[linkIDs[0]].split(',')[0]
                else:
                    df_main['LeftID'][ix] = linkIDs[0]
                if exonClass[1] == 'LASS':
                    df_main['RightID'][ix] = metaLASS[linkIDs[1]].split(',')[0]
                elif exonClass[1] == 'RASS':
                    df_main['RightID'][ix] = metaRASS[linkIDs[1]]
                else:
                    df_main['RightID'][ix] = linkIDs[1]
        elif df_main['ExonClass'][ix] in ['Cassette', 'LASS', 'RASS'] or df_main['MetadataMUTEX'][ix] != 'None':
            if df_main['MetadataLASS'][ix] != 'None' and df_main['MetadataRASS'][ix] == 'None':
                df_main['LeftID'][ix] = df_main['MetadataLASS'][ix]
                df_main['RightID'][ix] = df_main['MetadataLASS'][ix].split(',')[0]
            if df_main['MetadataLASS'][ix] == 'None' and df_main['MetadataRASS'][ix] != 'None':
                df_main['LeftID'][ix] = df_main['MetadataRASS'][ix].split(',')[0]
                df_main['RightID'][ix] = df_main['MetadataRASS'][ix]
            if df_main['MetadataLASS'][ix] == 'None' and df_main['MetadataRASS'][ix] == 'None':
                df_main['LeftID'][ix] = df_main['ExonID'][ix]
                df_main['RightID'][ix] = df_main['ExonID'][ix]
            if df_main['MetadataLASS'][ix] != 'None' and df_main['MetadataRASS'][ix] != 'None':
                sys.exit('Error in LeftID/RightID processing: ExonID is both LASS and RASS')
        else:
            dropList.append(ix)
    df_main.drop(dropList, inplace=True)

    for ix in df_main.index:
        if '|' not in str(df_main['ExonStart'][ix]) and '|' not in str(df_main['ExonEnd'][ix]):
            df_main['ExonLength'][ix] = str(int(df_main['ExonEnd'][ix]) - int(df_main['ExonStart'][ix]) + 1)
        elif '|' in str(df_main['ExonStart'][ix]):
            es = [int(x) for x in df_main['ExonStart'][ix].split('|')]
            exonLength = []
            for a in es:
                exonLength.append(str(int(df_main['ExonEnd'][ix]) - a + 1))
            df_main['ExonLength'][ix] = '|'.join(exonLength)
        elif '|' in str(df_main['ExonEnd'][ix]):
            ee = [int(x) for x in df_main['ExonEnd'][ix].split('|')]
            exonLength = []
            for a in ee:
                exonLength.append(str(a - int(df_main['ExonStart'][ix]) + 1))
            df_main['ExonLength'][ix] = '|'.join(exonLength)
        else:
            sys.exit('Error in calculating exon length: exon is both LASS and RASS')

        if df_main['LeftID'][ix] == 'None' and df_main['RightID'][ix] == 'None':
            pass
        elif df_main['LeftID'][ix] != 'None' and df_main['RightID'][ix] != 'None':
            LID = df_main['LeftID'][ix].split(',')
            RID = df_main['RightID'][ix].split(',')
            if len(LID) > 1:
                df_main['LeftPSIstr'][ix] = addPSI(LID, 'LEFT')
            elif len(LID) == 1:
                df_main['LeftPSIstr'][ix] = LPSI[LID[0]]
            else:
                sys.exit('Error: LID is empty')
            if len(RID) > 1:
                df_main['RightPSIstr'][ix] = addPSI(RID, 'RIGHT')
            elif len(RID) == 1:
                df_main['RightPSIstr'][ix] = RPSI[RID[0]]
            else:
                sys.exit('Error: RID is empty')
        else:
            sys.exit('Error in getting LeftPSIstr/RightPSIstr: an ID is missing')

        if df_main['MetadataLASS'][ix] != 'None':
            df_main['LASS'][ix] = 'Yes'
        if df_main['MetadataRASS'][ix] != 'None':
            df_main['RASS'][ix] = 'Yes'
        if df_main['MetadataLINKED'][ix] != 'None':
            df_main['LINKED'][ix] = 'Yes'
    df_main.to_csv(fname, sep=',')


if __name__ == '__main__':
    # =========================================================================================
    # Parameters
    exonList = args.exons
    exonListFile = exonList.split('/')[-1]
    configPreset = args.cfg
    numProcesses = args.p
    if numProcesses > 25:
        numProcesses = 25
    gtfPath = args.gtf
    minJC = args.min
    makeNewExonList = args.new
    cwd = os.getcwd()
    outputDir = cwd + '/tempdir_' + time.strftime('%m%d%y_%H%M%S') + '_' + exonListFile.split('.')[0]
    dictDir = outputDir + '/dict'
    bulkDir = outputDir + '/bulkOutput/'
    classifyDir = outputDir + '/classifyOutput/'
    classifyInput = bulkDir + '/y_exLR.csv'
    newExonListDir = outputDir + '/newExonList/'
    newExonListFile = newExonListDir + exonListFile.rsplit('.', 1)[0] \
                      + time.strftime('_%m%d%y_%H%M%S.') + exonListFile.rsplit('.', 1)[1]
    if not os.path.exists(exonList):
        sys.exit('Cannot find list of exons: ' + exonList)
    if not os.path.exists(gtfPath):
        sys.exit('Cannot find GTF file: ' + gtfPath)
    if not os.path.exists(newExonListDir):
        os.makedirs(newExonListDir)
    if not os.path.exists(bulkDir):
        os.makedirs(bulkDir)
    if not os.path.exists(classifyDir):
        os.makedirs(classifyDir)
    if not os.path.exists(dictDir):
        os.makedirs(dictDir)
    print1 = time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Run started')
        
    # =========================================================================================
    # Create ExonID:Location and ExonID:Strand dictionaries
    exid2loc = {}
    exid2strand = {}
    with open(exonList, 'r') as f:
        for line in f:
            if line[0] == '#':
                pass
            else:
                attr = line.split('\t')
                exid2loc[attr[0]] = attr[1]
                exid2strand[attr[0]] = attr[2].strip()
    pickle.dump(exid2loc, open(dictDir + '/exid2loc.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(exid2strand, open(dictDir + '/exid2strand.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    print2 = time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Created exid2loc and exid2strand dictionaries')

    # =========================================================================================
    # Snaptron bulk query on exon list input
    qstr = 'python3 ./bin/rsrBulk.py --preset ' + configPreset + ' --input ./input.tab --queries-per-batch 50'
    inputlist = []
    counter = 0
    batchsize = 200
    f1 = open('z_leftpsi.csv', 'w')
    f2 = open('z_rightpsi.csv', 'w')
    f3 = open('z_avgpsi.csv', 'w')
    f4 = open('z_leftjc.csv', 'w')
    f5 = open('z_rightjc.csv', 'w')
    f6 = open('z_totaljc.csv', 'w')
    f1.close()
    f2.close()
    f3.close()
    f4.close()
    f5.close()
    f6.close()
    f7 = open('y_inLR.csv', 'w')
    f8 = open('y_exLR.csv', 'w')
    f7.close()
    f8.close()

    def splitOutput(fileToSplit):
        regex = r"\d"
        leftpsi_tempfile = open('leftpsi_temp.csv', 'w')
        rightpsi_tempfile = open('rightpsi_temp.csv', 'w')
        avgpsi_tempfile = open('avgpsi_temp.csv', 'w')
        leftjc_tempfile = open('leftjc_temp.csv', 'w')
        rightjc_tempfile = open('rightjc_temp.csv', 'w')
        totaljc_tempfile = open('totaljc_temp.csv', 'w')

        with open(fileToSplit, 'r') as f:
            for line in f:
                line = line.split(',')
                if line[1] == 'LeftPSI':
                    del line[1]
                    leftpsi_tempfile.write(','.join(line))
                elif line[1] == 'RightPSI':
                    del line[1]
                    rightpsi_tempfile.write(','.join(line))
                elif line[1] == 'AvgPSI':
                    del line[1]
                    avgpsi_tempfile.write(','.join(line))
                elif line[1] == 'LeftJunctionCount':
                    del line[1]
                    leftjc_tempfile.write(','.join(line))
                elif line[1] == 'RightJunctionCount':
                    del line[1]
                    rightjc_tempfile.write(','.join(line))
                elif line[1] == 'TotalJunctionCount':
                    del line[1]
                    totaljc_tempfile.write(','.join(line))

        leftpsi_tempfile.close()
        rightpsi_tempfile.close()
        avgpsi_tempfile.close()
        leftjc_tempfile.close()
        rightjc_tempfile.close()
        totaljc_tempfile.close()
        rightpsi_tempfile = open('rightpsi_temp.csv', 'r')
        avgpsi_tempfile = open('avgpsi_temp.csv', 'r')
        leftjc_tempfile = open('leftjc_temp.csv', 'r')
        rightjc_tempfile = open('rightjc_temp.csv', 'r')
        totaljc_tempfile = open('totaljc_temp.csv', 'r')
        leftpsi_file = open('z_leftpsi.csv', 'a')
        rightpsi_file = open('z_rightpsi.csv', 'a')
        avgpsi_file = open('z_avgpsi.csv', 'a')
        leftjc_file = open('z_leftjc.csv', 'a')
        rightjc_file = open('z_rightjc.csv', 'a')
        totaljc_file = open('z_totaljc.csv', 'a')

        with open('leftpsi_temp.csv', 'r')  as leftpsi_tempfile:
            for leftpsi_templine in leftpsi_tempfile:
                leftpsi_templine = leftpsi_templine.split(',')
                rightpsi_templine = rightpsi_tempfile.readline().split(',')
                avgpsi_templine = avgpsi_tempfile.readline().split(',')
                leftjc_templine = leftjc_tempfile.readline().split(',')
                rightjc_templine = rightjc_tempfile.readline().split(',')
                totaljc_templine = totaljc_tempfile.readline().split(',')
                if re.match(regex, str(leftpsi_templine[1][0])) != None:
                    if re.match(regex, str(rightpsi_templine[1][0])) != None:
                        leftpsi_file.write(','.join(leftpsi_templine))
                        rightpsi_file.write(','.join(rightpsi_templine))
                        avgpsi_file.write(','.join(avgpsi_templine))
                        leftjc_file.write(','.join(leftjc_templine))
                        rightjc_file.write(','.join(rightjc_templine))
                        totaljc_file.write(','.join(totaljc_templine))
        rightpsi_tempfile.close()
        avgpsi_tempfile.close()
        leftjc_tempfile.close()
        rightjc_tempfile.close()
        totaljc_tempfile.close()
        os.remove('leftpsi_temp.csv')
        os.remove('rightpsi_temp.csv')
        os.remove('avgpsi_temp.csv')
        os.remove('leftjc_temp.csv')
        os.remove('rightjc_temp.csv')
        os.remove('totaljc_temp.csv')
        leftpsi_file.close()
        rightpsi_file.close()
        avgpsi_file.close()
        leftjc_file.close()
        rightjc_file.close()
        totaljc_file.close()

    with open(exonList, 'r') as f:
        for line in f:
            inputlist.append(line)
            counter = counter + 1
            if counter % batchsize == 0:
                with open('./input.tab', 'w') as exons:
                    exons.writelines(inputlist)
                inputlist = []
                subprocess.run(qstr, shell=True)
                splitOutput('./bulkoutput.csv')
                resetfile = open('./bulkoutput.csv', 'w')
                resetfile.close()
                print('Exons processed = ' + str(counter) + time.strftime(", %b %d %Y %H:%M:%S"))
        with open('./input.tab', 'w') as exons:
            exons.writelines(inputlist)
        inputlist = []
        subprocess.run(qstr, shell=True)
        splitOutput('./bulkoutput.csv')
        resetfile = open('./bulkoutput.csv', 'w')
        resetfile.close()
        print('Total exons processed = ' + str(counter) + time.strftime(", %b %d %Y %H:%M:%S"))
    os.remove('./input.tab')
    os.remove('./bulkoutput.csv')

    LPSI_filtered = open(cwd + '/a_LPSIfiltered.csv', 'w')
    RPSI_filtered = open(cwd + '/a_RPSIfiltered.csv', 'w')
    leftpsi_file = open(cwd + '/z_leftpsi.csv', 'r')
    rightpsi_file = open(cwd + '/z_rightpsi.csv', 'r')
    leftjc_file = open(cwd + '/z_leftjc.csv', 'r')
    rightjc_file = open(cwd + '/z_rightjc.csv', 'r')
    leftjc = {}
    rightjc = {}
    for line in leftpsi_file:
        LPSI_list = line.strip().split(',')
        RPSI_list = rightpsi_file.readline().strip().split(',')
        LJC_temp = leftjc_file.readline().strip().split(',', 1)
        RJC_temp = rightjc_file.readline().strip().split(',', 1)
        leftjc[LJC_temp[0]] = LJC_temp[1]
        rightjc[RJC_temp[0]] = RJC_temp[1]
        LJC_list = [LJC_temp[0]]
        LJC_list.extend(LJC_temp[1].split(','))
        RJC_list = [RJC_temp[0]]
        RJC_list.extend(RJC_temp[1].split(','))
        for ix in range(1, len(LPSI_list)):
            if int(LJC_list[ix]) < minJC:
                LPSI_list[ix] = '-1'
            if int(RJC_list[ix]) < minJC:
                RPSI_list[ix] = '-1'
        LPSI_filtered.write(','.join(LPSI_list) + '\n')
        RPSI_filtered.write(','.join(RPSI_list) + '\n')
    LPSI_filtered.close()
    RPSI_filtered.close()
    leftpsi_file.close()
    rightpsi_file.close()
    leftjc_file.close()
    rightjc_file.close()
    pickle.dump(leftjc, open(dictDir + '/leftjc.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(rightjc, open(dictDir + '/rightjc.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)
    # leftjc = pickle.load(open(dictDir + '/leftjc.pickle', 'rb'))
    # rightjc = pickle.load(open(dictDir + '/rightjc.pickle', 'rb'))
    shutil.move(cwd + '/a_LPSIfiltered.csv', bulkDir + '/a_LPSIfiltered.csv')
    shutil.move(cwd + '/a_RPSIfiltered.csv', bulkDir + '/a_RPSIfiltered.csv')
    shutil.move(cwd + '/z_leftjc.csv', bulkDir + '/z_leftjc.csv')
    shutil.move(cwd + '/z_rightjc.csv', bulkDir + '/z_rightjc.csv')
    shutil.move(cwd + '/y_exLR.csv', bulkDir + '/y_exLR.csv')
    shutil.move(cwd + '/z_totaljc.csv', bulkDir + '/z_totaljc.csv')
    shutil.move(cwd + '/z_leftpsi.csv', bulkDir + '/z_leftpsi.csv')
    shutil.move(cwd + '/z_rightpsi.csv', bulkDir + '/z_rightpsi.csv')
    shutil.move(cwd + '/z_avgpsi.csv', bulkDir + '/z_avgpsi.csv')
    shutil.move(cwd + '/y_inLR.csv', bulkDir + '/y_inLR.csv')

    with open(bulkDir + '/run_parameters.txt', 'w') as f:
        f.write('Snaptron query saved to: ' + bulkDir)
        f.write('Snaptron query completed on: ' + time.strftime('%m/%d/%Y at %H:%M:%S UTC\n'))
        f.write('exonList = ' + exonList + '\n')
        f.write('configPreset = ' + configPreset + '\n')
    with open('time.log', 'w') as f:
        f.write(print1 + '\n')
        f.write(print2 + '\n')
        f.write(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Snaptron bulk query complete\n'))
    print(print1)
    print(print2)
    print(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Snaptron bulk query complete'))
    
    # =========================================================================================
    # Classify exons from Snaptron output
    with open('time.log', 'a') as f:
        f.write(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Classify exons from Snaptron output\n'))
    print(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Classify exons from Snaptron output'))

    exid2loc = pickle.load(open(dictDir + '/exid2loc.pickle', 'rb'))
    chrList = []
    with open(classifyInput, 'r') as f:
        for line in f:
            chr = exid2loc[line.split(',')[0]].split(':')[0]
            if chr not in chrList:
                chrList.append(chr)

    sdirList = []
    for x in chrList:
        sdirList.append([classifyDir + '/' + x, x])
    for x in sdirList:
        if not os.path.exists(x[0]):
            os.makedirs(x[0])
    for x in sdirList:
        with open(x[0] + '/rsr-cleanup.log', 'w') as f:
            f.write('rsr-cleanup timeline:\n')
        with open(x[0] + '/y_input_v1.csv', 'w') as input_v1:
            with open(classifyInput, 'r') as f:
                for line in f:
                    chr = exid2loc[line.split(',')[0]].split(':')[0]
                    if chr == x[1]:
                        input_v1.write(line)

    PSImax = {}
    PSImin = {}
    with open(bulkDir + '/a_RPSIfiltered.csv', 'r') as R_File:
        with open(bulkDir + '/a_LPSIfiltered.csv', 'r') as L_File:
            for line in L_File:
                LF_line = [x for x in line.strip().split(',') if x != '']
                RF_line = [x for x in R_File.readline().strip().split(',') if x != '']
                exonID = LF_line[0]
                LPSI_list = [float(x) for x in LF_line[1:]]
                RPSI_list = [float(x) for x in RF_line[1:]]
                L5high = sorted(LPSI_list, reverse=True)[0:5]
                R5high = sorted(RPSI_list, reverse=True)[0:5]
                L5high.extend(R5high)
                PSImax[exonID] = ','.join(str(x) for x in L5high)
                L5low = [777 if x == -1 else x for x in LPSI_list]
                R5low = [777 if x == -1 else x for x in RPSI_list]
                L5low = sorted(L5low)[0:5]
                R5low = sorted(R5low)[0:5]
                L5low.extend(R5low)
                L5low = [-1 if x == 777 else x for x in L5low]
                PSImin[exonID] = ','.join(str(x) for x in L5low)
    LPSI = {}
    RPSI = {}
    with open(bulkDir + '/z_rightpsi.csv', 'r') as R_File:
        with open(bulkDir + '/z_leftpsi.csv', 'r') as L_File:
            for line in L_File:
                LF_line = [x for x in line.strip().split(',') if x != '']
                RF_line = [x for x in R_File.readline().strip().split(',') if x != '']
                exonID = LF_line[0]
                LPSI[exonID] = ','.join(LF_line[1:])
                RPSI[exonID] = ','.join(RF_line[1:])
    pickle.dump(PSImax, open(dictDir + '/PSImax.pickle', "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(PSImin, open(dictDir + '/PSImin.pickle', "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(LPSI, open(dictDir + '/LPSI.pickle', "wb"), protocol=pickle.HIGHEST_PROTOCOL)
    pickle.dump(RPSI, open(dictDir + '/RPSI.pickle', "wb"), protocol=pickle.HIGHEST_PROTOCOL)

    with open('time.log', 'a') as f:
        f.write(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Created PSImax/PSImin/LPSI/RPSI dictionaries\n'))
    print(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Created PSImax/PSImin/LPSI/RPSI dictionaries'))

    # =========================================================================================
    # First multithreading step
    manager = Manager()
    dict_1a = manager.dict()
    dict_1b = manager.dict()
    dict_1a.update(LPSI)
    dict_1b.update(RPSI)
    pool = multiprocessing.Pool(processes=numProcesses)
    for x in sdirList:
        pool.apply_async(worker1, args=(x[0], dictDir, dict_1a, dict_1b,))
    pool.close()
    pool.join()
    df6List = []
    dfcList = []

    for x in sdirList:
        df6List.append(x[0] + '/df6_REFINE.csv')
        dfcList.append(x[0] + '/df_constitutive.csv')

    with open(classifyDir + '/df_constitutive_merged.csv', 'w') as dfcmerge:
        with open(dfcList[0], 'r') as f:
            dfcmerge.write(f.readline())
        for x in dfcList:
            with open(x, 'r') as f:
                next(f)
                for line in f:
                    dfcmerge.write(line)
    with open(classifyDir + '/df_constitutive_merged.csv', 'r') as f:
        df_main = pd.read_csv(f, sep=',', index_col=0)
    df_main = df_main.reset_index(drop=True)
    df_main.to_csv(classifyDir + '/df_constitutive_merged.csv', sep=',')

    with open(classifyDir + '/df6_merged.csv', 'w') as df6merge:
        if os.path.isfile(df6List[0]):
            with open(df6List[0], 'r') as f:
                df6merge.write(f.readline())
        for x in df6List:
            if os.path.isfile(x):
                with open(x, 'r') as f:
                    next(f)
                    for line in f:
                        df6merge.write(line)
    if os.path.getsize(classifyDir + '/df6_merged.csv') == 0:
        sys.exit('Error: No alternative exons detected')

    with open(classifyDir + '/df6_merged.csv', 'r') as f:
        df_main = pd.read_csv(f, sep=',', index_col=0)
    df_main = df_main.reset_index(drop=True)
    df_main.to_csv(classifyDir + '/df6_merged.csv', sep=',')

    with open('time.log', 'a') as f:
        f.write(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || First multithreading step complete\n'))
    print(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || First multithreading step complete'))
    
    # =========================================================================================
    # Add exon metadata
    with open(classifyDir + '/df6_merged.csv', 'r') as f:
        df_main = pd.read_csv(f, sep=',', index_col=0)
    df_main = df_main.reset_index(drop=True)
    df_main['LeftID'] = 'None'
    df_main['RightID'] = 'None'
    df_main['LeftPSIstr'] = 'None'
    df_main['RightPSIstr'] = 'None'
    df_main['ExonLength'] = 'None'
    df_main['LASS'] = 'No'
    df_main['RASS'] = 'No'
    df_main['LINKED'] = 'No'

    ixList = []
    a = len(df_main.index)
    b = int(a / 5000)
    floor = 0
    fcount = 1
    fname = ''
    for x in range(b):
        ceiling = (x + 1) * 5000
        fname = 'temp' + str(fcount) + '.csv'
        ixList.append([floor, ceiling, fname])
        floor = ceiling + 1
        fcount = fcount + 1
    ixList.append([floor, a - 1, 'temp' + str(fcount) + '.csv'])
    df_main['Colon'] = ':'
    df_main['Dash'] = '-'
    df_main['ExonLocation'] = df_main['ExonChr'].map(str) + df_main['Colon'] + df_main['ExonStart'].astype(str) + \
                              df_main['Dash'] + df_main['ExonEnd'].astype(str)
    df_main['ExonBoundary'] = df_main['ExonChr'].map(str) + df_main['Colon'] + df_main['1_ELS'].astype(str) + \
                              df_main['Dash'] + df_main['1_ERE'].astype(str)
    df_main = df_main.drop(['Colon', 'Dash'], axis=1)

    metaLASS = {}
    metaRASS = {}
    for ix in range(len(df_main.index)):
        exonClass = df_main['ExonID'][ix].split('_')[0]
        if exonClass == 'LASS':
            metaLASS[df_main['ExonID'][ix]] = df_main['MetadataLASS'][ix]
        if exonClass == 'RASS':
            metaRASS[df_main['ExonID'][ix]] = df_main['MetadataRASS'][ix]
    LPSI = pickle.load(open(dictDir + '/LPSI.pickle', 'rb'))
    RPSI = pickle.load(open(dictDir + '/RPSI.pickle', 'rb'))

    with open('time.log', 'a') as f:
        f.write(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Created metaLASS/metaRASS dictionaries\n'))
    print(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Created metaLASS/metaRASS dictionaries'))
    
    # =========================================================================================
    # Second multithreading step
    dict_2a = manager.dict()
    dict_2b = manager.dict()
    dict_2c = manager.dict()
    dict_2d = manager.dict()
    dict_2a.update(metaLASS)
    dict_2b.update(metaRASS)
    dict_2c.update(LPSI)
    dict_2d.update(RPSI)
    pool = multiprocessing.Pool(processes=numProcesses)
    for x in ixList:
        dfsubset = df_main.iloc[(df_main.index >= x[0]) & (df_main.index <= x[1])]
        pool.apply_async(worker2, args=(dfsubset, x[2], dict_2a, dict_2b, dict_2c, dict_2d,))
    pool.close()
    pool.join()

    with open('df7temp.csv', 'w') as df7temp:
        with open(ixList[0][2], 'r') as header:
            df7temp.write(header.readline())
        for x in ixList:
            with open(x[2], 'r') as f:
                next(f)
                for line in f:
                    df7temp.write(line)
            os.remove(x[2])
    with open('df7temp.csv', 'r') as f:
        df_main = pd.read_csv(f, sep=',', index_col=0)
    os.remove('df7temp.csv')
    df_main = df_main.reset_index(drop=True)
    df_main.rename(columns={'MetadataMUTEX': 'MUTEX', 'Cassette': 'CASSETTE'}, inplace=True)

    with open('time.log', 'a') as f:
        f.write(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Second multithreading step complete\n'))
    print(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Second multithreading step complete'))

    # =========================================================================================
    # Generate new exon list
    if makeNewExonList == 1:
        keepList = []
        a = df_main['LeftID'].tolist()
        a.extend(df_main['RightID'].tolist())
        for idStr in a:
            idList = idStr.split(',')
            if len(idList) == 0:
                keepList.append(idList[0])
            else:
                for id in idList:
                    keepList.append(id)
        keepList = list(set(keepList))
        dfOldInput = pd.read_table(exonList, sep='\t', index_col=0, names=['Location', 'Strand'])
        dfNewInput = dfOldInput.ix[keepList]
        dfNewInput.sort_index(inplace=True)
        dfNewInput.to_csv(newExonListFile, sep='\t', header=False)

        with open('time.log', 'a') as f:
            f.write(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || New exon list generated\n'))
        print(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || New exon list generated'))

    # =========================================================================================
    # Create genesym2attr and boundary2gene dictionaries
    parseGTF = 1
    if parseGTF == 1:
        geneBed = open('geneLocation.bed', 'w')
        name2id = {}
        with open(gtfPath, 'r') as gtf:
            count = 1
            for line in gtf:
                if line[0:2] == '##':
                    pass
                else:
                    elem = line.split('\t')
                    if elem[2] == 'gene':
                        chr = elem[0]
                        start = elem[3]
                        end = elem[4]
                        strand = elem[6]
                        gene_id = ''
                        gene_type = ''
                        gene_name = ''
                        attributes = re.split(r'[ ;]', elem[8].replace('"', ''))
                        for ix in range(len(attributes)):
                            if attributes[ix] == 'gene_id':
                                gene_id = attributes[ix + 1]
                            if attributes[ix] == 'gene_type':
                                gene_type = attributes[ix + 1]
                            elif attributes[ix] == 'gene_name':
                                gene_name = attributes[ix + 1]
                        geneBed.write('\t'.join([chr, start, end, gene_id, '1', strand]) + '\n')
                        name2id[gene_id] = gene_name + ',' + gene_type
        geneBed.close()
        pickle.dump(name2id, open(dictDir + '/genesym2attr.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    df_main['TempScore'] = '1'
    boundaryBed = df_main[['ExonChr', '1_ELS', '1_ERE', 'ExonID', 'TempScore', 'ExonStrand']]
    boundaryBed = boundaryBed.sort_values(['1_ELS', 'ExonChr'], ascending=[True, True])
    boundaryBed.to_csv('df7b2b.bed', sep='\t', header=False, index=False)
    subprocess.run('bedtools sort -i df7b2b.bed > df7b2b_sorted.bed', shell=True)
    subprocess.run('bedtools sort -i geneLocation.bed > geneLocation_sorted.bed', shell=True)
    subprocess.run('bedtools closest -s -d -k 10 -a df7b2b_sorted.bed \
                                                 -b geneLocation_sorted.bed > closestGene.bed', shell=True)
    boundary2gene = {}
    with open('closestGene.bed', 'r') as f:
        exID_track = ''
        dist_track = -1
        for line in f:
            a = line.split('\t')
            chr = a[0]
            exID = a[3]
            boundaryStart = int(a[1])
            boundaryEnd = int(a[2])
            geneID = a[9]
            dist = int(a[12])
            geneType = name2id[geneID].split(',')[1]
            if exID != exID_track:
                exID_track = exID
                dist_track = -1
                if dist < 10000:
                    dist_track = dist
                    boundary2gene[chr + ':' + str(boundaryStart) + '-' + str(boundaryEnd)] = geneID
            else:
                if geneType == 'protein_coding' and dist <= dist_track:
                    dist_track = dist
                    boundary2gene[chr + ':' + str(boundaryStart) + '-' + str(boundaryEnd)] = geneID
    pickle.dump(boundary2gene, open(dictDir + '/boundary2gene.pickle', 'wb'), protocol=pickle.HIGHEST_PROTOCOL)

    with open('time.log', 'a') as f:
        f.write(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Create genesym2attr and boundary2gene dictionaries\n'))
    print(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Create genesym2attr and boundary2gene dictionaries'))

    # =========================================================================================    
    # Add gene attribute columns
    boundary2gene = pickle.load(open(dictDir + '/boundary2gene.pickle', 'rb'))
    genesym2attr = pickle.load(open(dictDir + '/genesym2attr.pickle', 'rb'))
    df_main['GeneID'] = df_main['ExonBoundary'].map(boundary2gene)
    df_main['GeneAttributes'] = df_main['GeneID'].map(genesym2attr)
    df_main['GeneID'] = df_main['GeneID'].replace(np.nan, 'N/A')
    df_main['GeneAttributes'] = df_main['GeneAttributes'].replace(np.nan, 'N/A,N/A')
    df_main[['GeneSymbol', 'GeneType']] = df_main['GeneAttributes'].str.split(',', expand=True)
    df_main = df_main.drop(['GeneAttributes'], axis=1)
    df_main['LJCID'] = df_main['LeftID'].str.split(',', expand=True, n=1)[0]
    df_main['RJCID'] = df_main['RightID'].str.split(',', expand=True, n=1)[0]
    df_main['LeftJCstr'] = df_main['LJCID'].map(leftjc)
    df_main['RightJCstr'] = df_main['RJCID'].map(rightjc)
    os.remove('geneLocation.bed')
    os.remove('df7b2b.bed')
    os.remove('closestGene.bed')
    os.remove('geneLocation_sorted.bed')
    os.remove('df7b2b_sorted.bed')
    df_main = df_main[['ExonID', 'LeftID', 'RightID', 'CASSETTE', 'LASS', 'RASS', 'LINKED', 'MUTEX',
                       'GeneID', 'GeneSymbol', 'GeneType', 'ExonLocation', 'ExonBoundary', 'ExonStrand', 'ExonLength',
                       'LeftJCstr', 'RightJCstr', 'LeftPSIstr', 'RightPSIstr']]

    dropList = []
    for ix in df_main.index:
        if df_main['CASSETTE'][ix] == 'No' and df_main['LINKED'][ix] == 'No' and df_main['MUTEX'][ix] == 'None':
            dropList.append(ix)
    df_main.drop(dropList, inplace=True)
    df_main['MUTEX'] = df_main['MUTEX'].replace('None', 'No')
    df_main.loc[df_main['LASS'] == 'Yes', 'CASSETTE'] = 'No'
    df_main.loc[df_main['RASS'] == 'Yes', 'CASSETTE'] = 'No'
    df_main.to_csv(args.out, sep='\t', index=False)

    with open('time.log', 'a') as f:
        f.write(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Run Complete\n'))
    print(time.strftime('%-m/%-d/%Y || %-I:%M%p UTC || Run Complete'))
    shutil.move(cwd + '/time.log', outputDir + '/time.log')

    # Delete temporary output directory
    if args.deldir == 1:
        shutil.rmtree(outputDir)