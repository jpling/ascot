def classifyExons(argInfile, argSteps, sdir, exid2loc, exid2strand, PSImax, PSImin, LPSI, RPSI):
    import pandas as pd
    pd.options.mode.chained_assignment = None
    import re
    import time

    chrstr = [x for x in re.split(r'[/\\]', sdir) if x != ''][-1] + '.'

    steps = [0, 0, 0, 0, 0, 0]
    for x in range(argSteps):
        steps[x] = 1

    # ============================================
    clst = open(sdir + '/rsr-classify.log', 'w')
    clst.write('rsr-classify.py run started on ' + time.strftime('%m/%d/%y, %H:%M:%S') + '\n')

    # Step0: Load y_exLR, Setup df_main, report constitutive
    if steps[0] == 1:
        # Inclusion/Exclusion Junction Coordinates (No PSI info in these files)
        # y_exLR = 1st/2nd/3rd highest EXCLUSION junctions (start/end/percentage) for left and right
        with open(sdir + '/' + argInfile, 'r') as f:
            df_main = pd.read_csv(f, sep=',', index_col=False,
                                  names=['ExonID',
                                         '1_ELS', '1_ELE', '1_ELP',
                                         '2_ELS', '2_ELE', '2_ELP',
                                         '3_ELS', '3_ELE', '3_ELP',
                                         '1_ERS', '1_ERE', '1_ERP',
                                         '2_ERS', '2_ERE', '2_ERP',
                                         '3_ERS', '3_ERE', '3_ERP'])
        df_main['ExonChr'] = df_main['ExonID'].map(exid2loc)
        df_main[['ExonChr', 'ExonStart', 'ExonEnd']] = df_main['ExonChr'].str.split(r':|-', expand=True)
        df_main['ExonStrand'] = df_main['ExonID'].map(exid2strand)
        df_main['ExonClass'] = 'Unclassified'
        df_main['MergedCount'] = 0
        df_main['MetadataLASS'] = 'None'
        df_main['MetadataRASS'] = 'None'

        # Add columns for detecting alt splice sites that should be linked together
        # LASS = Left Alt Splice Sites, RASS = Right Alt Splice Sites
        df_main['LASS'] = df_main[['ExonChr', '1_ELS', 'ExonEnd', '1_ERE', 'ExonStrand']].astype(str).apply(
            lambda x: ':'.join(x), axis=1)
        df_main['RASS'] = df_main[['ExonChr', '1_ELS', 'ExonStart', '1_ERE', 'ExonStrand']].astype(str).apply(
            lambda x: ':'.join(x), axis=1)

        df_main['PSImax'] = df_main['ExonID'].map(PSImax)
        df_main['PSImin'] = df_main['ExonID'].map(PSImin)
        df_main[['LPSImax1', 'LPSImax2', 'LPSImax3', 'LPSImax4', 'LPSImax5',
                 'RPSImax1', 'RPSImax2', 'RPSImax3', 'RPSImax4', 'RPSImax5']] = df_main['PSImax'].str.split(',',
                                                                                                            expand=True)
        df_main[['LPSImin1', 'LPSImin2', 'LPSImin3', 'LPSImin4', 'LPSImin5',
                 'RPSImin1', 'RPSImin2', 'RPSImin3', 'RPSImin4', 'RPSImin5']] = df_main['PSImin'].str.split(',',
                                                                                                            expand=True)
        df_constitutive = df_main.copy()
        df_main = df_main.drop(['PSImax', 'PSImin',
                                'LPSImax2', 'LPSImax3', 'LPSImax4', 'LPSImax5',
                                'RPSImax2', 'RPSImax3', 'RPSImax4', 'RPSImax5',
                                'LPSImin2', 'LPSImin3', 'LPSImin4', 'LPSImin5',
                                'RPSImin2', 'RPSImin3', 'RPSImin4', 'RPSImin5'], axis=1)
        df_main = df_main.reset_index(drop=True)
        df_main.to_csv(sdir + '/df1_exLR.csv', sep=',')

        keepList = []
        minPSI = 90
        for ix in df_constitutive.index:
            if str(df_constitutive['1_ELP'][ix]) == '999':
                if str(df_constitutive['1_ERP'][ix]) == '999':
                    keepList.append(ix)
                elif float(df_constitutive['RPSImin1'][ix]) >= minPSI:
                    if float(df_constitutive['RPSImin5'][ix]) >= minPSI:
                        keepList.append(ix)
            elif str(df_constitutive['1_ERP'][ix]) == '999':
                if str(df_constitutive['1_ELP'][ix]) == '999':
                    keepList.append(ix)
                elif float(df_constitutive['LPSImin1'][ix]) >= minPSI:
                    if float(df_constitutive['LPSImin5'][ix]) >= minPSI:
                        keepList.append(ix)
            elif float(df_constitutive['LPSImin1'][ix]) >= minPSI:
                if float(df_constitutive['LPSImin5'][ix]) >= minPSI:
                    if float(df_constitutive['RPSImin1'][ix]) >= minPSI:
                        if float(df_constitutive['RPSImin5'][ix]) >= minPSI:
                            keepList.append(ix)

        df_constitutive = df_constitutive.loc[keepList]
        df_constitutive = df_constitutive.reset_index(drop=True)
        df_constitutive = df_constitutive[['ExonID', 'ExonChr', 'ExonStart', 'ExonEnd',
                                           '1_ELS', '1_ERE', 'ExonStrand']]
        df_constitutive['LeftPSIstr'] = 'None'
        df_constitutive['RightPSIstr'] = 'None'
        for ix in df_constitutive.index:
            df_constitutive['LeftPSIstr'][ix] = LPSI[df_constitutive['ExonID'][ix]]
            df_constitutive['RightPSIstr'][ix] = RPSI[df_constitutive['ExonID'][ix]]
        df_constitutive.to_csv(sdir + '/df_constitutive.csv')

    # Step1: Create LASS/RASS Records
    if steps[1] == 1:
        with open(sdir + '/df1_exLR.csv', 'r') as f:
            df_main = pd.read_csv(f, sep=',', index_col=0)

        # Make copy of df_main to append LASS/RASS records
        df_LASSRASS = df_main.copy()
        PSIthreshold = 25  # Minimum LPSI (for LASS) or RPSI (for RASS) for an exon to be merged
        # ------------------------------------------------
        # ***************** LASS SECTION *****************
        # ------------------------------------------------
        # Merge associated alt splice site exons into single cassette records
        df_main = df_main.sort_values(by=['LASS', 'LPSImax1'], ascending=[True, False])
        df_main = df_main.reset_index(drop=True)
        ExonList = [df_main['ExonID'][0]]  # Track ExonIDs to merge
        RowSlice = df_main.loc[[0, 1]]  # Need two rows to stop conversion to series
        RowSlice.index = ['a', 'b']
        RecordCounter = 1  # Tracks new cassette exon records derived from merged records
        ex2ix = {}  # Dictionary to convert ExonID to df_main index
        for ix in range(len(df_main.index)):
            ex2ix[df_main['ExonID'][ix]] = ix
        LRex2ix = {}  # Dictionary to convert ExonID to df_LASSRASS index
        for ix in range(len(df_LASSRASS.index)):
            LRex2ix[df_LASSRASS['ExonID'][ix]] = ix

        for ix in range(1, len(df_main.index)):
            if df_main['LASS'][ix] == RowSlice['LASS']['a']:
                if df_main['LPSImax1'][ix] >= PSIthreshold:
                    ExonList.append(df_main['ExonID'][ix])
                else:
                    eid = df_main['ExonID'][ix]
                    df_LASSRASS['ExonClass'][LRex2ix[eid]] = 'Filtered:LowPSI:LASS'
            else:
                if len(ExonList) > 1:
                    clst.write('LASS recorded at index ' + str(ix) + time.strftime(', %m/%d/%y, %H:%M:%S') + '\n')

                    # Append record to df_main
                    ExonStartSites = []
                    ELElist = []
                    for eid in ExonList:
                        exst = int(df_main['ExonStart'][ex2ix[eid]])
                        ExonStartSites.append(exst)
                        ELElist.append(exst - 1)
                        df_LASSRASS['ExonClass'][LRex2ix[eid]] = 'Merged:' + 'LASS_' + chrstr + str(
                            RecordCounter).zfill(7)
                        df_LASSRASS['MergedCount'][LRex2ix[eid]] = df_LASSRASS['MergedCount'][LRex2ix[eid]] + 1
                    RowSlice['ExonStart']['a'] = '|'.join(str(x) for x in ExonStartSites)
                    RowSlice['ExonID']['a'] = 'LASS_' + chrstr + str(RecordCounter).zfill(7)
                    RowSlice['MetadataLASS']['a'] = ','.join(ExonList)
                    RowSlice['ExonClass']['a'] = 'LASS'

                    looplist = ['1_ELE', '2_ELE', '3_ELE']
                    for x in looplist:
                        if int(RowSlice[x]['a']) not in ELElist:
                            if x != '1_ELE':
                                RowSlice['1_ELE']['a'] = RowSlice[x]['a']
                                RowSlice['1_ELP']['a'] = 222
                                RowSlice[x]['a'] = 'copied to 1_ELE'
                            break
                        if x == looplist[-1]:
                            RowSlice['1_ELE']['a'] = -1  # unable to find true 1_ELE
                            RowSlice['1_ELP']['a'] = -111

                    df_LASSRASS = df_LASSRASS.append(RowSlice.loc['a'])
                    RecordCounter = RecordCounter + 1
                    # Setup next record
                    ExonList = [df_main['ExonID'][ix]]
                    RowSlice.loc['a'] = df_main.loc[ix]
                    RowSlice.index = ['a', 'b']
                    if df_main['LPSImax1'][ix] < PSIthreshold:
                        RowSlice['LASS']['a'] = 'LowPSI'
                else:
                    # Setup next record, do not append
                    ExonList = [df_main['ExonID'][ix]]
                    RowSlice.loc['a'] = df_main.loc[ix]
                    RowSlice.index = ['a', 'b']
                    if df_main['LPSImax1'][ix] < PSIthreshold:
                        RowSlice['LASS']['a'] = 'LowPSI'

        # ------------------------------------------------
        # ***************** RASS SECTION *****************
        # ------------------------------------------------
        # Merge associated alt splice site exons into single cassette records
        df_main = df_main.sort_values(by=['RASS', 'RPSImax1'], ascending=[True, False])
        df_main = df_main.reset_index(drop=True)
        ExonList = [df_main['ExonID'][0]]  # Track ExonIDs to merge
        RowSlice = df_main.loc[[0, 1]]  # Need two rows to stop conversion to series
        RowSlice.index = ['a', 'b']
        RecordCounter = 1  # Tracks new cassette exon records derived from merged records
        ex2ix = {}  # Dictionary to convert ExonID to df_main index
        for ix in range(len(df_main.index)):
            ex2ix[df_main['ExonID'][ix]] = ix
        LRex2ix = {}  # Dictionary to convert ExonID to df_LASSRASS index
        for ix in range(len(df_LASSRASS.index)):
            LRex2ix[df_LASSRASS['ExonID'][ix]] = ix

        for ix in range(1, len(df_main.index)):
            if df_main['RASS'][ix] == RowSlice['RASS']['a']:
                if df_main['RPSImax1'][ix] >= PSIthreshold:
                    ExonList.append(df_main['ExonID'][ix])
                else:
                    eid = df_main['ExonID'][ix]
                    df_LASSRASS['ExonClass'][LRex2ix[eid]] = 'Filtered:LowPSI:RASS'
            else:
                if len(ExonList) > 1:
                    clst.write('RASS recorded at index ' + str(ix) + time.strftime(', %m/%d/%y, %H:%M:%S') + '\n')

                    # Append record to df_main
                    ExonEndSites = []
                    ERSlist = []
                    for eid in ExonList:
                        exen = int(df_main['ExonEnd'][ex2ix[eid]])
                        ExonEndSites.append(exen)
                        ERSlist.append(exen + 1)
                        df_LASSRASS['ExonClass'][LRex2ix[eid]] = 'Merged:' + 'RASS_' + chrstr + str(
                            RecordCounter).zfill(7)
                        df_LASSRASS['MergedCount'][LRex2ix[eid]] = df_LASSRASS['MergedCount'][LRex2ix[eid]] + 1
                    RowSlice['ExonEnd']['a'] = '|'.join(str(x) for x in ExonEndSites)
                    RowSlice['ExonID']['a'] = 'RASS_' + chrstr + str(RecordCounter).zfill(7)
                    RowSlice['MetadataRASS']['a'] = ','.join(ExonList)
                    RowSlice['ExonClass']['a'] = 'RASS'

                    looplist = ['1_ERS', '2_ERS', '3_ERS']
                    for x in looplist:
                        if int(RowSlice[x]['a']) not in ERSlist:
                            if x != '1_ERS':
                                RowSlice['1_ERS']['a'] = RowSlice[x]['a']
                                RowSlice['1_ERP']['a'] = 222
                                RowSlice[x]['a'] = 'copied to 1_ERS'
                            break
                        if x == looplist[-1]:
                            RowSlice['1_ERS']['a'] = -1  # unable to find true 1_ERS
                            RowSlice['1_ERP']['a'] = -111

                    df_LASSRASS = df_LASSRASS.append(RowSlice.loc['a'])
                    RecordCounter = RecordCounter + 1
                    # Setup next record
                    ExonList = [df_main['ExonID'][ix]]
                    RowSlice.loc['a'] = df_main.loc[ix]
                    RowSlice.index = ['a', 'b']
                    if df_main['RPSImax1'][ix] < PSIthreshold:
                        RowSlice['RASS']['a'] = 'LowPSI'
                else:
                    # Setup next record, do not append
                    ExonList = [df_main['ExonID'][ix]]
                    RowSlice.loc['a'] = df_main.loc[ix]
                    RowSlice.index = ['a', 'b']
                    if df_main['RPSImax1'][ix] < PSIthreshold:
                        RowSlice['RASS']['a'] = 'LowPSI'

        df_LASSRASS = df_LASSRASS.drop(['LASS', 'RASS'], axis=1)
        df_main = df_LASSRASS.reset_index(drop=True)

        df_main = df_main.reset_index(drop=True)
        df_main.to_csv(sdir + '/df2_LASSRASS.csv', sep=',')

    # Step2: Find Cassette Exons
    if steps[2] == 1:
        with open(sdir + '/df2_LASSRASS.csv', 'r') as f:
            df_main = pd.read_csv(f, sep=',', index_col=0)
        df_main['Cassette'] = 'No'
        df_main = df_main.reset_index(drop=True)

        for ix in range(len(df_main.index)):
            if ix % 1000 == 0:
                clst.write('CASSETTE: ' + str(ix) + time.strftime(', %m/%d/%y, %H:%M:%S') + '\n')
            if df_main['ExonClass'][ix] in ['Unclassified', 'LASS', 'RASS']:
                if df_main['1_ELS'][ix] == df_main['1_ERS'][ix] and df_main['1_ELE'][ix] == df_main['1_ERE'][ix]:
                    df_main['Cassette'][ix] = 'Yes'
                    df_main['ExonClass'][ix] = 'Cassette'
        df_main = df_main.reset_index(drop=True)
        df_main.to_csv(sdir + '/df3_CASSETTE.csv', sep=',')

    # Step3: Find Linked Exons
    if steps[3] == 1:
        with open(sdir + '/df3_CASSETTE.csv', 'r') as f:
            df_main = pd.read_csv(f, sep=',', index_col=0)
        df_main['LINKLEFT'] = df_main[['ExonChr', '1_ELS', '1_ELE', 'ExonStrand']].astype(str).apply(
            lambda x: ':'.join(x), axis=1)
        df_main['LINKRIGHT'] = df_main[['ExonChr', '1_ERS', '1_ERE', 'ExonStrand']].astype(str).apply(
            lambda x: ':'.join(x), axis=1)
        df_main['MetadataLINKED'] = 'None'

        df_main = df_main.sort_values(['ExonChr', '1_ERS', '1_ELS', '1_ERE'], ascending=[True, True, True, True])
        # Sorting needed to ensure order is correct so that no record is inadvertently skipped in comparison to 'LINKRIGHT'
        df_main = df_main.reset_index(drop=True)

        # -------------------------------------------------------------
        # Loop over each LINKLEFT row and search for match in LINKRIGHT
        # -------------------------------------------------------------
        RowSlice = df_main.loc[[0, 1]]
        RowSlice.index = ['a', 'b']
        RecordCounter = 1
        ex2ix = {}
        for ix in range(len(df_main.index)):
            ex2ix[df_main['ExonID'][ix]] = ix
        # Create dict to set index bounds for LINKRIGHT match search
        # Inclusive so if x = chrbounds['chr1'].split(','), range( x[0], x[1]+1 )
        chrbounds = {}
        chrid = df_main['ExonChr'][0]
        chrix = [0, 1]
        for ix in range(len(df_main.index)):
            if df_main['ExonChr'][ix] != chrid:
                chrix[1] = ix - 1
                chrbounds[chrid] = ','.join(str(x) for x in chrix)
                chrid = df_main['ExonChr'][ix]
                chrix[0] = ix
            if ix == len(df_main.index) - 1:
                chrix[1] = ix
                chrbounds[chrid] = ','.join(str(x) for x in chrix)

        PSIthreshold = 10
        for ix in range(len(df_main.index)):
            if df_main['ExonClass'][ix] in ['Unclassified', 'LASS', 'RASS']:
                if df_main['Cassette'][ix] == 'No':
                    # Need [df_main['1_ELE'][ix] != -1] and especially [df_main['1_ERS'][ix] != -1]
                    # because 1_ERS=-1 messes up the purpose of 1_ERS sorting and results in slow searching
                    if df_main['1_ELP'][ix] != 999 and df_main['1_ELE'][ix] != -1 and df_main['1_ERS'][ix] != -1:
                        if df_main['LPSImax1'][ix] > PSIthreshold:
                            leftquery = df_main['LINKLEFT'][ix]
                            srange = [int(x) for x in chrbounds[df_main['ExonChr'][ix]].split(',')]
                            for ix2 in range(ix, srange[1] + 1):
                                leftExonEnd = df_main['ExonEnd'][ix]
                                if type(leftExonEnd) is str:
                                    leftExonEnd = [int(x) for x in df_main['ExonEnd'][ix].split('|')]
                                    leftExonEnd = max(leftExonEnd)
                                if int(df_main['1_ERS'][ix2]) > leftExonEnd:
                                    break
                                elif df_main['ExonClass'][ix2] in ['Unclassified', 'LASS', 'RASS']:
                                    if df_main['Cassette'][ix2] == 'No':
                                        if df_main['1_ERP'][ix2] != 999 and df_main['1_ERS'][ix2] != -1:
                                            if df_main['RPSImax1'][ix2] > PSIthreshold:
                                                v1 = df_main['ExonStart'][ix2]
                                                if type(v1) is str:
                                                    v1 = [int(x) for x in v1.split('|')]
                                                    v1 = min(v1)
                                                v2 = df_main['ExonEnd'][ix]
                                                if type(v2) is str:
                                                    v2 = [int(x) for x in v2.split('|')]
                                                    v2 = max(v2)
                                                if v1 > v2:
                                                    if df_main['LINKRIGHT'][ix2] == leftquery:
                                                        RowSlice.loc['a'] = df_main.loc[ix]
                                                        RowSlice['MetadataLINKED']['a'] = df_main['ExonID'][
                                                                                              ix] + ',' + \
                                                                                          df_main['ExonID'][ix2]
                                                        df_main['MetadataLINKED'][
                                                            ix] = 'Merged:LINKED_' + chrstr + str(
                                                            RecordCounter).zfill(7)
                                                        df_main['MetadataLINKED'][
                                                            ix2] = 'Merged:LINKED_' + chrstr + str(
                                                            RecordCounter).zfill(7)
                                                        RowSlice['ExonEnd']['a'] = df_main['ExonEnd'][ix2]
                                                        RowSlice['1_ERS']['a'] = df_main['1_ELS'][ix]
                                                        RowSlice['1_ERE']['a'] = df_main['1_ELE'][ix]
                                                        RowSlice['ExonID']['a'] = 'LINKED_' + chrstr + str(
                                                            RecordCounter).zfill(7)
                                                        RowSlice['ExonClass']['a'] = 'LINKED'
                                                        df_main = df_main.append(RowSlice.loc['a'])
                                                        clst.write('LINKED recorded at ' + str(ix) + ', ' +
                                                                   str(RowSlice['ExonChr']['a']) + ':' +
                                                                   str(RowSlice['ExonStart']['a']) + '-' +
                                                                   str(RowSlice['ExonEnd']['a']) +
                                                                   time.strftime(', %m/%d/%y, %H:%M:%S') + '\n')
                                                        RecordCounter = RecordCounter + 1

        # Drop duplicate LINKED records
        dropList = []
        prev = ''
        df_main = df_main.reset_index(drop=True)
        df_main['LINKDUPLICATED'] = df_main[['ExonChr', 'ExonStart', 'ExonEnd', 'ExonStrand']].astype(str).apply(
            lambda x: ':'.join(x), axis=1)
        for ix in range(len(df_main.index)):
            if df_main['ExonClass'][ix] == 'LINKED':
                if df_main['LINKDUPLICATED'][ix] == prev:
                    dropList.append(ix)
                else:
                    prev = df_main['LINKDUPLICATED'][ix]
        df_main.drop(df_main.index[dropList], inplace=True)

        df_main = df_main.drop(['LINKLEFT', 'LINKRIGHT', 'LINKDUPLICATED'], axis=1)
        df_main = df_main.reset_index(drop=True)
        df_main.to_csv(sdir + '/df4_LINKED.csv', sep=',')

    # Step4: Find Mutually Exclusive Regions
    if steps[4] == 1:
        with open(sdir + '/df4_LINKED.csv', 'r') as f:
            df_main = pd.read_csv(f, sep=',', index_col=0)
        df_main['MetadataMUTEX'] = 'None'
        df_main['MUTEX'] = df_main[['ExonChr', '1_ELS', '1_ERE', 'ExonStrand']].astype(str).apply(
            lambda x: ':'.join(x), axis=1)
        df_main = df_main.sort_values(by=['MUTEX', 'ExonStart'], ascending=[True, True])
        df_main = df_main.reset_index(drop=True)

        RecordCounter = 1
        MUTEXid = 'empty'
        ixList = []
        PSIthreshold = 30

        for ix in range(len(df_main.index)):
            if ix % 1000 == 0:
                clst.write('MUTEX: ' + str(ix) + time.strftime(', %m/%d/%y, %H:%M:%S') + '\n')
            if df_main['MUTEX'][ix] == MUTEXid:
                cond1 = df_main['ExonClass'][ix] not in ['Cassette', 'LASS', 'RASS', 'LINKED']
                cond2 = df_main['MetadataLASS'][ix] == 'None'
                cond3 = df_main['MetadataRASS'][ix] == 'None'
                if cond1 and cond2 and cond3:
                    maxPSI = min([int(df_main['LPSImax1'][ix]), int(df_main['RPSImax1'][ix])])
                    if maxPSI > PSIthreshold:
                        ixList.append(ix)
            else:
                cond1 = df_main['ExonClass'][ix] not in ['Cassette', 'LASS', 'RASS', 'LINKED']
                cond2 = df_main['MetadataLASS'][ix] == 'None'
                cond3 = df_main['MetadataRASS'][ix] == 'None'
                if cond1 and cond2 and cond3:
                    if len(ixList) > 1:
                        for a in range(len(ixList)-1):
                            for b in list(range(len(ixList)))[a:]:
                                mutex1 = int(df_main['1_ELE'][ixList[a]]) == int(df_main['ExonStart'][ixList[b]]) - 1
                                mutex2 = int(df_main['1_ELE'][ixList[b]]) == int(df_main['ExonStart'][ixList[a]]) - 1
                                mutex3 = int(df_main['1_ERS'][ixList[a]]) == int(df_main['ExonEnd'][ixList[b]]) + 1
                                mutex4 = int(df_main['1_ERS'][ixList[b]]) == int(df_main['ExonEnd'][ixList[a]]) + 1
                                mutex5 = int(df_main['ExonEnd'][ixList[a]]) < int(df_main['ExonStart'][ixList[b]])
                                if mutex1 and mutex2 and mutex3 and mutex4 and mutex5:
                                    df_main['MetadataMUTEX'][ixList[a]] = 'MUTEX_' + chrstr + str(RecordCounter).zfill(7)
                                    df_main['MetadataMUTEX'][ixList[b]] = 'MUTEX_' + chrstr + str(RecordCounter).zfill(7)
                                    RecordCounter = RecordCounter + 1
                        MUTEXid = df_main['MUTEX'][ix]
                        ixList = [ix]
                    else:
                        # Setup next record
                        MUTEXid = df_main['MUTEX'][ix]
                        ixList = [ix]
            if ix == len(df_main.index) - 1:
                if len(ixList) > 1:
                    for a in range(len(ixList) - 1):
                        for b in list(range(len(ixList)))[a:]:
                            mutex1 = int(df_main['1_ELE'][ixList[a]]) == int(df_main['ExonStart'][ixList[b]]) - 1
                            mutex2 = int(df_main['1_ELE'][ixList[b]]) == int(df_main['ExonStart'][ixList[a]]) - 1
                            mutex3 = int(df_main['1_ERS'][ixList[a]]) == int(df_main['ExonEnd'][ixList[b]]) + 1
                            mutex4 = int(df_main['1_ERS'][ixList[b]]) == int(df_main['ExonEnd'][ixList[a]]) + 1
                            mutex5 = int(df_main['ExonEnd'][ixList[a]]) < int(df_main['ExonStart'][ixList[b]])
                            if mutex1 and mutex2 and mutex3 and mutex4 and mutex5:
                                df_main['MetadataMUTEX'][ixList[a]] = 'MUTEX_' + chrstr + str(RecordCounter).zfill(7)
                                df_main['MetadataMUTEX'][ixList[b]] = 'MUTEX_' + chrstr + str(RecordCounter).zfill(7)
                                RecordCounter = RecordCounter + 1
        df_main = df_main.drop('MUTEX', axis=1)
        df_main = df_main.reset_index(drop=True)
        df_main.to_csv(sdir + '/df5_MUTEX.csv', sep=',')

    # Step5: Refine List
    if steps[5] == 1:
        with open(sdir + '/df5_MUTEX.csv', 'r') as f:
            df_main = pd.read_csv(f, sep=',', index_col=0)
        df_main = df_main.filter(['ExonID', 'ExonChr', 'ExonStart', 'ExonEnd', 'ExonStrand',
                                  '1_ELS', '1_ERE', 'LPSImax1', 'RPSImax1', 'LPSImin1', 'RPSImin1',
                                  'ExonClass', 'Cassette',
                                  'MetadataLASS', 'MetadataRASS', 'MetadataLINKED', 'MetadataMUTEX'], axis=1)
        df_main.to_csv(sdir + '/df6_REFINE.csv', sep=',')

    # ============================================
    clst.close()


def writeNewInput(wd, ic, dl):
    prevFilename = wd + '/y_input_v' + str(ic) + '.csv'
    newFilename = wd + '/y_input_v' + str(ic + 1) + '.csv'
    nf = open(newFilename, 'w')
    with open(prevFilename, 'r') as f:
        for line in f:
            sp = line.split(',')
            if sp[0] not in dl:
                nf.write(line)
    nf.close()


def writeTimestamp(wd, checkcount, ct):
    import time
    with open(wd + '/rsr-cleanup.log', 'a') as f:
        dt = time.time() - ct
        f.write('Checkpoint ' + str(checkcount) + ', ' + time.strftime("%H:%M:%S") +
                ', runtime: ' + str(round(dt, 4)) + ' sec\n')
        
