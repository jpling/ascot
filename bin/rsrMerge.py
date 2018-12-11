import argparse
from argparse import RawTextHelpFormatter
import configparser
import os
import numpy as np
import re
import pandas as pd
pd.options.mode.chained_assignment = None
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
============================================================
 Merge samples based on config, takes input from rsrMain.py 
============================================================
""")
parser.add_argument('--cfglinked', action='store',
                    required=True,
                    help='Linked config preset')
parser.add_argument('--cfgall', action='store',
                    required=True,
                    help='Original config preset with all samples')
parser.add_argument('--i', action='store',
                    required=True,
                    help='Input exon metadata file')
parser.add_argument('--o', action='store',
                    default='merged_exon_psi.tsv',
                    help='Output filename')
parser.add_argument('--min', action='store',
                    type=int,
                    default=30,
                    help='Minimum junction count')
parser.add_argument('--f', action='store',
                    type=float,
                    default=0.5,
                    help='PSI fraction cutoff')
args = parser.parse_args()

# ==========================================================
linkedPreset = args.cfglinked
allPreset = args.cfgall
exondata = args.i
output = args.o
minJC = args.min

# ==========================================================
cwd = os.getcwd()
cfg = configparser.ConfigParser()
with open(cwd + '/cfg/master_config.ini', 'r') as f:
    cfg.read_file(f)
sLabels = cfg[linkedPreset]['sample_labels'].split('|')
linkedSamples = cfg[linkedPreset]['snaptron_samples'].split(',')
sColumns = []
for a in linkedSamples:
    sColumns.append([int(x) for x in a.split('-')])

# ==========================================================
# Merge sample PSI
# ==========================================================

df_main = pd.read_csv(exondata, sep='\t')
# Available metadata:
#   'ExonID', 'CASSETTE', 'LASS', 'RASS', 'LINKED', 'MUTEX',
#   'GeneID', 'GeneSymbol', 'GeneType',
#   'ExonLocation', 'ExonBoundary', 'ExonStrand', 'ExonLength',
#   'LeftJCstr', 'RightJCstr', 'LeftPSIstr', 'RightPSIstr'
columns = list(map(int, cfg[allPreset]['snaptron_samples'].split(',')))
df_lpsi = df_main[['ExonID', 'LeftPSIstr']]
df_rpsi = df_main[['ExonID', 'RightPSIstr']]
df_ljc = df_main[['ExonID', 'LeftJCstr']]
df_rjc = df_main[['ExonID', 'RightJCstr']]

df_lpsi[columns] = df_lpsi['LeftPSIstr'].str.split(',', expand=True).astype(float)
print('Split LeftPSIstr')
df_rpsi[columns] = df_rpsi['RightPSIstr'].str.split(',', expand=True).astype(float)
print('Split RightPSIstr')
df_ljc[columns] = df_ljc['LeftJCstr'].str.split(',', expand=True).astype(float)
print('Split LeftJCstr')
df_rjc[columns] = df_rjc['RightJCstr'].str.split(',', expand=True).astype(float)
print('Split RightJCstr')
df_lpsi.drop(['ExonID', 'LeftPSIstr'], axis=1, inplace=True)
df_rpsi.drop(['ExonID', 'RightPSIstr'], axis=1, inplace=True)
df_ljc.drop(['ExonID', 'LeftJCstr'], axis=1, inplace=True)
df_rjc.drop(['ExonID', 'RightJCstr'], axis=1, inplace=True)

PJL = df_lpsi * df_ljc.values
PJR = df_rpsi * df_rjc.values

for i in range(len(sLabels)):
    print('SamplePSI: ' + sLabels[i])
    print(sColumns[i])
    PJL[sLabels[i]] = PJL[sColumns[i]].sum(axis=1)
    PJR[sLabels[i]] = PJR[sColumns[i]].sum(axis=1)
    df_ljc[sLabels[i]] = df_ljc[sColumns[i]].sum(axis=1)
    df_rjc[sLabels[i]] = df_rjc[sColumns[i]].sum(axis=1)
PJL.drop(columns, axis=1, inplace=True)
PJR.drop(columns, axis=1, inplace=True)
df_ljc.drop(columns, axis=1, inplace=True)
df_rjc.drop(columns, axis=1, inplace=True)

ljc_filter = df_ljc.copy()
ljc_filter[ljc_filter < minJC] = -1
ljc_filter[ljc_filter >= minJC] = 1
rjc_filter = df_rjc.copy()
rjc_filter[rjc_filter < minJC] = -1
rjc_filter[rjc_filter >= minJC] = 1
df_ljc = df_ljc * ljc_filter.values
df_rjc = df_rjc * rjc_filter.values

newLPSI = PJL / df_ljc.values
newRPSI = PJR / df_rjc.values
newLPSI[newLPSI < 0] = np.inf
newRPSI[newRPSI < 0] = np.inf
newLPSI = newLPSI.fillna(value=np.inf)
newRPSI = newRPSI.fillna(value=np.inf)
maxLRPSI = pd.concat([newLPSI, newRPSI]).max(level=0)
maxLRPSI[maxLRPSI < 20] = np.NaN
newAPSI = (newLPSI + newRPSI.values) / 2
newAPSI[newAPSI == np.inf] = -1

deltaPSI = newLPSI - newRPSI.values
deltaPSI = deltaPSI.fillna(value=np.inf)
deltaPSI = deltaPSI.abs()
deltaPSI[deltaPSI == np.inf] = np.NaN
deltaPSI = deltaPSI / maxLRPSI.values
dprQuantiles = deltaPSI.quantile([0.01, 0.25, 0.5, 0.75, 0.99], axis=1)
dprQuantiles = dprQuantiles.transpose()
dprQuantiles.columns = ['deltaPSIq10', 'deltaPSIq25', 'deltaPSIq50', 'deltaPSIq75', 'deltaPSIq90']

df_main = pd.concat([df_main, dprQuantiles, newAPSI], axis=1)
c = ['ExonID', 'CASSETTE', 'LASS', 'RASS', 'LINKED', 'MUTEX', 'ExonStrand', 'ExonLength',
     'GeneType', 'GeneID', 'GeneSymbol', 'ExonLocation', 'ExonBoundary']
c.extend(['deltaPSIq10', 'deltaPSIq25', 'deltaPSIq50', 'deltaPSIq75', 'deltaPSIq90'])
c.extend(sLabels)
df_main = df_main[c]

# Cleanup
dpr_cutoff = args.f
emptyGeneSymbol = df_main[df_main['GeneSymbol'].isnull()].index.tolist()
df_main.drop(emptyGeneSymbol, inplace=True)
deltaPSIfilter = df_main[df_main['deltaPSIq50'] >= dpr_cutoff].index.tolist()
df_main.drop(deltaPSIfilter, inplace=True)
df_main.insert(loc=2, column='Alternative Splice Site Group', value='No')
df_main.loc[(df_main['LASS'] == 'Yes') & (df_main['ExonStrand'] == '+'), 'Alternative Splice Site Group'] = '3\''
df_main.loc[(df_main['LASS'] == 'Yes') & (df_main['ExonStrand'] == '-'), 'Alternative Splice Site Group'] = '5\''
df_main.loc[(df_main['RASS'] == 'Yes') & (df_main['ExonStrand'] == '+'), 'Alternative Splice Site Group'] = '5\''
df_main.loc[(df_main['RASS'] == 'Yes') & (df_main['ExonStrand'] == '-'), 'Alternative Splice Site Group'] = '3\''
df_main.drop(['deltaPSIq10', 'deltaPSIq25', 'deltaPSIq50', 'deltaPSIq75', 'deltaPSIq90', 'LASS', 'RASS'], axis=1, inplace=True)
df_main.rename(index=str, columns={
    'ExonID':'exon_id',
    'CASSETTE':'cassette_exon',
    'Alternative Splice Site Group':'alternative_splice_site_group',
    'LINKED':'linked_exons',
    'MUTEX':'mutually_exclusive_exons',
    'ExonStrand':'exon_strand',
    'ExonLength':'exon_length',
    'GeneType':'gene_type',
    'GeneID':'gene_id',
    'GeneSymbol':'gene_symbol',
    'ExonLocation':'exon_location',
    'ExonBoundary':'exon_boundary'
    }, inplace=True)
df_main.to_csv(output, sep='\t', index=False)
