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
 Calculate NAUC(normalized area under curve) gene coverage
============================================================
""")
parser.add_argument('--cfgall', action='store',
                    required=True,
                    help='Original config preset with all samples')
parser.add_argument('--cfglinked', action='store',
                    required=True,
                    help='Linked config preset')
parser.add_argument('--gc', action='store',
                    required=True,
                    help='Snaptron gene coverage normalized file')
parser.add_argument('--gco', action='store',
                    default='merged_gene_coverage.tsv',
                    help='Gene coverage output')
parser.add_argument('--gtf', action='store',
                    default='gencode.vM15.basic.annotation.gtf',
                    help='Reference GTF')
args = parser.parse_args()

# ==========================================================
linkedPreset = args.cfglinked
allPreset = args.cfgall
genecov = args.gc
nauc_output = args.gco
reference_gtf = args.gtf

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
# Merge gene coverage
# ==========================================================

# Parse gencode GTF, make name2id dictionary
geneid2genesymbol = {}
with open(reference_gtf, 'r') as gtf:
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
                        gene_id = attributes[ix + 1].split('.')[0]
                    if attributes[ix] == 'gene_type':
                        gene_type = attributes[ix + 1]
                    if attributes[ix] == 'gene_name':
                        gene_name = attributes[ix + 1]
                geneid2genesymbol[gene_id] = gene_name


df_genecov = pd.read_csv(genecov, sep='\t')
columns = df_genecov.columns.values.tolist()[1:]

# Convert sColumns and columns from int to str
for x in range(len(sColumns)):
    for y in range(len(sColumns[x])):
        sColumns[x][y] = str(sColumns[x][y])
for x in range(len(columns)):
    columns[x] = str(columns[x])

scalingfactor = 1000
for i in range(len(sLabels)):
    print('GeneCov: ' + sLabels[i])
    print(sColumns[i])
    df_genecov[sLabels[i]] = df_genecov[sColumns[i]].sum(axis=1)/(scalingfactor*len(sColumns[i]))
df_genecov.drop(columns, axis=1, inplace=True)

df_genecov[['gene_id', 'gid_version']] = df_genecov['gene_id'].str.split('.', expand=True)
df_genecov.drop(['gid_version'], axis=1, inplace=True)
df_genecov['gene_symbol'] = df_genecov['gene_id'].map(geneid2genesymbol)
col = list(df_genecov.columns.values)
col.remove('gene_symbol')
col.insert(1, 'gene_symbol')
df_genecov = df_genecov[col]
df_genecov.to_csv(nauc_output, sep='\t', index=False)
