import argparse
from argparse import RawTextHelpFormatter
import configparser
import gzip
import pandas as pd
import shutil
import subprocess
import sys
import os
import urllib.request

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
=================================================================
 Calculate NAUC (gene coverage) for a given snaptron compilation
=================================================================
""")
parser.add_argument('--a', action='store',
                    required=True,
                    help='Dataset config with all samples')
parser.add_argument('--c', action='store',
                    required=True,
                    help='Dataset config with linked samples')
parser.add_argument('--o', action='store',
                    default='ascot_nauc_output.tsv',
                    help='Output filename, default \'ascot_nauc_output.tsv\'')
args = parser.parse_args()

# =======================================================================
# Parameters
all_samples = args.a.upper()
config_linked = args.c.upper()
nauc_output = args.o
comp2species = {
    'supermouse':'mouse',
    'ct_m_s':'mouse',
    'encode1159':'human',
    'gtex':'human',
    'srav2':'human'
    }
cwd = os.getcwd()
cfg = configparser.ConfigParser()
with open(cwd + '/cfg/master_config.ini', 'r') as f:
    cfg.read_file(f)
if all_samples not in cfg.sections():
    sys.exit('--a Config section not found')
if config_linked not in cfg.sections():
    sys.exit('--c Config section not found')
dsrc = cfg[config_linked]['snaptron_datasrc'].split(',')
species = comp2species[dsrc[0]]

# =======================================================================
# Check for GTFs in /gtf
# Human: gencode.v25.basic.annotation.gtf (hg38)
# Mouse: gencode.vM15.basic.annotation.gtf (mm10)
if not os.path.exists('./gtf'):
    os.makedirs('./gtf')
if species == 'human':
    gtf = './gtf/gencode.v25.basic.annotation.gtf'
    if not os.path.exists(gtf):
        print('Human GTF required, GENCODE v25 basic gene annotation not found')
        print('Downloading \'gencode.v25.basic.annotation.gtf\' to /gtf ...')
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_25/gencode.v25.basic.annotation.gtf.gz'
        urllib.request.urlretrieve(url, gtf + '.gz')  
        with gzip.open(gtf + '.gz', 'rb') as fi:
            with open('./gtf/gencode.v25.basic.annotation.gtf', 'wb') as fo:
                shutil.copyfileobj(fi, fo)
        os.remove(gtf + '.gz')
    if not os.path.exists(gtf):
        sys.exit('Unable to download \'gencode.v25.basic.annotation.gtf\' to /gtf')
elif species == 'mouse':
    gtf = './gtf/gencode.vM15.basic.annotation.gtf'
    if not os.path.exists(gtf):
        print('Mouse GTF required, GENCODE vM15 basic gene annotation not found')
        print('Downloading \'gencode.vM15.basic.annotation.gtf\' to /gtf ...')
        url = 'ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M15/gencode.vM15.basic.annotation.gtf.gz'
        urllib.request.urlretrieve(url, gtf + '.gz')
        with gzip.open(gtf + '.gz', 'rb') as fi:
            with open('./gtf/gencode.vM15.basic.annotation.gtf', 'wb') as fo:
                shutil.copyfileobj(fi, fo)
        os.remove(gtf + '.gz')
    if not os.path.exists(gtf):
        sys.exit('Unable to download \'gencode.vM15.basic.annotation.gtf\' to /gtf')
else:
    sys.exit('Species must be \'human\' or \'mouse\'')

# =======================================================================
# Check for snaptron gene coverage file in /nauc
if not os.path.exists('./nauc'):
    os.makedirs('./nauc')
gcli = []
for snaptron_compilation in dsrc:
    nauc = './nauc/' + snaptron_compilation + '_gene_coverage_normalized.tsv'
    if not os.path.exists(nauc):
        print('Snaptron gene coverage file not found in /nauc (' + snaptron_compilation + '_gene_coverage_normalized.tsv)')
        print('Downloading \'' + snaptron_compilation + '_gene_coverage_normalized.tsv\' to /nauc ...')
        url = 'http://snaptron.cs.jhu.edu/data/' + snaptron_compilation + '/gene_coverage_normalized.tsv.bgz'
        urllib.request.urlretrieve(url, nauc + '.bgz')
        with gzip.open(nauc + '.bgz', 'rb') as fi:
            with open(nauc, 'wb') as fo:
                shutil.copyfileobj(fi, fo)
        os.remove(nauc + '.bgz')
    if not os.path.exists(nauc):
        sys.exit('Unable to download \'' + snaptron_compilation + '_gene_coverage_normalized.tsv\' to /gtf')
    if len(dsrc) > 1:    
        gcli.append(pd.read_csv(nauc, sep='\t', index_col=0, header=0))

print('Processing gene coverage file(s)')
if len(dsrc) > 1:
    nauc = './nauc/' + cfg[config_linked]['snaptron_datasrc'] + '_gene_coverage_normalized.tsv'
    if not os.path.exists(nauc):
        df_gcli = pd.concat(gcli, axis=1)
        df_gcli.index.name = 'gene_id'
        df_gcli.to_csv(nauc, sep='\t')

# =======================================================================
subprocess.run('python3 ./bin/rsrNAUC.py' +
               ' --cfgall ' + all_samples +
               ' --cfglinked ' + config_linked +
               ' --gc ' + nauc +
               ' --gco ' + nauc_output +
               ' --gtf ' + gtf, shell=True)
