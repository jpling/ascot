import argparse
from argparse import RawTextHelpFormatter
import configparser
import gzip
import shutil
import subprocess
import sys
import os
import urllib.request

parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
========================================================
 Run ASCOT scripts to identify exons and calculate PSIs
========================================================
""")
parser.add_argument('--i', action='store',
                    required=True,
                    help='Exon list input file')
parser.add_argument('--a', action='store',
                    required=True,
                    help='Dataset config with all samples')
parser.add_argument('--c', action='store',
                    required=True,
                    help='Dataset config with linked samples')
parser.add_argument('--o', action='store',
                    default='ascot_psi_output.tsv',
                    help='Output filename, default \'ascot_psi_output.tsv\'')
parser.add_argument('--p', action='store',
                    type=int,
                    default=2,
                    help='Processes to use')
parser.add_argument('--min', action='store',
                    type=int,
                    default=15,
                    help='Minimum junctions for PSI calculations')
parser.add_argument('--f', action='store',
                    type=float,
                    default=0.8,
                    help='PSI fraction cutoff')
parser.add_argument('--deletefirstoutput', action='store',
                    type=int,
                    default=1,
                    help='Delete initial output')
parser.add_argument('--deleteoutputdir', action='store',
                    type=int,
                    default=1,
                    help='Delete output directory')                        
args = parser.parse_args()

# =======================================================================
# Parameters
exon_list = args.i
all_samples = args.a.upper()
config_linked = args.c.upper()
psi_output = args.o
processes = args.p
if not os.path.isfile('./software/snaptron/qs'):
    sys.exit('Snaptron not installed:\n  cd ./software/snaptron\n  make')
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
# Check for GENCODE gtfs (v25/vM15)
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
subprocess.run('python3 ./bin/rsrMain.py' +
               ' --exons ' + exon_list +
               ' --p ' + str(processes) +
               ' --deldir ' + str(args.deleteoutputdir) +
               ' --gtf ' + gtf +
               ' --cfg ' + all_samples +
               ' --out ' + './unlinked_output.tsv', shell=True)
subprocess.run('python3 ./bin/rsrMerge.py' +
               ' --cfgall ' + all_samples +
               ' --cfglinked ' + config_linked +
               ' --i ' + './unlinked_output.tsv' +
               ' --o ' + psi_output +
               ' --min ' + str(args.min) + 
               ' --f ' + str(args.f), shell=True)

if args.deletefirstoutput == 1:
    os.remove('./unlinked_output.tsv')
os.rmdir('../snaptron_tmp')
