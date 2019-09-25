from __future__ import absolute_import, division, print_function
import sys
import os
import re
import configparser
import argparse
from argparse import RawTextHelpFormatter
import codecs
import numpy as np
import pandas as pd
import itertools
import csv
import tempfile
from builtins import bytes
import math
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

# %%%%%%%%%%%%%% Adjustable parameters %%%%%%%%%%%%%%
chromosome_list = [
    'chr1' , 'chr2' , 'chr3' , 'chr4' , 'chr5' , 'chr6' , 'chr7' , 'chr8' , 'chr9' , 'chr10',
    'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
    'chr21', 'chr22', 'chrX' , 'chrY' , 'chrM'
    ]  # only coordinates with these chromosomes will be processed

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Global variables
cwd = os.getcwd()
os.environ["CUR_DIR"] = cwd

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Inclusion/Exclusion Junction Reporting
if not os.path.exists('y_inLR.csv'):
    b = open('y_inLR.csv','w')
    b.close()
if not os.path.exists('y_exLR.csv'):
    b = open('y_exLR.csv','w')
    b.close()

# Parse arguments
parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter, description="""
  ====================================
   RSR - Find PSI for cassette exons.
  ====================================

  REQUIRES --input parameter; others are optional
  
  Input file is a text file with a list of queries, one per line.  Example syntax:
  
ExIDv1-0001	chr1:118512175-118512222	+
ExIDv1-0002	chr1:118512705-118512728	+
ExIDv1-0003	chr1:118523321-118523347	+
ExIDv1-0004	chr1:118531978-118532085	+

  (Reminder: the fields are seprated by tabs, not spaces)
  
  Snaptron parameters are stored in the config file:
    -Samples (snaptron_samples)
    -Datasource (snaptron_datasrc)
    -Query filter (snaptron_filter)
""")
parser.add_argument('--input', action='store',
                    metavar = 'path',
                    required=True,
                    help='TSV file with queries, one query per line: exid(tab)chr:begin-end(tab)strand')
parser.add_argument('--query-temp', action='store',
                    metavar='path',
                    required=False,
                    help='Path to file where first and second pass queries should be stored; by '
                         'default queries are stored in temporary files deleted upon exit.')
parser.add_argument('--cfg', action='store',
                    default='master',
                    dest='config_id',
                    metavar = '',
                    help='Config file identifier\n  e.g. -cfg 123456 for /cfg/123456_config.ini')
parser.add_argument('--preset', action='store',
                    default='MAIN',
                    dest='config_section',
                    metavar = '',
                    help='Optional: choose a preset section in config instead of \'MAIN\'\n  e.g. -preset SUPERMOUSE')
parser.add_argument('--inclusion-fraction', action='store',
                    type=float, default=0.7,
                    help='A junction is considered an inclusion junction if it accounts for at least '
                         'this fraction of the splicing events incident on the exon.')
parser.add_argument('--queries-per-batch', action='store',
                    type=int, default=50,
                    help='Number of queries to send to server in single "bulk" batch.  50 is maxmimum supported.')
parser.add_argument('--coverage-sum', action='store',
                    type=int, default=10,
                    help='Ignore junctions with less than this many total reads supporting it across selected samples.')
parser.add_argument('--max-length', action='store',
                    type=int, default=200000,
                    help='Ignore junctions spanning an intron longer than this.')
parser.add_argument('--keep-intermediates', action='store_true',
                    default=False,
                    help='Keep intermediate files containing queries and responses')
parser.add_argument('--verbose', action='store_true',
                    default=False,
                    help='Be talkative')
parser.add_argument('--debug', action='store_true',
                    default=False,
                    dest='error_log',
                    help='Record errors in \'error.log\'')
args = parser.parse_args()
error_log = args.error_log


def error(msg):
    if error_log:
        with open('error.log', 'w') as f:
            f.write(msg + '\n')
    raise RuntimeError(msg)


# Check config directory
dir_path = cwd + "/cfg"
if not os.path.exists(dir_path) or not os.listdir(dir_path):
    error('Error: no config directory')

if not any(map(lambda x: x.endswith('.ini'), os.listdir(dir_path))):
    error('no .ini files in config directory')

# Open config files
userconfig = configparser.ConfigParser()
with codecs.open(cwd + '/cfg/' + args.config_id + '_config.ini', 'r', encoding='utf-8') as f:
    userconfig.read_file(f)
masterconfig = configparser.ConfigParser()
with codecs.open(cwd + '/cfg/master_config.ini', 'r', encoding='utf-8') as f:
    masterconfig.read_file(f)

# Parse
dont_filter_by_sample_id = None
if 'dont_filter_by_sample_id' in userconfig[args.config_section]:
    dont_filter_by_sample_id = userconfig[args.config_section]['dont_filter_by_sample_id']
linked_samples = userconfig[args.config_section]['snaptron_samples'] #samples separated by comma/dash
unlinked_samples = re.sub(r'-', ',', linked_samples) #samples: replace dashes with commas
qs_cov_filt = b'coverage_sum>=%d' % args.coverage_sum
qs_len_filt = b'length<=%d' % args.max_length
rail_ids = unlinked_samples.split(',')
unlinked_samples = bytes(unlinked_samples, 'utf-8')
if dont_filter_by_sample_id == "1":
    sys.stdout.write("not filtering by sample id\n")
    unlinked_samples = bytes("", 'utf-8')


def mkdir_quiet(dr):
    """ Create directories needed to ensure 'dr' exists; no complaining """
    import errno
    if not os.path.isdir(dr):
        try:
            os.makedirs(dr)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


temp_dir = None
if args.query_temp is not None:
    temp_dir = args.query_temp
    mkdir_quiet(temp_dir)
else:
    temp_dir = tempfile.mkdtemp()


def parse_row(row):
    if len(row) != 3:
        error('Malformed row of --input file: %s' % str(row))
    exid, rang, strand = row
    if strand not in '-+.':
        error('Malformed stramd in --input file: %s' % str(row))
    rang_toks = re.split('[:-]', rang)
    if len(rang_toks) != 3:
        error('Malformed range in --input file: %s' % str(row))
    chrom, st, en = rang_toks
    st, en = int(st), int(en)
    assert en >= st
    if chrom not in chromosome_list:
        error('exon_chr not in chromosome_list, ExonID: ' + exid)
    return exid, chrom, st, en, strand


def compose_junction_client_bulk_query(fn):
    with open(fn, 'r') as fh:
        for i, row in enumerate(csv.reader(fh, delimiter='\t', quoting=csv.QUOTE_NONE)):
            rw = parse_row(row)
            if rw is None:
                break
            exid, chrom, st, en, strand = rw
            region = bytes('%s:%d-%d' % (chrom, st-1, en+1), 'utf-8')
            filters = b'&'.join([qs_cov_filt, qs_len_filt, b'strand=' + bytes(strand, 'utf-8')])
            exid, chrom, strand = map(lambda x: bytes(x, 'utf-8'), [exid, chrom, strand])
            yield [exid, chrom, st, en, strand, region, filters, unlinked_samples, i]


def grouper(n, iterable):
    it = iter(iterable)
    while True:
       chunk = list(itertools.islice(it, n))
       if not chunk:
           return
       yield chunk


def junction_query_batches(fn, queries_per_batch, verbose=False):
    query_fn = os.path.join(temp_dir, 'pass1.tsv')
    result_fn = os.path.join(temp_dir, 'pass1.tsv.snap_results.tsv')
    for i, chunk in enumerate(grouper(queries_per_batch, compose_junction_client_bulk_query(fn))):
        assert not os.path.exists(query_fn)
        assert not os.path.exists(result_fn)
        exon_infos = {}
        with open(query_fn, 'wb') as ofh:
            ofh.write(b'\t'.join([b'region', b'filters', b'samples', b'group']) + b'\n')
            for ln in chunk:
                exid, chrom, st, en, strand, region, filters, unlinked_samples, group_id = ln
                exon_infos[group_id] = {'exon_id': exid,
                                       'chromosome': chrom,
                                       'start': st,
                                       'end': en,
                                       'strand': strand}
                ofh.write(b'\t'.join([region, filters, unlinked_samples, bytes(str(group_id), 'utf-8')]) + b'\n')
        cmd = ['python', './software/snaptron/qs', '--bulk-query-file=' + query_fn,
           '--datasrc=' + userconfig[args.config_section]['snaptron_datasrc']]
        cmd = ' '.join(cmd)
        if verbose:
            print(cmd)
        ret = os.system(cmd)
        if ret != 0:
            raise RuntimeError('Snaptron query script returned %d' % ret)
        with open(result_fn, 'rb') as rfh:
            df = pd.read_csv(rfh, sep='\t')
        if not args.keep_intermediates:
            os.remove(query_fn)
            os.remove(result_fn)
        else:
            print (query_fn, temp_dir, i)
            os.rename(query_fn, os.path.join(temp_dir, 'pass1_queries_batch%d.tsv' % (i+1)))
            os.rename(result_fn, os.path.join(temp_dir, 'pass1_responses_batch%d.tsv' % (i+1)))
        yield df, exon_infos


def handle_junction_query(df, exon_info):
    """
    Use results from the first-pass junction query to construct second-pass
    exclusion junction queries.
    """
    assert df.shape[0] > 0
    exon_id, exon_chr, exon_start, exon_end, strand = \
        exon_info['exon_id'], exon_info['chromosome'], exon_info['start'], \
        exon_info['end'], exon_info['strand']

    inclusion_left, inclusion_right = None, None
    EJL_pass = 1  # 0 = EJL fraction/total too small; -1 = empty dataframe
    EJR_pass = 1  # 0 = EJR fraction/total too small; -1 = empty dataframe

    # Slice df_pm1 for exon left/right exon junctions and find index of max coverage_sum
    df_EJL = df.loc[df['end'] == exon_start - 1]
    df_EJR = df.loc[df['start'] == exon_end + 1]
    index_EJL_max, index_EJR_max = None, None
    EJLsum, EJRsum = 0, 0
    if df_EJL.empty and df_EJR.empty:
        return None

    if not df_EJL.empty:
        EJLsum = df_EJL['coverage_sum'].sum()
        index_EJL_max = df_EJL['coverage_sum'].idxmax()
    if not df_EJR.empty:
        EJRsum = df_EJR['coverage_sum'].sum()
        index_EJR_max = df_EJR['coverage_sum'].idxmax()

    # Inclusion/Exclusion Junction Reporting
    inLRstr = ""
    inContinue = 0    
    
    # Find inclusion junctions, left and right of exon
    # Left:    
    if df_EJL.empty:
        EJL_pass = -1
    elif df_EJL.get_value(index_EJL_max, 'coverage_sum') / EJLsum < args.inclusion_fraction:
        EJL_pass = 0
    else:
        inclusion_left = [
            bytes(df_EJL.get_value(index_EJL_max, 'chromosome'), 'utf-8'),
            int(df_EJL.get_value(index_EJL_max, 'start')),
            int(df_EJL.get_value(index_EJL_max, 'end')),
            bytes(df_EJL.get_value(index_EJL_max, 'strand'), 'utf-8')
        ]
        assert inclusion_left[0] == exon_chr

        # Inclusion/Exclusion Junction Reporting
        inLRstr = inLRstr + exon_id.decode("utf-8") + ','
        df_EJL = df_EJL.sort_values('coverage_sum', ascending=False)
        csSum = df_EJL['coverage_sum'].sum()
        inexCounter = 0
        for x in df_EJL.index.values:
            inLRstr = inLRstr + str(df_EJL.get_value(x, 'start')) + ','
            inLRstr = inLRstr + str(df_EJL.get_value(x, 'end')) + ','
            inLRstr = inLRstr + str(int(round(df_EJL.get_value(x, 'coverage_sum')*100/csSum))) + ','
            inexCounter = inexCounter + 1
            if inexCounter >= 3:
                break
        while inexCounter < 3:
            inLRstr = inLRstr + '-1,-1,-1,'
            inexCounter = inexCounter + 1
        inContinue = 1

    # Right:
    if df_EJR.empty:
        EJR_pass = -1
    elif df_EJR.get_value(index_EJR_max, 'coverage_sum') / EJRsum < args.inclusion_fraction:
        EJR_pass = 0
    else:
        inclusion_right = [
            bytes(df_EJR.get_value(index_EJR_max, 'chromosome'), 'utf-8'),
            int(df_EJR.get_value(index_EJR_max, 'start')),
            int(df_EJR.get_value(index_EJR_max, 'end')),
            bytes(df_EJR.get_value(index_EJR_max, 'strand'), 'utf-8')
        ]
        assert inclusion_right[0] == exon_chr
        
        # Inclusion/Exclusion Junction Reporting
        if inContinue == 1:
            df_EJR = df_EJR.sort_values('coverage_sum', ascending=False)
            csSum = df_EJR['coverage_sum'].sum()
            inexCounter = 0
            for x in df_EJR.index.values:
                inLRstr = inLRstr + str(df_EJR.get_value(x, 'start')) + ','
                inLRstr = inLRstr + str(df_EJR.get_value(x, 'end')) + ','
                inLRstr = inLRstr + str(int(round(df_EJR.get_value(x, 'coverage_sum')*100/csSum))) + ','
                inexCounter = inexCounter + 1
                if inexCounter >= 3:
                    break
            while inexCounter < 3:
                inLRstr = inLRstr + '-1,-1,-1,'
                inexCounter = inexCounter + 1
            with open('y_inLR.csv', 'a') as inex:
                inex.write(inLRstr + '\n')

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    # Snaptron query for boundary junctions
    boundary_coordinates = None
    if EJL_pass == 1 and EJR_pass == 1:
        st_str, en_str = bytes(str(inclusion_left[1]), 'utf-8'), bytes(str(inclusion_right[2]), 'utf-8')
        boundary_coordinates = b''.join([exon_chr, b':', st_str, b'-', en_str])
    elif EJL_pass == 1 and EJR_pass != 1:  # Query only single bp of left boundary
        st_str, en_str = bytes(str(inclusion_left[1]), 'utf-8'), bytes(str(inclusion_left[1]), 'utf-8')
        boundary_coordinates = b''.join([exon_chr, b':', st_str, b'-', en_str])
    elif EJL_pass != 1 and EJR_pass == 1:  # Query only single bp of right boundary
        st_str, en_str = bytes(str(inclusion_right[2]), 'utf-8'), bytes(str(inclusion_right[2]), 'utf-8')
        boundary_coordinates = b''.join([exon_chr, b':', st_str, b'-', en_str])

    exon_info['ejl_pass'] = EJL_pass
    exon_info['ejr_pass'] = EJR_pass
    exon_info['inclusion_left'] = inclusion_left
    exon_info['inclusion_right'] = inclusion_right
    return boundary_coordinates


def handle_junction_query_batch(df, exon_infos, batchi, verbose=False):
    query2_fn = os.path.join(temp_dir, 'pass2.tsv')
    result2_fn = os.path.join(temp_dir, 'pass2.tsv.snap_results.tsv')
    with open(query2_fn, 'wb') as ofh:
        ofh.write(b'\t'.join([b'region', b'filters', b'samples', b'group']) + b'\n')
        for group_id in set(df.Group):
            assert group_id in exon_infos
            exon_info = exon_infos[group_id]
            region = handle_junction_query(df.loc[df.Group == group_id], exon_info)
            if region is None:
                continue
            filters = b'&'.join([qs_cov_filt, qs_len_filt, b'strand=' + exon_info['strand']])
            ofh.write(b'\t'.join([region, filters, unlinked_samples, bytes(str(group_id), 'utf-8')]) + b'\n')

    cmd = ['python', './software/snaptron/qs', '--bulk-query-file=' + query2_fn,
           '--datasrc=' + userconfig[args.config_section]['snaptron_datasrc']]
    cmd = ' '.join(cmd)
    if verbose:
        print(cmd)
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError('Snaptron query script returned %d' % ret)
    nlines = len(open(result2_fn).readlines())
    result_df = None
    if nlines > 1:
        result_df = pd.read_csv(result2_fn, sep='\t')
    if not args.keep_intermediates:
        os.remove(query2_fn)
        os.remove(result2_fn)
    else:
        os.rename(query2_fn, os.path.join(temp_dir, 'pass2_queries_batch%d.tsv' % batchi))
        os.rename(result2_fn, os.path.join(temp_dir, 'pass2_responses_batch%d.tsv' % batchi))
    return result_df


# Iterate over batches of query exons
batchi = 0
avg_psis, tot_junc_counts = [], []
for df, exon_infos in junction_query_batches(args.input, args.queries_per_batch, verbose=args.verbose):
    batchi += 1
    # Iterate over Snaptron groups, each corresponding to a single query exon
    df2_full = handle_junction_query_batch(df, exon_infos, batchi, verbose=args.verbose)
    if df2_full is None:
        continue
    for group_id in set(df2_full.Group):

        df2 = df2_full.loc[df2_full.Group == group_id]
        exclusion_left, exclusion_right = [], []
        exon_info = exon_infos[group_id]
        EJL_pass, EJR_pass = exon_info['ejl_pass'], exon_info['ejr_pass']
        inclusion_left, inclusion_right = exon_info['inclusion_left'], exon_info['inclusion_right']
        exon_id = exon_info['exon_id']
        
        # Inclusion/Exclusion Junction Reporting
        exLRstr = ""
        exLRcontinue = 1
        
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Slice df_boundary for boundary left junctions, record exclusion junctions
        if EJL_pass == 1:
            df_BJL = df2[df2['start'] == inclusion_left[1]]
            df_BJL = df_BJL[df_BJL['end'] != inclusion_left[2]]  # remove inclusion_left
            df_BJL = df_BJL.sort_values('coverage_sum', ascending=False).reset_index(drop=True)  # sort and reset index

            if df_BJL.empty:
                # Constitutive Left Junction
                incLstr = ',' + str(inclusion_left[1]) + ',' + str(inclusion_left[2]) + ',999'
                exLRstr = exLRstr + exon_id.decode("utf-8") + incLstr + ',-1,-1,-1,-1,-1,-1,'
            else:
                for x in range(len(df_BJL.index)):
                    exclusion_left.append([
                        df_BJL.get_value(x, 'chromosome'),
                        df_BJL.get_value(x, 'start'),
                        df_BJL.get_value(x, 'end'),
                        df_BJL.get_value(x, 'strand')
                    ])
                
                # Inclusion/Exclusion Junction Reporting
                exLRstr = exLRstr + exon_id.decode("utf-8") + ','
                csSum = df_BJL['coverage_sum'].sum()
                inexCounter = 0
                for x in df_BJL.index.values:
                    exLRstr = exLRstr + str(df_BJL.get_value(x, 'start')) + ','
                    exLRstr = exLRstr + str(df_BJL.get_value(x, 'end')) + ','
                    exLRstr = exLRstr + str(int(round(df_BJL.get_value(x, 'coverage_sum')*100/csSum))) + ','
                    inexCounter = inexCounter + 1
                    if inexCounter >= 3:
                        break
                while inexCounter < 3:
                    exLRstr = exLRstr + '-1,-1,-1,'
                    inexCounter = inexCounter + 1
        else: # i.e. EJL_pass = 0 or -1
            exLRcontinue = 0

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # Slice df_boundary for boundary right junctions, record exclusion junctions
        if EJR_pass == 1:
            df_BJR = df2[df2['end'] == inclusion_right[2]]
            df_BJR = df_BJR[df_BJR['start'] != inclusion_right[1]]  # remove inclusion_right
            df_BJR = df_BJR.sort_values('coverage_sum', ascending=False).reset_index(drop=True)  # sort and reset index

            if exLRcontinue == 1:
                if df_BJR.empty:
                    # Constitutive Right Junction
                    incRstr = str(inclusion_right[1]) + ',' + str(inclusion_right[2]) + ',999'
                    exLRstr = exLRstr + incRstr + ',-1,-1,-1,-1,-1,-1,'
                    with open('y_exLR.csv', 'a') as inex:
                        inex.write(exLRstr + '\n')
                else:
                    for x in range(len(df_BJR.index)):
                        exclusion_right.append([
                            df_BJR.get_value(x, 'chromosome'),
                            df_BJR.get_value(x, 'start'),
                            df_BJR.get_value(x, 'end'),
                            df_BJR.get_value(x, 'strand')
                        ])

                    # Inclusion/Exclusion Junction Reporting
                    csSum = df_BJR['coverage_sum'].sum()
                    inexCounter = 0
                    for x in df_BJR.index.values:
                        exLRstr = exLRstr + str(df_BJR.get_value(x, 'start')) + ','
                        exLRstr = exLRstr + str(df_BJR.get_value(x, 'end')) + ','
                        exLRstr = exLRstr + str(int(round(df_BJR.get_value(x, 'coverage_sum')*100/csSum))) + ','
                        inexCounter = inexCounter + 1
                        if inexCounter >= 3:
                            break
                    while inexCounter < 3:
                        exLRstr = exLRstr + '-1,-1,-1,'
                        inexCounter = inexCounter + 1
                    with open('y_exLR.csv', 'a') as inex:
                        inex.write(exLRstr + '\n')

        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        # B. Calculate PSI and total junction counts
        # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        # Initiate inclusion/exclusion junction dataframes (and b2b dataframe)
        df_LeftInclusion = pd.DataFrame({"rail_id": rail_ids})
        df_RightInclusion = pd.DataFrame({"rail_id": rail_ids})
        df_LeftExclusion = pd.DataFrame({"rail_id": rail_ids})
        df_RightExclusion = pd.DataFrame({"rail_id": rail_ids})
        df_b2b = pd.DataFrame({"rail_id": rail_ids})

        # Fill out boundary-to-boundary (b2b) dataframe
        # Used to remove duplicate counting of b2b junction in TJC
        if EJL_pass == 1 and EJR_pass == 1:
            match_start = df2['start'] == inclusion_left[1]
            match_end = df2['end'] == inclusion_right[2]
            df_match = df2[match_start & match_end].reset_index(drop=True)
            dict_samplecounts = {}
            if df_match.shape[0] > 0:
                dict_str = df_match['samples'][0]
                if dict_str[0] == ',':
                    dict_str = dict_str[1:]
                dict_samplecounts = {x.split(':')[0]: x.split(':')[1] for x in dict_str.split(',')}
            df_b2b['b2b'] = df_b2b["rail_id"].map(dict_samplecounts.get).fillna(value=0)
            df_b2b = df_b2b.apply(pd.to_numeric).drop('rail_id', 1)
            df_b2b = df_b2b.sum(axis=1)

        # Fill out df_LeftInclusion
        if EJL_pass == 1:
            match_start = df2['start'] == inclusion_left[1]
            match_end = df2['end'] == inclusion_left[2]
            df_match = df2[match_start & match_end].reset_index(drop=True)
            dict_str = df_match['samples'][0]
            if dict_str[0] == ',':
                dict_str = dict_str[1:]
            dict_samplecounts = {x.split(':')[0]: x.split(':')[1] for x in dict_str.split(',')}
            df_LeftInclusion[exon_id] = df_LeftInclusion["rail_id"].map(dict_samplecounts.get).fillna(value=0)
            df_LeftInclusion = df_LeftInclusion.apply(pd.to_numeric).drop('rail_id', 1)
            df_LeftInclusion = df_LeftInclusion.sum(axis=1)

        # Fill out df_RightInclusion
        if EJR_pass == 1:
            match_start = df2['start'] == inclusion_right[1]
            match_end = df2['end'] == inclusion_right[2]
            df_match = df2[match_start & match_end].reset_index(drop=True)
            dict_str = df_match['samples'][0]
            if dict_str[0] == ',':
                dict_str = dict_str[1:]
            dict_samplecounts = {x.split(':')[0]: x.split(':')[1] for x in dict_str.split(',')}
            df_RightInclusion[exon_id] = df_RightInclusion["rail_id"].map(dict_samplecounts.get).fillna(value=0)
            df_RightInclusion = df_RightInclusion.apply(pd.to_numeric).drop('rail_id', 1)
            df_RightInclusion = df_RightInclusion.sum(axis=1)

        # Fill out df_LeftExclusion
        if EJL_pass == 1:
            for x in range(len(exclusion_left)):
                match_start = df2['start'] == exclusion_left[x][1]
                match_end = df2['end'] == exclusion_left[x][2]
                df_match = df2[match_start & match_end].reset_index(drop=True)
                dict_str = df_match['samples'][0]
                if dict_str[0] == ',':
                    dict_str = dict_str[1:]
                dict_samplecounts = {x.split(':')[0]: x.split(':')[1] for x in dict_str.split(',')}
                df_LeftExclusion[exon_id + b'_' + bytes(str(x), 'utf-8')] = df_LeftExclusion["rail_id"].map(dict_samplecounts.get)

            df_LeftExclusion = df_LeftExclusion.fillna(value=0).apply(pd.to_numeric)
            df_tempL = df_LeftExclusion
            df_LeftExclusion = df_LeftExclusion.drop('rail_id', 1)
            df_LeftExclusion = df_LeftExclusion.sum(axis=1)

        # Fill out df_RightExclusion
        if EJR_pass == 1:
            for x in range(len(exclusion_right)):
                match_start = df2['start'] == exclusion_right[x][1]
                match_end = df2['end'] == exclusion_right[x][2]
                df_match = df2[match_start & match_end].reset_index(drop=True)
                dict_str = df_match['samples'][0]
                if dict_str[0] == ',':
                    dict_str = dict_str[1:]
                dict_samplecounts = {x.split(':')[0]: x.split(':')[1] for x in dict_str.split(',')}
                df_RightExclusion[exon_id + b'_' + bytes(str(x), 'utf-8')] = df_RightExclusion["rail_id"].map(dict_samplecounts.get)

            df_RightExclusion = df_RightExclusion.fillna(value=0).apply(pd.to_numeric)
            df_RightExclusion = df_RightExclusion.drop('rail_id', 1)
            df_RightExclusion = df_RightExclusion.sum(axis=1)

        # Calculate PSI [Left/Right/Avg] and JunctionCount [Left/Right/Total]
        if EJL_pass == 1 and EJR_pass == 1:
            df_LeftJunctionCount = df_LeftInclusion.add(df_LeftExclusion).astype(int)
            df_RightJunctionCount = df_RightInclusion.add(df_RightExclusion).astype(int)
            df_TotalJunctionCount = df_LeftJunctionCount.add(df_RightJunctionCount).subtract(df_b2b).astype(int) #.subtract(df_b2b) to remove boundary to boundary junction
            df_LeftPSI = df_LeftInclusion.div(df_LeftJunctionCount)
            df_RightPSI = df_RightInclusion.div(df_RightJunctionCount)
            df_AvgPSI = df_LeftPSI.add(df_RightPSI).multiply(0.5)

            # Link '-' delimited samples together
            df_PSIandTJC = pd.DataFrame({"rail_id": rail_ids})
            df_PSIandTJC['LPSI'] = df_LeftPSI
            df_PSIandTJC['RPSI'] = df_RightPSI
            df_PSIandTJC['APSI'] = df_AvgPSI
            df_PSIandTJC['LJC'] = df_LeftJunctionCount
            df_PSIandTJC['RJC'] = df_RightJunctionCount
            df_PSIandTJC['TJC'] = df_TotalJunctionCount
            df_PSIandTJC = df_PSIandTJC.set_index('rail_id')
            LinkedSamples = linked_samples.split(',')
            LPSI_newcol = []
            RPSI_newcol = []
            APSI_newcol = []
            LJC_newcol = []
            RJC_newcol = []
            TJC_newcol = []
            for x in range(len(LinkedSamples)):
                linked_rail_ids = LinkedSamples[x].split('-')
                if len(linked_rail_ids) > 1:
                    totalLPSI = 0
                    totalRPSI = 0
                    totalAPSI = 0
                    totalLJC = 0
                    totalRJC = 0
                    totalTJC = 0
                    
                    for y in linked_rail_ids:
                        addLPSI = df_PSIandTJC.get_value(y, 'LPSI')
                        addRPSI = df_PSIandTJC.get_value(y, 'RPSI')
                        addAPSI = df_PSIandTJC.get_value(y, 'APSI')
                        addLJC = df_PSIandTJC.get_value(y, 'LJC')
                        addRJC = df_PSIandTJC.get_value(y, 'RJC')
                        addTJC = df_PSIandTJC.get_value(y, 'TJC')
                        
                        #If statements handle if rail_id occurs multiple times in snaptron_samples
                        if type(addLPSI) is np.ndarray:
                            addLPSI = addLPSI[0]
                        if type(addRPSI) is np.ndarray:
                            addRPSI = addRPSI[0]
                        if type(addAPSI) is np.ndarray:
                            addAPSI = addAPSI[0]
                        if type(addLJC) is np.ndarray:
                            addLJC = addLJC[0]
                        if type(addRJC) is np.ndarray:
                            addRJC = addRJC[0]
                        if type(addTJC) is np.ndarray:
                            addTJC = addTJC[0]

                        if math.isnan(addLPSI):
                            addLPSI = 0
                        if math.isnan(addRPSI):
                            addRPSI = 0
                        if math.isnan(addAPSI):
                            addAPSI = 0
                        
                        totalLPSI = totalLPSI + addLPSI * addLJC
                        totalRPSI = totalRPSI + addRPSI * addRJC
                        totalAPSI = totalAPSI + addAPSI * addTJC
                        totalLJC = totalLJC + addLJC
                        totalRJC = totalRJC + addRJC
                        totalTJC = totalTJC + addTJC
                    
                    if totalLJC == 0:
                        totalLPSI = 0
                    else:
                        totalLPSI = totalLPSI / totalLJC
                    if totalRJC == 0:
                        totalRPSI = 0
                    else:
                        totalRPSI = totalRPSI / totalRJC
                    if totalTJC == 0:
                        totalAPSI = 0
                    else:
                        totalAPSI = totalAPSI / totalTJC
                    
                    LPSI_newcol.extend([totalLPSI])
                    RPSI_newcol.extend([totalRPSI])
                    APSI_newcol.extend([totalAPSI])
                    LJC_newcol.extend([totalLJC])
                    RJC_newcol.extend([totalRJC])
                    TJC_newcol.extend([totalTJC])
                else:
                    addLPSI = df_PSIandTJC.get_value(linked_rail_ids[0], 'LPSI')
                    addRPSI = df_PSIandTJC.get_value(linked_rail_ids[0], 'RPSI')
                    addAPSI = df_PSIandTJC.get_value(linked_rail_ids[0], 'APSI')
                    addLJC = df_PSIandTJC.get_value(linked_rail_ids[0], 'LJC')
                    addRJC = df_PSIandTJC.get_value(linked_rail_ids[0], 'RJC')                    
                    addTJC = df_PSIandTJC.get_value(linked_rail_ids[0], 'TJC')
                    if type(addLPSI) is np.ndarray:
                        addLPSI = addLPSI[0]
                    if type(addRPSI) is np.ndarray:
                        addRPSI = addRPSI[0]
                    if type(addAPSI) is np.ndarray:
                        addAPSI = addAPSI[0]
                    if type(addLJC) is np.ndarray:
                        addLJC = addLJC[0]
                    if type(addRJC) is np.ndarray:
                        addRJC = addRJC[0]
                    if type(addTJC) is np.ndarray:
                        addTJC = addTJC[0]                        
                    LPSI_newcol.extend([addLPSI])
                    RPSI_newcol.extend([addRPSI])
                    APSI_newcol.extend([addAPSI])
                    LJC_newcol.extend([addLJC])
                    RJC_newcol.extend([addRJC])
                    TJC_newcol.extend([addTJC])
                    
            df_LeftPSI = pd.Series(LPSI_newcol)
            df_RightPSI = pd.Series(RPSI_newcol)
            df_AvgPSI = pd.Series(APSI_newcol)
            df_LeftJunctionCount = pd.Series(LJC_newcol)
            df_RightJunctionCount = pd.Series(RJC_newcol)
            df_TotalJunctionCount = pd.Series(TJC_newcol)

            # Convert to csv
            LeftJunctionCount = df_LeftJunctionCount.to_frame().T.to_csv(index=False, header=False)
            RightJunctionCount = df_RightJunctionCount.to_frame().T.to_csv(index=False, header=False)
            TotalJunctionCount = df_TotalJunctionCount.to_frame().T.to_csv(index=False, header=False)

            LeftPSI = df_LeftPSI.fillna(value=0).multiply(100).round(3).astype(float).to_frame().T.to_csv(index=False, header=False)
            RightPSI = df_RightPSI.fillna(value=0).multiply(100).round(3).astype(float).to_frame().T.to_csv(index=False, header=False)
            AvgPSI = df_AvgPSI.fillna(value=0).multiply(100).round(3).astype(float).to_frame().T.to_csv(index=False, header=False)
        elif EJL_pass != 1 and EJR_pass == 1:
            df_RightJunctionCount = df_RightInclusion.add(df_RightExclusion).astype(int)
            df_RightPSI = df_RightInclusion.div(df_RightJunctionCount)
            RightJunctionCount = df_RightJunctionCount.to_frame().T.to_csv(index=False, header=False)
            RightPSI = df_RightPSI.fillna(value=0).multiply(100).round(3).astype(float).to_frame().T.to_csv(index=False, header=False)
            if EJL_pass == 0:
                LeftJunctionCount = 'ExonJunctionLeft<args.inclusion_fraction\n'
                TotalJunctionCount = 'ExonJunctionLeft<args.inclusion_fraction\n'
                LeftPSI = 'ExonJunctionLeft<args.inclusion_fraction\n'
                AvgPSI = 'ExonJunctionLeft<args.inclusion_fraction\n'
            elif EJL_pass == -1:
                LeftJunctionCount = 'ExonJunctionLeftEmptyDataframe\n'
                TotalJunctionCount = 'ExonJunctionLeftEmptyDataframe\n'
                LeftPSI = 'ExonJunctionLeftEmptyDataframe\n'
                AvgPSI = 'ExonJunctionLeftEmptyDataframe\n'
        elif EJL_pass == 1 and EJR_pass != 1:
            df_LeftJunctionCount = df_LeftInclusion.add(df_LeftExclusion).astype(int)
            df_LeftPSI = df_LeftInclusion.div(df_LeftJunctionCount)
            LeftJunctionCount = df_LeftJunctionCount.to_frame().T.to_csv(index=False, header=False)
            LeftPSI = df_LeftPSI.fillna(value=0).multiply(100).round(3).astype(float).to_frame().T.to_csv(index=False, header=False)
            if EJR_pass == 0:
                RightJunctionCount = 'ExonJunctionRight<args.inclusion_fraction\n'
                TotalJunctionCount = 'ExonJunctionRight<args.inclusion_fraction\n'
                RightPSI = 'ExonJunctionRight<args.inclusion_fraction\n'
                AvgPSI = 'ExonJunctionRight<args.inclusion_fraction\n'
            elif EJR_pass == -1:
                RightJunctionCount = 'ExonJunctionRightEmptyDataframe\n'
                TotalJunctionCount = 'ExonJunctionRightEmptyDataframe\n'
                RightPSI = 'ExonJunctionRightEmptyDataframe\n'
                AvgPSI = 'ExonJunctionRightEmptyDataframe\n'

        avg_psis.append([exon_id, b'LeftPSI', bytes(LeftPSI.rstrip(), 'utf-8')])
        avg_psis.append([exon_id, b'RightPSI', bytes(RightPSI.rstrip(), 'utf-8')])
        avg_psis.append([exon_id, b'AvgPSI', bytes(AvgPSI.rstrip(), 'utf-8')])
        tot_junc_counts.append([exon_id, b'LeftJunctionCount', bytes(LeftJunctionCount.rstrip(), 'utf-8')])
        tot_junc_counts.append([exon_id, b'RightJunctionCount', bytes(RightJunctionCount.rstrip(), 'utf-8')])
        tot_junc_counts.append([exon_id, b'TotalJunctionCount', bytes(TotalJunctionCount.rstrip(), 'utf-8')])

if not os.path.exists('bulkoutput.csv'):
    f = open('bulkoutput.csv','w')
    f.close()

with open('bulkoutput.csv', 'a') as f:
    for avg_psi_line in avg_psis:
        f.write(b','.join(avg_psi_line).decode('utf-8') + '\n')
    for tot_junc_count_line in tot_junc_counts:
        f.write(b','.join(tot_junc_count_line).decode('utf-8') + '\n')
