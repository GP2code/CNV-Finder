from CNV_finder.data_methods import check_interval, create_test_set, fill_window_df
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Arguments for Running CNV-Finder ML Data Prep.')    
parser.add_argument('--interval_name', type=str, default=None, help='Name for NDD-related Gene region.')
parser.add_argument('--interval_file', type=str, default='ref_files/glist_hg38_intervals.csv', help='Gene or other feature intervals to analyze. Header is [NAME,CHR,START,STOP], one line per interval. Autosomes only.')
parser.add_argument('--chr', type=str, default=None, help='Chromsome for region of interest.')
parser.add_argument('--start', type=int, default=None, help='Starting position in base pairs (hg38).')
parser.add_argument('--stop', type=int, default=None, help='Stopping position in base pairs (hg38).')
parser.add_argument('--buffer', type=int, default=250000, help='Kilobase window around each interval, in bases.')
parser.add_argument('--min_gentrain', type=float, default=0.2, help='Minimum GenTrain Score threshold.')
parser.add_argument('--split_interval', type=int, default=5, help='Number of intervals to split region into (before overlap).')
parser.add_argument('--total_windows', type=int, default=31, help='Count of windows in region per sample (with overlap).')
parser.add_argument('--bim_file', type=str, default=None, help='PLINK .bim file following sample QC.')
parser.add_argument('--pvar_file', type=str, default=None, help='PLINK2 .pvar file following sample QC.')
parser.add_argument('--master_file', type=str, default='ref_files/master_key.txt', help='Master key for all available samples in data release (need IID and binary GDPR (0-no, 1-yes) column).')
parser.add_argument('--study_name', type=str, default='all', help='Subset or cohort of interest from larger data release.')
parser.add_argument('--metrics_path', type=str, default='ref_files/snp_metrics', help='Path to SNP metrics files/parquets with the format {metrics_path}/{sentrixbarcode}/snp_metrics_{sentrixbarcode}/Sample_ID={sample}.')
parser.add_argument('--out_path', type=str, default=None, help='Path to output reports with suggested format {Cohort}_{Gene or Interval Name}.')

# need to figure out training set creation
parser.add_argument('--training_ids', type=str, default=None, help='List of IDs used in training set with header IID.')
parser.add_argument('--testing_ids', type=str, default=None, help='List of IDs used in testing set with header IID.')
parser.add_argument('--create_training', action='store_true', help='Create new training set (functionality coming soon).')
parser.add_argument('--train_size', type=int, default=100, help='Number of samples to include in training set.')
parser.add_argument('--create_testing', action='store_true', help='Create new testing set with no overlaps in training set IIDs.')
parser.add_argument('--test_size', type=int, default=500, help='Number of samples to include in testing set')

args = parser.parse_args()

interval_name = args.interval_name
interval_file = args.interval_file
chr = args.chr
start_pos = args.start
stop_pos = args.stop
buffer = args.buffer
min_gentrain = args.min_gentrain
split_interval = args.split_interval
total_windows = args.total_windows
bim = args.bim_file
pvar = args.pvar_file
master_file = args.master_file
study_name = args.study_name
snp_metrics_path = args.metrics_path
out_path = args.out_path
create_train = args.create_training
create_test = args.create_testing
train_size = args.train_size
test_size = args.test_size
train_ids = args.training_ids
test_ids = args.testing_ids

# if submitted interval in interval file, find positions
# if empty dataframe, request new name/manual chr/positions
if interval_name:
    positions = check_interval(interval_name, interval_file)
    if len(positions) > 0:
        print(positions)
        chr = positions.CHR.values[0]
        start_pos = positions.START.values[0]
        stop_pos = positions.STOP.values[0]
    else:
        print('Interval name not found in interval reference file. Please enter a new interval name or manually enter chromosome with start and stop base pair positions for interval of interest.')

if create_test:
    test_set = True
    create_test_set(master_file, test_size, train_ids, snp_metrics_path, out_path, study_name)
    fill_window_df(f'{out_path}_testing_IDs.csv', chr, start_pos, stop_pos, out_path, split_interval, total_windows, buffer, test_set, min_gentrain, bim, pvar)
    
elif test_ids:
    test_set = True
    fill_window_df(test_ids, chr, start_pos, stop_pos, out_path, split_interval, total_windows, buffer, test_set, min_gentrain, bim, pvar)