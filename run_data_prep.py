from CNV_finder.data_methods import check_interval, create_test_set, make_window_df, fill_window_df
import pandas as pd
import numpy as np
import argparse
import multiprocessing


def main():
    parser = argparse.ArgumentParser(description='Arguments for Running CNV-Finder ML Data Prep.')    
    parser.add_argument('--interval_name', type=str, default=None, help='Name for NDD-related Gene region.')
    parser.add_argument('--interval_file', type=str, default='ref_files/glist_hg38_intervals.csv', help='Gene or other \
                        feature intervals to analyze. Header is [NAME,CHR,START,STOP], one line per interval. Autosomes only.')
    parser.add_argument('--chr', type=str, default=None, help='Chromsome for region of interest.')
    parser.add_argument('--start', type=int, default=None, help='Starting position in base pairs (hg38).')
    parser.add_argument('--stop', type=int, default=None, help='Stopping position in base pairs (hg38).')
    parser.add_argument('--buffer', type=int, default=250000, help='Kilobase window around each interval, in bases.')
    parser.add_argument('--min_gentrain', type=float, default=0.2, help='Minimum GenTrain Score threshold.')
    parser.add_argument('--split_interval', type=int, default=5, help='Number of intervals to split region into (before overlap).')
    parser.add_argument('--total_windows', type=int, default=31, help='Count of windows in region per sample (with overlap).')
    parser.add_argument('--bim_file', type=str, default=None, help='PLINK .bim file following sample QC.')
    parser.add_argument('--pvar_file', type=str, default=None, help='PLINK2 .pvar file following sample QC.')
    parser.add_argument('--master_file', type=str, default='ref_files/master_key.txt', help='Master key for all available samples in data \
                        release (need IID and binary GDPR (0-no, 1-yes) column).')
    parser.add_argument('--study_name', type=str, default='all', help='Subset or cohort of interest from larger data release.')
    parser.add_argument('--metrics_path', type=str, default='ref_files/snp_metrics', help='Path to SNP metrics files/parquets with the format \
                         {metrics_path}/{sentrixbarcode}/snp_metrics_{sentrixbarcode}/Sample_ID={sample}.')
    parser.add_argument('--out_path', type=str, default=None, help='Path to output reports with suggested format {Cohort}_{Gene or Interval Name}.')
    parser.add_argument('--cpus', type=int, default=8, help='Number of CPUs available for the job.')

    # need to figure out training set creation
    parser.add_argument('--training_ids', type=str, default=None, help='List of IDs used in training set with headers IID. May also include the \
                        interval where CNV was found in each sample. Otherwise interval/chromosome with position range must be included in proper labels.')
    parser.add_argument('--testing_ids', type=str, default=None, help='List of IDs used in testing set with header IID.')
    parser.add_argument('--create_training', action='store_true', help='Create new training set (functionality coming soon).')
    parser.add_argument('--train_size', type=int, default=100, help='Number of samples to include in training set.')
    parser.add_argument('--create_testing', action='store_true', help='Create new testing set with no overlaps in training set IIDs.')
    parser.add_argument('--test_size', type=int, default=500, help='Number of samples to include in testing set')

    args = parser.parse_args()

    interval_name = args.interval_name
    interval_file = args.interval_file
    cpus = args.cpus
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
    train_df = args.training_ids
    test_df = args.testing_ids

    # if submitted interval in interval file, find positions
    # if empty dataframe, request new name/manual entry of chr/positions
    if interval_name:
        positions = check_interval(interval_name, interval_file)
        if len(positions) > 0:
            print(positions)
            chr = positions.CHR.values[0]
            start_pos = positions.START.values[0]
            stop_pos = positions.STOP.values[0]
        elif not chr or not start_pos or not stop_pos:
            print('Interval name not found in interval reference file. Please enter a new interval name or manually enter chromosome with start and stop base pair positions for interval of interest.')
        else:
            print('Interval name not found in interval reference file. Added to custom reference file with base pair positions you provided.')
            new_intervals = pd.DataFrame({'NAME': chr,'CHR': chr,'START': start_pos,'STOP': stop_pos})
            new_intervals.to_csv('ref_files/custom_intervals.csv', mode='a')

    if create_test:
        cnv_exists = np.nan
        
        if not test_df:
            create_test_set(master_file, test_size, train_df, snp_metrics_path, out_path, study_name)

            # Read in testing IDs
            test_df = pd.read_csv(f'{out_path}_testing_IDs.csv')

        # Create df to hold windows that span interval of interest
        window_df = make_window_df(chr, start_pos, stop_pos, split_interval, total_windows, buffer)
        all_samples = pd.DataFrame(columns = ['START', 'STOP', 'dosage_interval', 'dosage_gene', 'del_dosage', 'dup_dosage', 'ins_dosage', 'avg_baf', 'avg_lrr', \
                                             'std_baf', 'std_lrr', 'iqr_baf', 'iqr_lrr', 'cnv_range_count', 'IID', 'CHR', 'window', 'CNV_exists'])
        all_samples.to_csv(f'{out_path}_samples_windows.csv', index = False)

        # Parallelize creation of df that holds all samples with aggregated feature in each window
        with multiprocessing.Pool(cpus) as pool:
            pool.map(fill_window_df, [(out_path, row.IID, row.snp_metrics_path, window_df, cnv_exists, chr, start_pos, stop_pos, buffer, min_gentrain, bim, pvar) for index, row in test_df.iterrows()])
        

    if create_train:
        train_df = pd.read_csv(train_df)

        for sample in train_df.IID:
            code = sample.split('_')[0]
            train_df['snp_metrics_path'] = f'{snp_metrics_path}/{code}/snp_metrics_{code}/Sample_ID={sample}'

        if "Interval" in train_df.columns:
            intervals = np.unique(train_df.Interval)
        elif interval_name:
            intervals = [interval_name]

        for interval_name in intervals:
            positions = check_interval(interval_name, interval_file)
            if len(positions) == 0:
                positions = check_interval(interval_name, 'ref_files/custom_intervals.csv')

                if len(positions) == 0 and chr and start_pos and stop_pos:
                    print('Interval name not found in interval reference file. Added to custom reference file with base pair positions you provided.')
                    new_intervals = pd.DataFrame({'NAME': chr,'CHR': chr,'START': start_pos,'STOP': stop_pos})
                    new_intervals.to_csv('ref_files/custom_intervals.csv', mode='a')
                elif len(positions) == 0:
                    print('Interval name not found in interval reference file. Please enter a new interval name or manually enter chromosome with start and stop base pair positions for interval of interest.')
                    break
            else:
                print(positions)
                chr = positions.CHR.values[0]
                start_pos = positions.START.values[0]
                stop_pos = positions.STOP.values[0]
                
            # Create df to hold windows that span interval of interest
            window_df = make_window_df(chr, start_pos, stop_pos, split_interval, total_windows, buffer)
        window_df = make_window_df(chr, start_pos, stop_pos, split_interval, total_windows, buffer)
        all_samples = pd.DataFrame(columns = ['START', 'STOP', 'dosage_interval', 'dosage_gene', 'del_dosage', 'dup_dosage', 'ins_dosage', 'avg_baf', 'avg_lrr', \
                                             'std_baf', 'std_lrr', 'iqr_baf', 'iqr_lrr', 'cnv_range_count', 'IID', 'CHR', 'window', 'CNV_exists'])
        all_samples.to_csv(f'{out_path}_samples_windows.csv', index = False)

        # Parallelize creation of df that holds all samples with aggregated feature in each window
        with multiprocessing.Pool(cpus) as pool:
            pool.map(fill_window_df, [(out_path, row.IID, row.snp_metrics_path, window_df, cnv_exists, chr, start_pos, stop_pos, buffer, min_gentrain, bim, pvar) for index, row in train_df.iterrows()])
        

if __name__ == "__main__":
    main()