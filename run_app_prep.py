from locale import strcoll
from CNV_finder.data_methods import check_interval, create_app_ready_file, generate_pred_cnvs
import multiprocessing
import pandas as pd
import argparse


def main():
    parser = argparse.ArgumentParser(description='Arguments for Running CNV-Finder ML Data Prep.')
    parser.add_argument('--interval_name', type=str, default=None, help='Name for NDD-related Gene region.')
    parser.add_argument('--interval_file', type=str, default='ref_files/glist_hg38_intervals.csv', help='Gene or other feature intervals to analyze. Header is [NAME,CHR,START,STOP], one line per interval. Autosomes only.')
    parser.add_argument('--chr', type=str, default=None, help='Chromsome for region of interest.')
    parser.add_argument('--start', type=int, default=None, help='Starting position in base pairs (hg38).')
    parser.add_argument('--stop', type=int, default=None, help='Stopping position in base pairs (hg38).')
    parser.add_argument('--buffer', type=int, default=250000, help='Kilobase window around each interval, in bases.')
    parser.add_argument('--min_gentrain', type=float, default=0.2, help='Minimum GenTrain Score threshold.')
    parser.add_argument('--bim_file', type=str, default=None, help='PLINK .bim file following sample QC.')
    parser.add_argument('--pvar_file', type=str, default=None, help='PLINK2 .pvar file following sample QC.')
    parser.add_argument('--test_set_ids', type=str, default=None, help='Path to file with testing set IIDs and SNP metrics path (headers: IID, snp_metrics_path).')
    parser.add_argument('--test_set_windows', type=str, default=None, help='Path to file with windowed test set used in training.')
    parser.add_argument('--test_set_results', type=str, default=None, help='Path to file with model results on test set.')
    parser.add_argument('--probability', type=float, default=0.8, help='Probability threshold of model predictions that must be met to be included in app.')
    parser.add_argument('--out_path', type=str, default=None, help='Path to output app-ready report with suggested format Cohort_Gene or Interval Name.')
    parser.add_argument('--make_app_ready', action='store_true', help='Create 1 app ready file including all testing samples with necessary info for app creation.')

    args = parser.parse_args()

    interval_name = args.interval_name
    interval_file = args.interval_file
    chr = args.chr
    start_pos = args.start
    stop_pos = args.stop
    buffer = args.buffer
    min_gentrain = args.min_gentrain
    bim = args.bim_file
    pvar = args.pvar_file
    test_set_ids = args.test_set_ids
    test_set_windows = args.test_set_windows
    test_set_results = args.test_set_results
    probability = args.probability
    out_path = args.out_path
    app_ready = args.make_app_ready

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

    if app_ready:
        above_probab = create_app_ready_file(test_set_ids, test_set_windows, test_set_results, out_path, probability)

        # non-parallelized method
        # for metrics in above_probab.snp_metrics_path:
            # generate_pred_cnvs(metrics, chr, start_pos, stop_pos, out_path, buffer, min_gentrain, bim, pvar)

        with multiprocessing.Pool() as pool:
            pool.map(generate_pred_cnvs, [(above_probab.snp_metrics_path, chr, start_pos, stop_pos, out_path, buffer, min_gentrain, bim, pvar) for index, row in above_probab.iterrows()])

if __name__ == "__main__":
    main()