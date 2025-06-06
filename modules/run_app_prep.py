import argparse
import multiprocessing
import pandas as pd
from cnv_finder.data_methods import check_interval, subset_metadata, create_app_ready_file, generate_pred_cnvs


def main():
    parser = argparse.ArgumentParser(
        description='Arguments for Running CNV-Finder App-Ready Data Prep.')
    parser.add_argument('--interval_name', type=str,
                        default=None, help='Name for NDD-related Gene region.')
    parser.add_argument('--interval_file', type=str, default='ref_files/glist_hg38_intervals.csv',
                        help='Gene or other feature intervals to analyze. Header is [NAME,CHR,START,STOP], one line per interval. Autosomes only.')
    parser.add_argument('--chrom', type=str, default=None,
                        help='Chromsome for region of interest.')
    parser.add_argument('--start', type=int, default=None,
                        help='Starting position in base pairs (hg38).')
    parser.add_argument('--stop', type=int, default=None,
                        help='Stopping position in base pairs (hg38).')
    parser.add_argument('--buffer', type=int, default=250000,
                        help='Kilobase window around each interval, in bases.')
    parser.add_argument('--min_gentrain', type=float,
                        default=0.2, help='Minimum GenTrain Score threshold.')
    parser.add_argument('--metadata_path', type=str, default='ref_files/NBA_metadata', help='Repeated SNP metadata')
    parser.add_argument('--bim_file', type=str, default=None,
                        help='PLINK .bim file following sample QC.')
    parser.add_argument('--pvar_file', type=str, default=None,
                        help='PLINK2 .pvar file following sample QC.')
    parser.add_argument('--test_set_ids', type=str, default=None,
                        help='Path to file with testing set IIDs and SNP metrics path (headers: IID, snp_metrics_path).')
    parser.add_argument('--test_set_windows', type=str, default=None,
                        help='Path to file with windowed test set used in training.')
    parser.add_argument('--test_set_results', type=str, default=None,
                        help='Path to file with model results on test set.')
    parser.add_argument('--probability', type=float, default=0.8,
                        help='Probability threshold of model predictions that must be met to be included in app.')
    parser.add_argument('--out_path', type=str, default=None,
                        help='Path to output app-ready report with suggested format Cohort_Gene or Interval Name.')
    parser.add_argument('--cpus', type=int, default=8,
                        help='Number of CPUs available for the job.')
    parser.add_argument('--make_app_ready', action='store_true',
                        help='Create 1 app ready file including all testing samples with necessary info for app creation.')

    # Define variables from argument flags
    args = parser.parse_args()

    interval_name = args.interval_name
    interval_file = args.interval_file
    cpus = args.cpus
    chrom = args.chrom
    start_pos = args.start
    stop_pos = args.stop
    buffer = args.buffer
    min_gentrain = args.min_gentrain
    metadata_path = args.metadata_path
    bim = args.bim_file
    pvar = args.pvar_file
    test_set_ids = args.test_set_ids
    test_set_windows = args.test_set_windows
    test_set_results = args.test_set_results
    probability = args.probability
    out_path = args.out_path
    app_ready = args.make_app_ready

    # Finds the chromosome, start, and stop positions for a submitted interval name
    if interval_name:
        chrom, start_pos, stop_pos = check_interval(
            interval_name, interval_file)
        if not chrom or not start_pos or not stop_pos:
            print('Interval name not found in interval reference file. Please enter a new interval name or manually enter chromosome with start and stop base pair positions for interval of interest.')

    # Subsets metadata file for relevant info
    if metadata_path:
        snp_info = subset_metadata(metadata_path, chrom, start_pos, stop_pos, buffer, min_gentrain)

    # Prepares app-ready files for samples with predicted values above a specified threshold
    if app_ready:
        above_probab = create_app_ready_file(
            test_set_ids, test_set_windows, test_set_results, out_path, probability)

        # Parallelizes the file creation process
        with multiprocessing.Pool(cpus) as pool:
            pool.map(generate_pred_cnvs, [(row.IID, row.snp_metrics_path, snp_info, chrom, start_pos, stop_pos,
                     out_path, buffer, min_gentrain, bim, pvar) for index, row in above_probab.iterrows()])


if __name__ == "__main__":
    main()
