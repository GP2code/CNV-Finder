import os
import shutil
import random
import pandas as pd
import numpy as np
import plotly.express as px
from scipy.stats import iqr


# Supress Pandas copy warning
pd.options.mode.chained_assignment = None


def check_interval(interval_name, interval_file='ref_files/glist_hg38_intervals.csv'):
    """
    Checks for the chromosome, hg38 start, and stop positions for a given gene name in a reference file.
    If the interval is not found, it adds a new entry to a custom reference file for future use.

    Arguments:
    interval_name (str): The name of the gene to look for. 
                         Handles alternate gene names, such as mapping 'PRKN' to 'PARK2'.
    interval_file (str): Path to the interval reference CSV file. Defaults to 
                         'ref_files/glist_hg38_intervals.csv'.

    Returns:
    tuple: A tuple containing:
           - chrom (str or None): The chromosome identifier if found, otherwise None.
           - start_pos (int or None): The starting base pair position if found, otherwise None.
           - stop_pos (int or None): The stopping base pair position if found, otherwise None.
    """
    interval_dir = os.path.dirname(os.path.abspath(interval_file))

    # Handle alternate names for genes, e.g., PRKN and PARK2 - can expand to other common NDD genes
    if interval_name == 'PRKN':
        interval_name = 'PARK2'

    interval_df = pd.read_csv(interval_file)
    positions = interval_df[interval_df.NAME == interval_name]
    if len(positions) > 0:
        chrom = positions.CHR.values[0]
        start_pos = positions.START.values[0]
        stop_pos = positions.STOP.values[0]
    else:
        # Add extracted chromosome, start, and stop positions to custom_interval file for future reference if not found
        chrom, start_pos, stop_pos = None, None, None
        print('Interval name not found in interval reference file. Added to custom reference file with base pair positions you provided.')
        new_intervals = pd.DataFrame(
            {'NAME': interval_name, 'CHR': chrom, 'START': start_pos, 'STOP': stop_pos})
        new_intervals.to_csv(f'{interval_dir}/custom_intervals.csv', mode='a')

    return chrom, start_pos, stop_pos

def merge_samples(tmp_dir, out_path):
    # Collect all file paths ending in .csv
    csv_files = [os.path.join(tmp_dir, f) for f in os.listdir(tmp_dir) if f.endswith(".csv")]
    
    # Read and concatenate
    df_all = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=False)
    
    # Save to output path
    df_all.to_csv(f'{out_path}_samples_windows.csv', index=False)

    # Remove tmp directory
    shutil.rmtree(tmp_dir)

def subset_metadata(metadata_path, chrom, start, stop, buffer, min_gentrain=0.2):
    ### Add functionality to input PLINK files for QC
    metadata = pd.read_parquet(f'{metadata_path}/CHROM={chrom}')

    snp_info = metadata[['snpID', 'POS', 'GenTrain_Score']][(
        metadata['POS'] >= (start - buffer)) & (metadata['POS'] <= (stop + buffer))]

    exclude_snps = snp_info[snp_info.GenTrain_Score.astype(
        float) <= min_gentrain].snpID.unique()
    snp_info = snp_info[~snp_info.snpID.isin(exclude_snps)][['snpID', 'POS']]

    return snp_info


def create_overlapping_windows(data, window_size, num_intervals):
    """
    Creates overlapping windows spanning the region of interest based on specified window size and number of intervals. 
    It prevents an overlap of zero and ensures that windows do not exceed the bounds of the region of interest.

    Arguments:
    data (list or array-like): Range of base pair positions over which windows are to be created.
    window_size (int): The size of each window.
    num_intervals (int): The number of intervals/windows to create.

    Returns:
    tuple: A tuple containing:
           - start (list): A list of starting points of each window.
           - stop (list): A list of stopping points of each window.
    """
    # Finds necessary overlap to reach the end of the data w/ specified # windows
    total_data_points = len(data)

    # Prevents overlap of 0
    overlap = max(1, int((total_data_points - window_size) /
                  max(1, num_intervals - 1)))
    print(f'{num_intervals} total windows overlapping by: {overlap} base pairs')

    start, stop = [], []
    start_index = 0

    while start_index + window_size <= total_data_points:
        end_index = start_index + window_size
        start.append(data[start_index])

        # Prevent out-of-bounds error
        stop.append(data[min(end_index, total_data_points - 1)])
        start_index += overlap

    return start, stop


def mean_within_interval(row, col_name, df):
    """
    Calculates the mean value of a specified column within a defined interval.

    Arguments:
    row (pandas.Series): A row from a DataFrame containing 'START' and 'STOP' 
                         columns to define the interval.
    col_name (str): The name of the column in the DataFrame for which the mean is
                    calculated.
    df (pandas.DataFrame): The DataFrame containing the data points and 
                           positions.

    Returns:
    float: The mean value of the specified column within the interval. 
           Returns NaN if no data points fall within the interval.
    """
    interval_mask = (df['POS'] >= row['START']) & (
        df['POS'] <= row['STOP'])
    return df.loc[interval_mask, col_name].mean()


def iqr_within_interval(row, col_name, df):
    """
    Calculates the interquartile range (IQR) of a specified column within a defined interval.

    Arguments:
    row (pandas.Series): A row from a DataFrame containing 'START' and 'STOP' 
                         columns to define the interval.
    col_name (str): The name of the column in the DataFrame for which the IQR is
                    calculated.
    df (pandas.DataFrame): The DataFrame containing the data points and positions.

    Returns:
    float: The IQR of the specified column within the interval. 
           Returns NaN if no data points fall within the interval.
    """
    interval_mask = (df['POS'] >= row['START']) & (
        df['POS'] <= row['STOP'])
    return iqr(df.loc[interval_mask, col_name], interpolation='midpoint')


def std_within_interval(row, col_name, df):
    """
    Calculates the standard deviation of a specified column within a defined interval.

    Arguments:
    row (pandas.Series): A row from a DataFrame containing 'START' and 'STOP' 
                         columns to define the interval.
    col_name (str): The name of the column in the DataFrame for which the standard
                    deviation is calculated.
    df (pandas.DataFrame): The DataFrame containing the data points and positions.

    Returns:
    float: The standard deviation of the specified column within the interval. 
           Returns NaN if no data points fall within the interval.
    """
    interval_mask = (df['POS'] >= row['START']) & (
        df['POS'] <= row['STOP'])
    return df.loc[interval_mask, col_name].std()


def dosage_full(row, col_name, df):
    """
    Calculates the dosage value of a specified column within a defined interval
    by summing the values within the interval defined by the 'START' and 'STOP'
    positions in the row and dividing by the total number of CNV candidates in
    the full DataFrame.

    Arguments:
    row (pandas.Series): A row from a DataFrame containing 'START' and 'STOP' 
                         columns to define the interval.
    col_name (str): The name of the column in the DataFrame for which the dosage is
                    calculated.
    df (pandas.DataFrame): The DataFrame containing the data points and positions.

    Returns:
    float: The dosage value of the specified column within the interval. 
           Returns 0 if the DataFrame is empty or no data points fall within the interval.
    """
    interval_mask = (df['POS'] >= row['START']) & (
        df['POS'] <= row['STOP'])
    calls = sum(df.loc[interval_mask, col_name])

    # Catches potential division by zero error
    if len(df) > 0:
        dosage = calls/len(df)
    else:
        dosage = 0
    return dosage


def create_train_set():
    """
    Placeholder for creating a training set.
    Currently can create training set with required columns manually and train
    using run_lstm_model.py and model_methods.py.
    """
    pass


def create_test_set(master_key, num_samples, training_file, snp_metrics_path, out_path, chrom, study_name='all', interval_name=None):
    """
    Creates a test set from a master key file by selecting samples that do not overlap with an existing training set.
    This function can filter samples based on a specified study or cohort name and checks for the 
    existence of SNP metrics/signal intensity files for each sample.

    Arguments:
    master_key (str): Path to the master key file (supports .txt or .csv).
    num_samples (int): The number of samples to include in the test set.
    training_file (str): Path to the training set file to avoid overlapping samples.
    snp_metrics_path (str): Path to the directory containing SNP metrics for the samples.
    out_path (str): The path and file name (without extension) for saving the output test set CSV.
    study_name (str, optional): The name of the study to filter the samples by. Defaults to 'all' for no filtering.
    interval_name (str, optional): The name of the gene/interval to filter training data by to catch overlaps. Defaults to None.

    Returns:
    None: Outputs a CSV file containing the selected test set with the columns 'IID' and 'snp_metrics_path'.
          Prints the number of samples successfully included in the test set.
    """

    # Read the master key file
    if master_key.endswith('.txt'):
        master = pd.read_csv(master_key, sep='\t', low_memory = False)
    elif master_key.endswith('.csv'):
        master = pd.read_csv(master_key)

    # Filter by study name if specified
    if study_name.lower() != 'all':
        full_samples_list = master[master.study == study_name]
    else:
        full_samples_list = master

    # Read and process the training file to avoid overlapping samples
    train_df = pd.read_csv(training_file)
    train_df.columns = map(str.lower, train_df.columns)

    # Checks training file for specific interval
    if "interval" in train_df.columns:
        train_interval = train_df[train_df.interval == interval_name]
        open_ids = full_samples_list[~full_samples_list.IID.isin(
            train_interval.iid)]
    else:
        open_ids = full_samples_list[~full_samples_list.IID.isin(train_df.IID)]

    # Select a subset of non-overlapping samples for the test set
    k = min(len(open_ids), num_samples)
    test_filenames = random.sample(sorted(open_ids.IID), k=k)

    # Create a DataFrame for the test set
    if 'label' in master.columns:
        test_set = master[['IID', 'label']][master['IID'].isin(test_filenames)]
        test_set.reset_index(drop=True, inplace=True)
    else:
        test_set = master.IID[master['IID'].isin(test_filenames)]
        test_set = pd.DataFrame(test_set).reset_index(drop=True)
        test_set.columns = ['IID']

    test_set['snp_metrics_path'] = 0
    remove = []

    # Verify existence of SNP metrics for each sample
    for i in range(len(test_set)):
        sample = test_set.IID.iloc[i]
        code = sample.split('_')[0]

        mfile1 = f'{snp_metrics_path}/{code}/{sample}/chromosome={chrom}'

        if os.path.isdir(mfile1):
            test_set['snp_metrics_path'].iloc[i] = mfile1
        else:
            remove.append(sample)  # remove these from test set

    # Finalize the test set by removing samples without valid SNP metrics
    test_set_final = test_set[~test_set.IID.isin(remove)]
    test_set_final.to_csv(f'{out_path}_testing_IDs.csv', index=False)
    print(f'{len(test_set_final)} of requested {num_samples} samples have necessary SNP metrics')


def make_window_df(start, stop, split_interval, window_count, buffer):
    """
    Creates a DataFrame representing genomic windows within a specified interval.
    This function divides a gene/interval into smaller overlapping windows
    with calculations based on the chosen number of splits and buffer size.

    Arguments:
    start (int): The starting position of the interval of interest.
    stop (int): The stopping position of the interval of interest.
    split_interval (int): The number of divisions for the interval to create 
                          windows.
    window_count (int): The number of overlapping windows to create.
    buffer (int): The buffer size to extend beyond the start and stop positions.

    Returns:
    pandas.DataFrame: A DataFrame containing 'START' and 'STOP' columns, each 
                      representing the beginning and end of a window within 
                      the specified range.
    """
    window_size = round(((stop+buffer) - (start-buffer))/split_interval)
    print(
        f'Interval of interest split {split_interval} times with window size of {window_size} base pairs')

    # Create intervals with no gaps between windows
    l = np.arange((start-buffer), (stop+buffer)+window_size, window_size)
    no_gaps = [value for value in l[1:-1] for _ in range(2)]
    no_gaps.insert(0, l[0])
    no_gaps.insert(len(no_gaps), l[-1])

    # Aggregate SNP metrics into features within each window
    window_df = pd.DataFrame({'START': no_gaps[::2]})
    window_df['STOP'] = window_df['START'] + window_size

    # Create overlapping windows
    data_range = range((start-buffer), (stop+buffer))
    start, stop = create_overlapping_windows(
        data_range, window_size, window_count)

    # Final DataFrame with windows
    window_df = pd.DataFrame({'START': start, 'STOP': stop})

    return window_df


def fill_window_df(sample_data):
    """
    Fills a DataFrame with metrics for genomic windows based on SNP data.
    It identifies CNV candidates based on BAllele Frequency (BAF) and 
    Log R Ratio (LRR) thresholds and calculates statistical metrics for
    these candidates across specified intervals.

    Arguments:
    sample_data (tuple): A tuple containing the following:
        - out_path (str): Path for saving the output CSV.
        - sample (str): Sample ID.
        - snp_metrics_file (str): Path to the file containing SNP metrics data.
        - split_interval (int): The number of divisions for creating windows.
        - total_windows (int): The number of overlapping windows to create.
        - cnv_exists (bool): Whether CNV data exists for this sample.
        - chrom (str): Chromosome identifier.
        - start (int): The start position of the interval.
        - stop (int): The stop position of the interval.
        - buffer (int): Buffer size for extending the start and stop positions.
        - min_gentrain (float): Minimum GenTrain score threshold for filtering SNPs.
        - bim_file (str or None): Path to a BIM file for filtering QCed SNPs.
        - pvar_file (str or None): Path to a PVAR file for filtering QCed SNPs.

    Returns:
    None: Outputs a CSV file containing genomic window metrics with columns such as 
          'START', 'STOP', 'dosage_interval', 'dosage_full', and statistical features 
          like 'avg_baf', 'std_lrr', and 'iqr_baf'. The CSV file also includes an 
          identifier for CNV presence and basic sample metadata.
    """
    out_path, sample, snp_metrics_file, snp_info, split_interval, total_windows, cnv_exists, chrom, start, stop, buffer, min_gentrain, bim_file, pvar_file = sample_data
    window_df = make_window_df(
        start, stop, split_interval, total_windows, buffer)
    metrics_df = pd.read_parquet(snp_metrics_file)

    # Filter data to included SNPs
    sample_interval = metrics_df.merge(snp_info, on='snpID', how='inner')

    # Identify CNV types based on BAF and LRR thresholds
    sample_interval['BAF_insertion'] = np.where((sample_interval['BAF'].between(
        0.65, 0.85, inclusive='neither')) | (sample_interval['BAF'].between(0.15, 0.35, inclusive='neither')), 1, 0)
    sample_interval['L2R_deletion'] = np.where(
        sample_interval['LRR'] < -0.2, 1, 0)
    sample_interval['L2R_duplication'] = np.where(
        sample_interval['LRR'] > 0.2, 1, 0)
    sample_interval['BAF_middle'] = np.where((sample_interval['BAF'] >= 0.15) & (
        sample_interval['BAF'] <= 0.85), 1, 0)

    sample_interval['ALT_pred'] = np.where(sample_interval['BAF_insertion'] == 1, '<INS>',
                                    np.where(sample_interval['L2R_deletion'] == 1, '<DEL>',
                                    np.where(sample_interval['L2R_duplication'] == 1, '<DUP>', '')))
    sample_interval['CNV_call'] = np.where(sample_interval['ALT_pred'] == '', 0,
                                    np.where(sample_interval['ALT_pred'] != '', 1, ''))

    # Extract CNV candidates
    pred_cnv = sample_interval[sample_interval['CNV_call'] == '1']

    # Calculate window-level metrics for ML features
    sample_interval = sample_interval.astype(
        {'BAF': 'float', 'LRR': 'float', 'CNV_call': 'int'})
    pred_cnv = pred_cnv.astype(
        {'BAF': 'float', 'LRR': 'float', 'CNV_call': 'int'})
    window_df['dosage_interval'] = window_df.apply(
        lambda row: mean_within_interval(row, 'CNV_call', sample_interval), axis=1)
    window_df['dosage_full'] = window_df.apply(
        lambda row: dosage_full(row, 'CNV_call', pred_cnv), axis=1)

    window_df['del_dosage_full'] = window_df.apply(
        lambda row: dosage_full(row, 'L2R_deletion', pred_cnv), axis=1)
    window_df['dup_dosage_full'] = window_df.apply(
        lambda row: dosage_full(row, 'L2R_duplication', pred_cnv), axis=1)
    window_df['ins_dosage_full'] = window_df.apply(
        lambda row: dosage_full(row, 'BAF_insertion', pred_cnv), axis=1)

    window_df['del_dosage_interval'] = window_df.apply(
        lambda row: mean_within_interval(row, 'L2R_deletion', pred_cnv), axis=1)
    window_df['dup_dosage_interval'] = window_df.apply(
        lambda row: mean_within_interval(row, 'L2R_duplication', pred_cnv), axis=1)
    window_df['ins_dosage_interval'] = window_df.apply(
        lambda row: mean_within_interval(row, 'BAF_insertion', pred_cnv), axis=1)

    window_df['avg_baf'] = window_df.apply(
        lambda row: mean_within_interval(row, 'BAF', pred_cnv), axis=1)
    window_df['avg_mid_baf'] = window_df.apply(lambda row: mean_within_interval(
        row, 'BAF', sample_interval[sample_interval['BAF_middle'] == 1]), axis=1)
    window_df['avg_lrr'] = window_df.apply(
        lambda row: mean_within_interval(row, 'LRR', pred_cnv), axis=1)

    window_df['std_baf'] = window_df.apply(
        lambda row: std_within_interval(row, 'BAF', pred_cnv), axis=1)
    window_df['std_mid_baf'] = window_df.apply(lambda row: std_within_interval(
        row, 'BAF', sample_interval[sample_interval['BAF_middle'] == 1]), axis=1)
    window_df['std_lrr'] = window_df.apply(
        lambda row: std_within_interval(row, 'LRR', pred_cnv), axis=1)

    window_df['iqr_baf'] = window_df.apply(
        lambda row: iqr_within_interval(row, 'BAF', pred_cnv), axis=1)
    window_df['iqr_mid_baf'] = window_df.apply(lambda row: iqr_within_interval(
        row, 'BAF', sample_interval[sample_interval['BAF_middle'] == 1]), axis=1)
    window_df['iqr_lrr'] = window_df.apply(
        lambda row: iqr_within_interval(row, 'LRR', pred_cnv), axis=1)

    # Additional metadata and final formatting
    window_df['cnv_range_count'] = len(pred_cnv)
    window_df['IID'] = sample
    window_df['CHR'] = chrom

    window_df.fillna(0, inplace=True)
    window_df['window'] = window_df.index
    window_df['CNV_exists'] = cnv_exists

    # window_df.to_csv(f'{out_path}_samples_windows.csv',
    #                  mode='a', header=None, index=False)
    window_df.to_csv(f'{out_path}/{sample}_windows.csv', index=False)


def create_app_ready_file(test_set_id_path, test_set_path, test_result_path, out_path, prob_threshold=0.8):
    """
    Creates an Streamlit-ready CSV file by merging and filtering test set outputs.
    It calculates the maximum IQR of LRR values, sums the number of CNVs in the region,
    and filters the data based on a specified probability threshold for prediction values.

    Arguments:
    test_set_id_path (str): Path to the CSV file containing test set IDs and labels.
    test_set_path (str): Path to the CSV file containing test set data with window metrics.
    test_result_path (str): Path to the CSV file containing prediction results.
    out_path (str): Path and filename (without extension) for the output CSV file.
    prob_threshold (float, optional): Probability threshold for filtering samples. 
                                      Defaults to 0.8.

    Returns:
    pandas.DataFrame: A DataFrame containing the filtered results with samples 
                      that meet the probability threshold criteria.
    """
    test_df = pd.read_csv(test_set_path)

    # Rename column if present
    if 'Unnamed: 0' in test_df.columns:
        test_df.rename(columns={'Unnamed: 0': 'window'}, inplace=True)

    # Calculate absolute IQR for LRR
    test_df['abs_iqr_lrr'] = abs(test_df['iqr_lrr'])
    max_iqr = test_df.groupby('IID')['abs_iqr_lrr'].max()

    test_results = pd.read_csv(test_result_path)
    label_path = pd.read_csv(test_set_id_path)

    # Merge all necessary data
    results = test_results.merge(label_path, on='IID', how='inner')
    full_results = results.merge(
        test_df[['IID', 'cnv_range_count']], on='IID', how='inner').drop_duplicates()
    full_results = full_results.merge(max_iqr, on='IID', how='inner')

    # Filter samples based on the probability threshold
    above_probab = full_results[full_results['Pred Values'] >= prob_threshold]
    above_probab.to_csv(f'{out_path}_app_ready.csv', index=False)

    return above_probab


# FIX
def generate_pred_cnvs(sample_data):
    """
    Generates files for samples with predicted CNVs for future plots in the Streamlit app.

    Arguments:
    sample_data (tuple): A tuple containing the following:
        - sample (str): Sample identifier.
        - metrics (str): Path to the file containing SNP metrics data in Parquet format.
        - chrom (str): Chromosome identifier.
        - start (int): The starting position of the interval.
        - stop (int): The stopping position of the interval.
        - out_path (str): Path for saving the output CSV.
        - buffer (int): Buffer size for extending the start and stop positions.
        - min_gentrain (float): Minimum GenTrain score threshold for filtering SNPs.
        - bim_file (str or None): Path to a BIM file for filtering SNPs.
        - pvar_file (str or None): Path to a PVAR file for filtering SNPs.

    Returns:
    None: Outputs a CSV file with columns such as 'snpID', 'chromosome', 'position', 
          'BAF', 'LRR', and predicted CNV types ('ALT_pred', 'CNV_call').
    """
    sample, metrics, snp_info, chrom, start, stop, out_path, buffer, min_gentrain, bim_file, pvar_file = sample_data
    out_dir = os.path.dirname(os.path.abspath(out_path))
    
    metrics_df = pd.read_parquet(metrics)

    # Filter data to included SNPs
    sample_interval = metrics_df.merge(snp_info, on='snpID', how='inner')

    # Identify CNV types based on BAF and LRR thresholds
    sample_interval['BAF_insertion'] = np.where((sample_interval['BAF'].between(
        0.65, 0.85, inclusive='neither')) | (sample_interval['BAF'].between(0.15, 0.35, inclusive='neither')), 1, 0)
    sample_interval['L2R_deletion'] = np.where(
        sample_interval['LRR'] < -0.2, 1, 0)
    sample_interval['L2R_duplication'] = np.where(
        sample_interval['LRR'] > 0.2, 1, 0)

    sample_interval['ALT_pred'] = np.where(sample_interval['BAF_insertion'] == 1, '<INS>',
                                           np.where(sample_interval['L2R_deletion'] == 1, '<DEL>',
                                                    np.where(sample_interval['L2R_duplication'] == 1, '<DUP>', '')))
    sample_interval['CNV_call'] = np.where(sample_interval['ALT_pred'] == '', 0,
                                           np.where(sample_interval['ALT_pred'] != '', 1, ''))

    # Save the results to a CSV file in the 'pred_cnvs' directory
    pred_path = f'{out_dir}/pred_cnvs'
    os.makedirs(pred_path, exist_ok=True)
    sample_interval.to_csv(
        f'{pred_path}/{sample}_full_interval.csv', index=False)


def plot_variants(df, x_col='BAF', y_col='LRR', gtype_col='GT', title='snp plot', opacity=1, midline=False, cnvs=None, xmin=None, xmax=None):
    """
    Plots an interactive scatter plot of genetic variants with customizable axes, colors, and features.

    Arguments:
    df (pandas.DataFrame): The DataFrame containing the data to be plotted.
    x_col (str, optional): The column name for the x-axis. Defaults to 'BAF'.
    y_col (str, optional): The column name for the y-axis. Defaults to 'LRR'.
    gtype_col (str, optional): The column name used for coloring the points. Defaults to 'GT'.
    title (str, optional): The title of the plot. Defaults to 'snp plot'.
    opacity (float, optional): The opacity level of the points (0 to 1). Defaults to 1.
    midline (bool, optional): Whether to add a midline representing the average trend. Defaults to False.
    cnvs (pandas.DataFrame, optional): A DataFrame with CNV data to overlay on the plot. Defaults to None.
    xmin (float, optional): The minimum x-axis limit. Defaults to None.
    xmax (float, optional): The maximum x-axis limit. Defaults to None.

    Returns:
    plotly.graph_objs._figure.Figure: The Plotly figure object representing the scatter plot.
    """
    d3 = px.colors.qualitative.D3

    cmap = {
        'AA': d3[0],
        'AB': d3[1],
        'BA': d3[1],
        'BB': d3[2],
        'NC': d3[3]
    }

    cmap_ALT = {
        '<INS>': d3[0],
        '<DEL>': d3[1],
        '<DUP>': d3[2],
        '<None>': d3[7]
    }

    # Set default x-axis limits if not provided
    if not xmin and not xmin:
        xmin, xmax = df[x_col].min(), df[x_col].max()

    ymin, ymax = df[y_col].min(), df[y_col].max()
    xlim = [xmin-.1, xmax+.1]
    ylim = [ymin-.1, ymax+.1]

    lmap = {'BAF': 'BAF', 'LRR': 'LRR'}
    smap = {'Control': 'circle', 'PD': 'diamond-open-dot'}

    # Choose color map based on genotype column
    if gtype_col == 'ALT_pred' or gtype_col == 'ALT':
        cmap_choice = cmap_ALT
    else:
        cmap_choice = cmap

    if isinstance(cnvs, pd.DataFrame):
        fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map=cmap_choice,
                         color_continuous_scale=px.colors.sequential.matter, width=650, height=497, labels=lmap, symbol_map=smap, hover_data=[gtype_col])
        fig.update_traces(opacity=opacity)
        fig.add_traces(px.scatter(cnvs, x=x_col, y=y_col, hover_data=[
                       gtype_col]).update_traces(marker_color="black").data)
    else:
        if gtype_col == None:
            fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_sequence=[
                             'grey'], width=650, height=497, labels=lmap, symbol_map=smap, opacity=opacity)
        else:
            fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map=cmap_choice,
                             width=650, height=497, labels=lmap, symbol_map=smap, opacity=opacity)

    if midline:
        # Calculate average y-value for each unique x-value
        unique_x = np.linspace(min(df[x_col]), max(df[x_col]), num=50)

        # Creates bins and calculate average y-value within each bin
        df['x_bin'] = pd.cut(df[x_col], bins=unique_x)
        grouped_df = df[[x_col, 'x_bin', y_col]].groupby(
            'x_bin').mean().reset_index()

        # Plot a midline
        fig.add_traces(px.line(grouped_df, x=x_col, y=y_col).update_traces(
            line=dict(color='red', width=3), name='Average Line').data)

    fig.update_xaxes(range=xlim, nticks=10, zeroline=False)
    fig.update_yaxes(range=ylim, nticks=10, zeroline=False)
    fig.update_layout(margin=dict(r=76, t=63, b=75))
    fig.update_layout(legend_title_text='CNV Range Class')
    fig.update_layout(title_text=f'<b>{title}<b>')

    return fig
