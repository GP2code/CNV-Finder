import pandas as pd
import os
import shutil
import numpy as np
import glob
import statsmodels.api as sm
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
import importlib
import plotly
import random
from numpy.core.numeric import NaN
from sklearn.preprocessing import MinMaxScaler

# Supress copy warning.
pd.options.mode.chained_assignment = None

# create chromosome interval here
def check_interval(interval_name, interval_file = 'ref_files/glist_hg38_intervals.csv'):
    # can catch for other common NDD genes
    if interval_name == 'PRKN':
        interval_name = 'PARK2'
        
    interval_df = pd.read_csv(interval_file)
    positions = interval_df[interval_df.NAME == interval_name]
    return positions

def find_interval(df, position):
    for start, stop in zip(df['POS'], df['POS_stop']):
        if start <= position <= stop:
            return f"{start}-{stop}"
    return None

def mean_within_interval(row, col_name, df):
    interval_mask = (df['position'] >= row['START']) & (df['position'] <= row['STOP'])
    return df.loc[interval_mask, col_name].mean()

def create_overlapping_windows(data, window_size, num_intervals):
    # Finds necessary overlap to reach the end of the data w/ specified # windows
    total_data_points = len(data)
    # Prevents overlap of 0
    overlap = max(1, int((total_data_points - window_size) / max(1, num_intervals - 1)))
    print(f'{num_intervals} total windows overlapping by: {overlap} base pairs')

    start, stop = [], []
    start_index = 0

    while start_index + window_size <= total_data_points:
        end_index = start_index + window_size
        start.append(data[start_index])

        # Prevent out of bounds error
        stop.append(data[min(end_index, total_data_points - 1)])
        start_index += overlap

    return start, stop

def iqr_within_interval(row, col_name, df):
    interval_mask = (df['position'] >= row['START']) & (df['position'] <= row['STOP'])
    return stats.iqr(df.loc[interval_mask, col_name], interpolation = 'midpoint')

def std_within_interval(row, col_name, df):
    interval_mask = (df['position'] >= row['START']) & (df['position'] <= row['STOP'])
    return df.loc[interval_mask, col_name].stdev()

def dosage_within_gene(row, col_name, df):
    interval_mask = (df['position'] >= row['START']) & (df['position'] <= row['STOP'])
    calls = sum(df.loc[interval_mask, col_name])

    # catch potential division by zero error
    if len(df) > 0:
        dosage = calls/len(df)
    else:
        dosage = 0
    return dosage

def create_train_set():
    # create file if doesn't exist, or append to existing file so you can add pos/neg cases
    # can potentially supply path to folder with all hand-checked samples
    # will accept output from app
    pass

def create_test_set(master_key, num_samples, training_file, snp_metrics_path, out_path, study_name = 'all'):
    # need master key, study name or GP2_all for whole GP2, need path to training set so no overlap with test set 
    # can later change so that it only doesn't repeat IIDs for the specific gene of interest
    if master_key.endswith('.txt'):
        master = pd.read_csv(master_key, sep='\t')
    elif master_key.endswith('.csv'):
        master = pd.read_csv(master_key)
        
    full_samples_list = master[master.GDPR == 0]

    if study_name.lower() != 'all':
        full_samples_list = full_samples_list[full_samples_list.study == study_name]

    train_df = pd.read_csv(training_file)
    
    open_ids = full_samples_list[~full_samples_list.IID.isin(train_df.IID)]
    k = min(len(open_ids), num_samples)
    test_filenames = random.sample(set(open_ids.IID), k=k)

    test_set = master[['IID', 'label']][master['IID'].isin(test_filenames)]
    test_set.reset_index(drop = True, inplace = True)
    test_set['snp_metrics_path'] = 0
    remove = []

    for i in range(len(test_set)):
        sample = test_set.IID.iloc[i]
        label = test_set.label.iloc[i]
        code = sample.split('_')[0]
        
        mfile1 = f'{snp_metrics_path}/{code}/snp_metrics_{code}/Sample_ID={sample}'
    
        if os.path.isdir(mfile1):
            test_set['snp_metrics_path'].iloc[i] = mfile1
        else:
            remove.append(sample) # remove these from test set 

    # later fix to account for missing snp metrics that result in < num_samples
    test_set_final = test_set[~test_set.IID.isin(remove)]
    test_set_final.to_csv(f'{out_path}_testing_IDs.csv', index = False)
    print(f'{len(test_set_final)} of requested {num_samples} samples have necessary SNP metrics')

def make_window_df(chr, start, stop, split_interval, window_count, buffer):
    window_size = round(((stop+buffer) - (start-buffer))/split_interval)
    print(f'Interval of interest split {split_interval} times with window size of {window_size} base pairs')
    
    # Create intervals with no gaps between windows
    l = np.arange((start-buffer), (stop+buffer)+window_size, window_size)
    no_gaps = [value for value in l[1:-1] for _ in range(2)]
    no_gaps.insert(0, l[0])
    no_gaps.insert(len(no_gaps), l[-1])
    
    # Will aggregate SNP metrics into features within each window
    window_df = pd.DataFrame({'START':no_gaps[::2]})
    window_df['STOP'] = window_df['START'] + window_size

    # Create overlapping windows
    data_range = range((start-buffer), (stop+buffer))
    start, stop = create_overlapping_windows(data_range, window_size, window_count)
    
    # Final dataframe with windows
    window_df = pd.DataFrame({'START':start, 'STOP': stop})

    # Dataframe that will hold all samples' windows_df
    all_samples = pd.DataFrame(columns=['START', 'STOP', 'dosage_interval', 'dosage_gene', 'std_baf', 'std_lrr', 'iqr_baf', 'iqr_lrr', 'CHR', 'IID'])

    return window_df, all_samples

def fill_window_df(data_file, chr, start, stop, out_path, split_interval = 5, total_window_count = 31, buffer=250000, test_set = True, min_gentrain=0.2, bim_file = None, pvar_file = None):
    samples_list = pd.read_csv(data_file)

    # get clean dataframe to fill
    window_df, all_samples = make_window_df(chr, start, stop, split_interval, total_window_count, buffer)
    
    for i in range(len(samples_list)):
        sample = samples_list.IID.iloc[i]
        code = sample.split('_')[0]
        snp_metrics_file = samples_list.snp_metrics_path.iloc[i]
        metrics_df = pd.read_parquet(snp_metrics_file)

        # may need to run one ancestry label at a time depending on how bim is organized 
        if bim_file or pvar_file:
            if os.path.isfile(bim_file):
                bim = pd.read_csv(bim_file, sep='\s+', header=None, names=['chr','id','pos','bp','a1','a2'], usecols=['id'])
                sample_df = metrics_df.loc[(metrics_df.snpID.isin(bim.id)) & (metrics_df.GenTrain_Score>=min_gentrain)]
            elif os.path.isfile(pvar_file):
                pvar = pd.read_csv(pvar_file, sep='\s+', header=None, names=['#CHROM','POS','ID','REF','ALT'], usecols=['ID'])
                sample_df = metrics_df.loc[(metrics_df.snpID.isin(pvar.ID)) & (metrics_df.GenTrain_Score>=min_gentrain)]
        else:
            sample_df = metrics_df.loc[(metrics_df.GenTrain_Score>=min_gentrain)]
    
        sample_df_interval = sample_df[['snpID', 'chromosome', 'position', 'BAlleleFreq', 'LogRRatio']][(sample_df['chromosome'] == chr) 
                        & (sample_df['position'] >= (start-buffer)) 
                        & (sample_df['position'] <= (stop+buffer))]
        
        # find predicted CNV type
        sample_df_interval['BAF_insertion'] = np.where((sample_df_interval['BAlleleFreq'].between(0.65, 0.85, inclusive='neither')) | (sample_df_interval['BAlleleFreq'].between(0.15, 0.35, inclusive='neither')), 1, 0)
        sample_df_interval['L2R_deletion'] = np.where(sample_df_interval['LogRRatio'] < -0.2, 1, 0)
        sample_df_interval['L2R_duplication'] = np.where(sample_df_interval['LogRRatio'] > 0.2, 1, 0)
        
        sample_df_interval['ALT_pred'] = np.where(sample_df_interval['BAF_insertion'] == 1, '<INS>', 
                                        np.where(sample_df_interval['L2R_deletion'] == 1, '<DEL>', 
                                        np.where(sample_df_interval['L2R_duplication'] == 1, '<DUP>', '')))
        sample_df_interval['CNV_call'] = np.where(sample_df_interval['ALT_pred'] == '', 0, 
                                        np.where(sample_df_interval['ALT_pred'] != '', 1, ''))

        # only where variants fall into CNV ranges
        pred_cnv = sample_df_interval[sample_df_interval['CNV_call'] == '1']

        # gather features for ML model
        sample_df_interval = sample_df_interval.astype({'BAlleleFreq':'float', 'LogRRatio':'float', 'CNV_call':'int'})
        pred_cnv = pred_cnv.astype({'BAlleleFreq':'float', 'LogRRatio':'float', 'CNV_call':'int'})
        window_df['dosage_interval'] = window_df.apply(lambda row: mean_within_interval(row, 'CNV_call', sample_df_interval), axis=1)
        window_df['dosage_gene'] = window_df.apply(lambda row: dosage_within_gene(row, 'CNV_call', pred_cnv), axis=1)
    
        # save for each type when working with all CNV types
        window_df['del_dosage'] = window_df.apply(lambda row: dosage_within_gene(row, 'L2R_deletion', pred_cnv), axis=1)
    
        window_df['avg_baf'] = window_df.apply(lambda row: mean_within_interval(row, 'BAlleleFreq', pred_cnv), axis=1)
        window_df['avg_lrr'] = window_df.apply(lambda row: mean_within_interval(row, 'LogRRatio', pred_cnv), axis=1)
        
        window_df['std_baf'] = window_df.apply(lambda row: mean_within_interval(row, 'BAlleleFreq', pred_cnv), axis=1)
        window_df['std_lrr'] = window_df.apply(lambda row: mean_within_interval(row, 'LogRRatio', pred_cnv), axis=1)
            
        window_df['iqr_baf'] = window_df.apply(lambda row: mean_within_interval(row, 'BAlleleFreq', pred_cnv), axis=1)
        window_df['iqr_lrr'] = window_df.apply(lambda row: mean_within_interval(row, 'LogRRatio', pred_cnv), axis=1)
        
        window_df['cnv_range_count'] = len(pred_cnv)
        window_df['IID'] = sample
        window_df['CHR'] = chr
    
        if not test_set:
            window_df['CNV_exists'] = 0  # figure out with create train set
    
        # restarts index with every new sample - will want index = True when saving to csv for window count
        window_df.fillna(0, inplace = True)
        all_samples = pd.concat([all_samples, window_df])

        # figure out a way to parallelize?

    all_samples.to_csv(f'{out_path}_samples_windows.csv')

def create_app_ready_file(test_set_id_path, test_set_path, test_result_path, out_path, prob_threshold = 0.8):
    test_df = pd.read_csv(test_set_path)
    test_df.rename(columns = {'Unnamed: 0': 'window'}, inplace = True)
    test_df['abs_iqr_lrr'] = abs(test_df['iqr_lrr'])
    max_iqr = test_df.groupby('IID')['abs_iqr_lrr'].max()

    test_results = pd.read_csv(test_result_path)
    label_path = pd.read_csv(test_set_id_path)

    # merge all necessary files
    results = test_results.merge(label_path, on = 'IID', how = 'inner')
    full_results = results.merge(test_df[['IID', 'cnv_range_count']], on = 'IID', how = 'inner').drop_duplicates()
    full_results = full_results.merge(max_iqr, on = 'IID', how = 'inner')

    # will only save samples above probability threshold
    above_probab = full_results[full_results['Pred Values'] >= prob_threshold]
    above_probab.to_csv(f'{out_path}_app_ready.csv', index = False)

    return above_probab

# def generate_pred_cnvs(above_probab_path, chr, start, stop, out_path, buffer = 250000, bim_file = None, pvar_file = None):
def generate_pred_cnvs(metrics, chr, start, stop, out_path, buffer = 250000, min_gentrain= 0.2, bim_file = None, pvar_file = None):
    out_dir = os.path.dirname(os.path.abspath(out_path))
    # samples = pd.read_csv(above_probab_path)

    # if parrallelize, remove for loop
    # for metrics in samples.snp_metrics_path:
    sample = metrics.split('/')[-1].split('=')[-1]
        
    metrics_df = pd.read_parquet(metrics)

    if bim_file or pvar_file:
        if os.path.isfile(bim_file):
            bim = pd.read_csv(bim_file, sep='\s+', header=None, names=['chr','id','pos','bp','a1','a2'], usecols=['id'])
            sample_df = metrics_df.loc[(metrics_df.snpID.isin(bim.id)) & (metrics_df.GenTrain_Score>=min_gentrain)]
        elif os.path.isfile(pvar_file):
            pvar = pd.read_csv(pvar_file, sep='\s+', header=None, names=['#CHROM','POS','ID','REF','ALT'], usecols=['ID'])
            sample_df = metrics_df.loc[(metrics_df.snpID.isin(pvar.ID)) & (metrics_df.GenTrain_Score>=min_gentrain)]
    else:
        sample_df = metrics_df.loc[(metrics_df.GenTrain_Score>=min_gentrain)]
    
    sample_df_interval = sample_df[['snpID', 'chromosome', 'position', 'BAlleleFreq', 'LogRRatio']][(sample_df['chromosome'] == chr) 
                        & (sample_df['position'] >= (start-buffer)) 
                        & (sample_df['position'] <= (stop+buffer))]
    
    # find predicted CNV type
    sample_df_interval['BAF_insertion'] = np.where((sample_df_interval['BAlleleFreq'].between(0.65, 0.85, inclusive='neither')) | (sample_df_interval['BAlleleFreq'].between(0.15, 0.35, inclusive='neither')), 1, 0)
    sample_df_interval['L2R_deletion'] = np.where(sample_df_interval['LogRRatio'] < -0.2, 1, 0)
    sample_df_interval['L2R_duplication'] = np.where(sample_df_interval['LogRRatio'] > 0.2, 1, 0)
    
    sample_df_interval['ALT_pred'] = np.where(sample_df_interval['BAF_insertion'] == 1, '<INS>', 
                                    np.where(sample_df_interval['L2R_deletion'] == 1, '<DEL>', 
                                    np.where(sample_df_interval['L2R_duplication'] == 1, '<DUP>', '')))
    sample_df_interval['CNV_call'] = np.where(sample_df_interval['ALT_pred'] == '', 0, 
                                    np.where(sample_df_interval['ALT_pred'] != '', 1, ''))

    # more simplistic path relies on organized out_path selection
    # pred_path = f'{out_dir}/pred_cnvs/{chr}_{start-buffer}_{stop+buffer}'
    pred_path = f'{out_dir}/pred_cnvs'
    os.makedirs(pred_path, exist_ok=True)
    sample_df_interval.to_csv(f'{pred_path}/{sample}_full_interval.csv', index = False)
        
        
def plot_variants(df, x_col='BAlleleFreq', y_col='LogRRatio', gtype_col='GT', title='snp plot', opacity = 1, cnvs = None, xmin= None, xmax = None):
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

    # gtypes_list = (df[gtype_col].unique())
    if not xmin and not xmin:
        xmin, xmax = df[x_col].min(), df[x_col].max()
    
    ymin, ymax = df[y_col].min(), df[y_col].max()
    xlim = [xmin-.1, xmax+.1]
    ylim = [ymin-.1, ymax+.1]

    lmap = {'BAlleleFreq':'BAF','LogRRatio':'LRR'}
    smap = {'Control':'circle','PD':'diamond-open-dot'}

    if gtype_col == 'ALT_pred' or gtype_col == 'ALT':
        cmap_choice = cmap_ALT
    else:
        cmap_choice = cmap 

    if isinstance(cnvs, pd.DataFrame):
        fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map = cmap_choice, color_continuous_scale=px.colors.sequential.matter, width=650, height=497, labels=lmap, symbol_map=smap, hover_data=[gtype_col])
        fig.update_traces(opacity = opacity)
        fig.add_traces(px.scatter(cnvs, x=x_col, y=y_col, hover_data=[gtype_col]).update_traces(marker_color="black").data)
    else:
        fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map=cmap_choice, width=650, height=497, labels=lmap, symbol_map=smap)

    fig.update_xaxes(range=xlim, nticks=10, zeroline=False)
    fig.update_yaxes(range=ylim, nticks=10, zeroline=False)
    
    fig.update_layout(margin=dict(r=76, t=63, b=75))

    # fig.update_layout(legend=dict(orientation="h", yanchor="bottom", y=1, xanchor="right", x=1))

    fig.update_layout(legend_title_text='CNV Range Class')

    out_dict = {
        'fig': fig,
        'xlim': xlim,
        'ylim': ylim
    }
    
    fig.update_layout(title_text=f'<b>{title}<b>')
    
    return fig

