import os
import sys
import subprocess
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import seaborn as sns
import random
from datetime import datetime
from PIL import Image
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode
from io import StringIO, BytesIO

from hold_method import plot_variants


st.set_page_config(
    page_title="CNV Prediction Evaluator",
    layout="wide",
    initial_sidebar_state="expanded"
)

# add to hold_data?
def generate_sample(cohort_samples):
    st.session_state['no_plot'] = False
    sample_options = cohort_samples[~cohort_samples.isin(st.session_state[samples_seen])]

    if len(sample_options) == 0:
        st.error("No more samples to run through!")
        st.session_state['no_plot'] = True
    else:
        st.session_state['sample_name'] = random.choice(list(sample_options)) # is list function necessary here
    # return no_plot

def plot_sample():
    sample_df_interval = pd.read_csv(f"CNV_app/data/{st.session_state['model_choice']}/{st.session_state['cohort_choice']}/{st.session_state['gene_choice']}/app/pred_cnvs/{st.session_state['sample_name']}_full_interval.csv")
    sample_df_interval['ALT_pred'].fillna('<None>', inplace = True)
    pred_cnv = sample_df_interval[sample_df_interval['CNV_call'] == 1]

    # st.dataframe(sample_df_interval)
    # st.dataframe(pred_cnv)

    # predictions only
    fig_lrr = plot_variants(pred_cnv[pred_cnv['ALT_pred']!='<INS>'], x_col='position', y_col='LogRRatio', gtype_col = 'ALT_pred', title=f'{st.session_state["gene_choice"]} Interval CNV Predictions Only')
    xmin, xmax = pred_cnv['position'].min(), pred_cnv['position'].max()
    st.plotly_chart(fig_lrr)

    fig_baf = plot_variants(pred_cnv[pred_cnv['ALT_pred']=='<INS>'], x_col='position', y_col='BAlleleFreq', gtype_col = 'ALT_pred', title=f'{st.session_state["gene_choice"]} Interval CNV Predictions Only', xmin=xmin, xmax=xmax)
    st.plotly_chart(fig_baf)

    # just CNV predicted variants among all variants
    fig_lrr_full = plot_variants(sample_df_interval[sample_df_interval['ALT_pred']!='<INS>'], x_col='position', gtype_col = 'ALT_pred', y_col='LogRRatio', title=f'{st.session_state["gene_choice"]} Interval Colored by CNV Type')
    st.plotly_chart(fig_lrr_full)
    
    fig_baf_full = plot_variants(sample_df_interval[(sample_df_interval['ALT_pred']!='<DEL>') & (sample_df_interval['ALT_pred']!='<DUP>')], x_col='position', gtype_col = 'ALT_pred', y_col='BAlleleFreq', title=f'{st.session_state["gene_choice"]} Interval Colored by CNV Type')
    st.plotly_chart(fig_baf_full)

    # B&W all samples
    bw_lrr_full = plot_variants(sample_df_interval, x_col='position', gtype_col = None, y_col='LogRRatio',midline=True, title=f'All Variants in {st.session_state["gene_choice"]} Interval with Average Line', opacity=0.3)
    st.plotly_chart(bw_lrr_full)
    
    bw_baf_full = plot_variants(sample_df_interval, x_col='position', gtype_col = None, y_col='BAlleleFreq', title=f'All Variants in {st.session_state["gene_choice"]} Interval', opacity=0.3)
    st.plotly_chart(bw_baf_full)

### Create sidebar options
# will expand on model selection - need to distinguish test results based on model version
st.sidebar.markdown('### Choose a Model')
models = ['Preliminary Deletion Model', 'Preliminary Duplication Model', 'Updated Deletion Model', 'Updated Duplication Model']
model_name = st.sidebar.selectbox(label = 'Model Selection', label_visibility = 'collapsed', options=models)
models_dict = {'Preliminary Deletion Model': 'del_prelim_model', 'Preliminary Duplication Model': 'dup_prelim_model',
                'Updated Deletion Model': 'del_updated_model', 'Updated Duplication Model': 'dup_updated_model'}

# split GP2 by cohort but also provide full release option
st.sidebar.markdown('### Choose a GP2 Cohort')
cohorts = next(os.walk(f'CNV_app/data/{models_dict[model_name]}'))[1]
cohort_name = st.sidebar.selectbox(label = 'Cohort Selection', label_visibility = 'collapsed', options=cohorts)

# will expand on this selection to all NDDs
st.sidebar.markdown('### Choose an NDD-Related Gene')

# Mix of disease-related genes
# Main focus for benchmarking
genes = ['PARK2', 'LINGO2', 'MAPT', 'SNCA']

# Complete list of explored genes
# genes = ['PARK2', 'APP', 'LINGO2', 'SNCA', 'CHRNA7', 'CYFIP1', 'CPNE4', 'TPCN1', 'PSEN1', 'CREB1', 'LRRK2', 'GBA', 'MAPT', 'ABI3', 'ABCA7', '22q_small', '22q',
#         'PARK7', 'ATP13A2', 'PINK1', 'EIF4G1', 'FBXO7', 'VPS35', 'PLA2G6', 'GIGYF2','DCTN1', 'DNAJC6', 'GCH1', 'GRN']

# PD genes of interest
# genes = ['PARK2', 'SNCA', 'LRRK2', 'GBA', 'MAPT','PARK7', 'ATP13A2', 'PINK1', 'EIF4G1', 'FBXO7', 
#           'VPS35', 'PLA2G6', 'GIGYF2','DCTN1', 'DNAJC6', 'GCH1', 'GRN']
gene_name = st.sidebar.selectbox(label = 'NDD-Related Gene Selection', label_visibility = 'collapsed', options=genes)

# can change into a function
option_change = False
if 'gene_choice' in st.session_state:
    if gene_name != st.session_state['gene_choice']:
        option_change = True
if 'model_choice' in st.session_state:
    if model_name != st.session_state['model_choice']:
        option_change = True
if 'cohort_choice' in st.session_state:
    if cohort_name != st.session_state['cohort_choice']:
        option_change = True

st.session_state['cohort_choice'] = cohort_name
st.session_state['model_choice'] = models_dict[model_name]
st.session_state['gene_choice'] = gene_name

### Main Page
# App development notes:
# main page: randomly generate sample plots - make download button for text file of sample names in each category
# second page: look through plots of samples in each category - will need to keep csv saved of samples in each so can reload; add login?
# second page: upload text file and check plots

st.title('Evaluation of CNV Predictions')
model_path = f'CNV_app/data/{models_dict[model_name]}/{cohort_name}/{gene_name}/app/GP2_{cohort_name}_{gene_name}_app_ready.csv'

if not os.path.isfile(model_path):
    st.error('No CNVs to display!')
else:
    model_results = pd.read_csv(model_path)

    with st.expander("Filter Displayed Samples"):
        # if st.session_state["model_choice"] == 'dup_prelim_model':
        #     probab_options = [0.6, 0.7, 0.8, 0.9, 1]
        # else:
        #     probab_options = [0.8, 0.9, 1]
        probab_options = [0.6, 0.7, 0.8, 0.9, 1]
        confidence = st.select_slider('Display samples with prediction probability of at least:', options=probab_options, value = 1)

        # adjust these thresholds 
        iqr_range = np.linspace(min(abs(model_results['abs_iqr_lrr'])), max(abs(model_results['abs_iqr_lrr'])), num=50)
        iqr_threshold = st.select_slider('Maximum Absolute Value LRR range threshold:', options=iqr_range, value = max(iqr_range))

        # potentially use standard deviations here
        min_cnv_count = min(model_results['cnv_range_count'])
        max_cnv_count = max(model_results['cnv_range_count'])
        lower_range = np.linspace(min_cnv_count, max_cnv_count - (0.5 * max_cnv_count), num=50, dtype=int)
        upper_range = np.linspace(max_cnv_count- (0.5 * max_cnv_count), max_cnv_count, num=50, dtype=int)
        selected_value_lower = min(lower_range)
        selected_value_higher = max(upper_range)

        # lower_range = range(40, 110, 10)
        # upper_range = range(100, 500, 50)
        # selected_value_lower = 80
        # selected_value_higher = 300

        lower_range_threshold = st.select_slider('Minimum threshold for CNV count:', options=lower_range, value = selected_value_lower)
        upper_range_threshold = st.select_slider('Maximum threshold for CNV count :', options=upper_range, value = selected_value_higher)

        exp1, exp2, exp3 = st.columns([1.5,1,1])
        st.session_state['threshold_submit']  = exp2.checkbox('Plot With Thresholds')

    if st.session_state['threshold_submit'] :
        # may need to add 'submit' button bc of delay with generating IIDs with these metrics
        threshold_results = model_results[(model_results['abs_iqr_lrr'] <= iqr_threshold) & (model_results['cnv_range_count'] >= lower_range_threshold) & 
                                        (model_results['cnv_range_count'] <= upper_range_threshold) & (model_results['Pred Values'] >= confidence)]
        cohort_samples = threshold_results['IID']
    else:
        # cohort_samples = model_results.IID[model_results['Pred Values'] >= confidence]
        cohort_samples = model_results.IID[model_results['Pred Values'] == 1]

    samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'
    if samples_seen not in st.session_state:
        st.session_state[samples_seen] = []
    if 'yes_choices' not in st.session_state:
        st.session_state['yes_choices'] = []
    if 'maybe_choices' not in st.session_state:
        st.session_state['maybe_choices'] = []
    if 'no_choices' not in st.session_state:
        st.session_state['no_choices'] = []
    if 'yes_gene' not in st.session_state:
        st.session_state['yes_gene'] = []
    if 'maybe_gene' not in st.session_state:
        st.session_state['maybe_gene'] = []
    if 'no_gene' not in st.session_state:
        st.session_state['no_gene'] = []
    if 'yes_type' not in st.session_state:
        st.session_state['yes_type'] = []
    if 'maybe_type' not in st.session_state:
        st.session_state['maybe_type'] = []
    if 'no_type' not in st.session_state:
        st.session_state['no_type'] = []
    if 'sample_name' not in st.session_state:
        generate_sample(cohort_samples)
    if option_change:
        generate_sample(cohort_samples)

    col1, col2 = st.columns([3,0.7])
    btn1, btn2, btn3, btn4 = st.columns([0.5,0.5,0.5,0.5])

    # Display any initial warnings
    if 'Artifact Warning' in model_results.columns and model_results['Artifact Warning'].iloc[0] == 1:
        col1.error('Please note that a high number of samples in this gene were predicted to have CNVs, which may indicate that an artifact or other array-based issue is displayed.')

    col1.markdown(f'##### _Would you consider Sample {st.session_state["sample_name"]} a structural variant?_')
    col2.markdown(f'Prediction probability of {str(round(model_results.loc[model_results.IID == st.session_state["sample_name"], "Pred Values"].iloc[0], 2))}')

    if len(st.session_state[samples_seen]) == 0:
        st.session_state[samples_seen].append(st.session_state['sample_name'])

    ## create on_change function for checkbox
    if len(st.session_state[samples_seen]) == 1 and st.session_state['threshold_submit']:
        st.session_state[samples_seen].append(st.session_state['sample_name'])

    if not st.session_state['no_plot']:
        yes = btn1.button('Yes', use_container_width = True)
        maybe = btn2.button('Maybe', use_container_width = True)
        no_btn = btn3.button('No', use_container_width = True)

        # eventually make it so that this button prompts multiple choice/text entry for custom annotation
        other_cnv = btn4.button('Other CNV', use_container_width = True)
        plot_sample()
    else:
        yes = btn1.button('Yes', disabled = True, use_container_width = True)
        maybe = btn2.button('Maybe', disabled = True, use_container_width = True)
        no_btn = btn3.button('No', disabled = True, use_container_width = True)
        other_cnv = btn4.button('Other CNV', disabled = True, use_container_width = True)

    if yes:
        samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'
        st.session_state[samples_seen].append(st.session_state['sample_name'])
        st.session_state['yes_choices'].append(st.session_state[samples_seen][-2])
        st.session_state['yes_gene'].append(st.session_state["gene_choice"])
        st.session_state['yes_type'].append(st.session_state["model_choice"].split('_')[0])
    elif maybe:
        samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'
        st.session_state[samples_seen].append(st.session_state['sample_name'])
        st.session_state['maybe_choices'].append(st.session_state[samples_seen][-2])
        st.session_state['maybe_gene'].append(st.session_state["gene_choice"])
        st.session_state['maybe_type'].append(st.session_state["model_choice"].split('_')[0])
    elif no_btn:
        samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'
        st.session_state[samples_seen].append(st.session_state['sample_name'])
        st.session_state['no_choices'].append(st.session_state[samples_seen][-2])
        st.session_state['no_gene'].append(st.session_state["gene_choice"])
        st.session_state['no_type'].append(st.session_state["model_choice"].split('_')[0])
    elif other_cnv:
        samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'
        st.session_state[samples_seen].append(st.session_state['sample_name'])
        st.session_state['yes_choices'].append(st.session_state[samples_seen][-2])
        st.session_state['yes_gene'].append(st.session_state["gene_choice"])
        
        cnv_type = st.session_state["model_choice"].split('_')[0]
        if cnv_type == 'del':
            st.session_state['yes_type'].append('dup')
        elif cnv_type == 'dup':
            st.session_state['yes_type'].append('del')

    side_btn1, side_btn2, side_btn3 = st.sidebar.columns([0.5, 1, 0.5])
    
    yes_report = pd.DataFrame({'Yes Samples': st.session_state['yes_choices'], 'Interval': st.session_state['yes_gene'], 'Type': st.session_state['yes_type']})
    maybe_report = pd.DataFrame({'Maybe Samples': st.session_state['maybe_choices'], 'Interval': st.session_state['maybe_gene'], 'Type': st.session_state['maybe_type']})
    no_report = pd.DataFrame({'No Samples': st.session_state['no_choices'], 'Interval': st.session_state['no_gene'], 'Type': st.session_state['no_type']})

    # Add download button and make choices into dataframes/dictionaries that include gene name where CNV was found/not found
    with st.sidebar.expander("View Reported Samples"):
        st.data_editor(yes_report,
                        hide_index=True,
                        use_container_width=True
                    )
        st.data_editor(maybe_report,
                        hide_index=True,
                        use_container_width=True
                    )
        st.data_editor(no_report,
                        hide_index=True,
                        use_container_width=True
                    )

    # data_editor in streamlit 1.32 version has download to csv option built-in
    save = side_btn2.button('Save Report')
    if save:
        # time stamps for file naming to prevent overwriting
        current_date = str(datetime.now().strftime("%Y-%m-%d"))
        current_time = str(datetime.now().strftime("%H-%M-%S"))

        yes_report.to_csv(f'CNV_app/data/yes_samples_{current_date}_{current_time}.csv', index = False)
        maybe_report.to_csv(f'CNV_app/data/maybe_samples_{current_date}_{current_time}.csv', index = False)
        no_report.to_csv(f'CNV_app/data/no_samples_{current_date}_{current_time}.csv', index = False)
