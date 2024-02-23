import os
import sys
import subprocess
import datetime
import numpy as np
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import seaborn as sns
import random
from PIL import Image
from st_aggrid import GridOptionsBuilder, AgGrid, GridUpdateMode, DataReturnMode
from io import StringIO, BytesIO

from hold_method import plot_clusters


st.set_page_config(
    page_title="CNV Prediction Evaluater",
    layout="wide",
    initial_sidebar_state="expanded"
)

# add to hold_data?
def generate_sample(cohort_samples):
    st.session_state['no_plot'] = False
    sample_options = cohort_samples[~cohort_samples.isin(st.session_state[f'{st.session_state["gene_choice"]}_sample_seen'])]

    if len(sample_options) == 0:
        st.error("No more samples to run through!")
        st.session_state['no_plot'] = True
    else:
        st.session_state['sample_name'] = random.choice(list(sample_options)) # is list function necessary here
    # return no_plot

def plot_sample():
    sample_df_interval = pd.read_csv(f"data/{st.session_state['cohort_choice']}/{st.session_state['model_choice']}/{st.session_state['gene_choice']}/pred_cnvs/{st.session_state['sample_name']}_{st.session_state['gene_choice']}_full_interval.csv")
    sample_df_interval['ALT_pred'].fillna('<None>', inplace = True)
    pred_cnv = sample_df_interval[sample_df_interval['CNV_call'] == 1]

    # st.dataframe(sample_df_interval)
    # st.dataframe(pred_cnv)

    # predictions only
    fig_lrr = plot_clusters(pred_cnv[pred_cnv['ALT_pred']!='<INS>'], x_col='position', y_col='LogRRatio', gtype_col = 'ALT_pred', title=f'{st.session_state["gene_choice"]} Interval Predictions Only')
    xmin, xmax = pred_cnv['position'].min(), pred_cnv['position'].max()
    st.plotly_chart(fig_lrr)

    fig_baf = plot_clusters(pred_cnv[pred_cnv['ALT_pred']=='<INS>'], x_col='position', y_col='BAlleleFreq', gtype_col = 'ALT_pred', title=f'{st.session_state["gene_choice"]} Interval Predictions Only', xmin=xmin, xmax=xmax)
    st.plotly_chart(fig_baf)

    # just CNV predicted variants among all variants
    fig_lrr_full = plot_clusters(sample_df_interval[sample_df_interval['ALT_pred']!='<INS>'], x_col='position', gtype_col = 'ALT_pred', y_col='LogRRatio', title=f'{st.session_state["gene_choice"]} Interval')
    st.plotly_chart(fig_lrr_full)
    
    fig_baf_full = plot_clusters(sample_df_interval[(sample_df_interval['ALT_pred']!='<DEL>') & (sample_df_interval['ALT_pred']!='<DUP>')], x_col='position', gtype_col = 'ALT_pred', y_col='BAlleleFreq', title=f'{st.session_state["gene_choice"]} Interval')
    st.plotly_chart(fig_baf_full)

### Create sidebar options
# add session state methods later for these selections to remember when exploring other pages (if more pages added)
# split GP2 by cohort but also provide full release option
st.sidebar.markdown('### Choose a GP2 Cohort')
cohorts = ['CATPD', 'CORIELL', 'GP2']
cohort_name = st.sidebar.selectbox(label = 'Cohort Selection', label_visibility = 'collapsed', options=cohorts)

# will expand on model selection - need to distinguish test results based on model version
st.sidebar.markdown('### Choose a Model')
models = ['Preliminary Deletion Model']
model_name = st.sidebar.selectbox(label = 'Model Selection', label_visibility = 'collapsed', options=models)

# will expand on this selection to all NDDs
st.sidebar.markdown('### Choose an NDD-Related Gene')
genes = ['PRKN', 'SNCA', '22q']
gene_name = st.sidebar.selectbox(label = 'NDD-Related Gene Selection', label_visibility = 'collapsed', options=genes)

models_dict = {'Preliminary Deletion Model': 'prelim_del_model'}


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
model_results = pd.read_csv(f'data/{cohort_name}/{models_dict[model_name]}/{gene_name}/{gene_name}_{cohort_name}_app_ready.csv')

with st.expander("Filter Displayed Samples"):
    confidence = st.select_slider('Display samples with prediction probability of at least:', options=[0.8, 0.9, 1], value = 1)

    # adjust these thresholds 
    iqr_range = np.linspace(min(abs(model_results['abs_iqr_lrr'])), max(abs(model_results['abs_iqr_lrr'])), num=50)
    iqr_threshold = st.select_slider('Maximum Absolute Value LRR range threshold:', options=iqr_range, value = max(iqr_range))

    # make options list a range, not discrete options
    if st.session_state['gene_choice'] == '22q':
        # use standard deviation for ranges instead of hard coded ranges
        lower_range = np.linspace(min(model_results['cnv_range_count']), max(model_results['cnv_range_count']-2500), num=50)
        upper_range = np.linspace(max(model_results['cnv_range_count']-2500), max(model_results['cnv_range_count']), num=50)
        selected_value_lower = min(lower_range)
        selected_value_higher = max(upper_range)
    else:
        lower_range = range(40, 110, 10)
        upper_range = range(100, 500, 50)
        selected_value_lower = 80
        selected_value_higher = 300

    lower_range_threshold = st.select_slider('Minimum threshold for CNV count:', options=lower_range, value = selected_value_lower)
    upper_range_threshold = st.select_slider('Maximum threshold for CNV count :', options=upper_range, value = selected_value_higher)

    exp1, exp2, exp3 = st.columns([1.5,1,1])
    threshold_submit  = exp2.checkbox('Plot With Thresholds')

if threshold_submit:
    # may need to add 'submit' button bc of delay with generating IIDs with these metrics
    threshold_results = model_results[(model_results['abs_iqr_lrr'] <= iqr_threshold) & (model_results['cnv_range_count'] >= lower_range_threshold) & 
                                    (model_results['cnv_range_count'] <= upper_range_threshold) & (model_results['Pred Values'] >= confidence)]
    cohort_samples = threshold_results['IID']
else:
    # cohort_samples = model_results.IID[model_results['Pred Values'] >= confidence]
    cohort_samples = model_results.IID[model_results['Pred Values'] == 1]

if f'{st.session_state["gene_choice"]}_sample_seen' not in st.session_state:
    st.session_state[f'{st.session_state["gene_choice"]}_sample_seen'] = []
if 'yes_choices' not in st.session_state:
    st.session_state['yes_choices'] = []
if 'maybe_choices' not in st.session_state:
    st.session_state['maybe_choices'] = []
if 'no_choices' not in st.session_state:
    st.session_state['no_choices'] = []
if 'sample_name' not in st.session_state:
    no_plot = generate_sample(cohort_samples)
if option_change:
    no_plot = generate_sample(cohort_samples)

col1, col2 = st.columns([3,0.7])
btn0, btn1, btn2, btn3, btn4 = st.columns([1,0.5,0.5,0.5,0.5])

yes = btn1.button('Yes')
maybe = btn2.button('Maybe')
no_btn = btn3.button('No')

if yes:
    st.session_state['yes_choices'].append(st.session_state['sample_name'])
    st.session_state[f'{st.session_state["gene_choice"]}_sample_seen'].append(st.session_state['sample_name'])
    no_plot = generate_sample(cohort_samples)
elif maybe:
    st.session_state['maybe_choices'].append(st.session_state['sample_name'])
    st.session_state[f'{st.session_state["gene_choice"]}_sample_seen'].append(st.session_state['sample_name'])
    no_plot = generate_sample(cohort_samples)
elif no_btn:
    st.session_state['no_choices'].append(st.session_state['sample_name'])
    st.session_state[f'{st.session_state["gene_choice"]}_sample_seen'].append(st.session_state['sample_name'])
    no_plot = generate_sample(cohort_samples)

col1.markdown(f'##### _Would you consider Sample {st.session_state["sample_name"]} a structural variant?_')
col2.markdown(f'Prediction probability of {str(round(model_results.loc[model_results.IID == st.session_state["sample_name"], "Pred Values"].iloc[0], 2))}')

if not st.session_state['no_plot']:
    plot_sample()

# Add download button and make choices into dataframes/dictionaries that include gene name where CNV was found/not found

with st.sidebar.expander("View Reported Samples"):
    st.data_editor(st.session_state['yes_choices'],
                    column_config={"value": st.column_config.TextColumn("'Yes' Choices")},
                    hide_index=True,
                    use_container_width=True
                )
    st.data_editor(st.session_state['maybe_choices'],
                    column_config={"value": st.column_config.TextColumn("'Maybe' Choices")},
                    hide_index=True,
                    use_container_width=True
                )
    st.data_editor(st.session_state['no_choices'],
                    column_config={"value": st.column_config.TextColumn("'No' Choices")},
                    hide_index=True,
                    use_container_width=True
                )

