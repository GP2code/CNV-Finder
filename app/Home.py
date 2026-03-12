import os
import random
import numpy as np
import pandas as pd
import streamlit as st
import plotly.express as px
import plotly.io as pio
from datetime import datetime

from variant_plots import plot_variants

# Supress Pandas copy warning
pd.options.mode.chained_assignment = None

# Streamlit configuration for main page
st.set_page_config(
    page_title="CNV Prediction Evaluator",
    layout="wide",
    initial_sidebar_state="expanded"
)

def generate_sample(cohort_samples):
    """
    Generates a random sample that has not been previously seen from a cohort of samples.
    If no new samples are available, an error message is displayed, and a flag 
    is set to indicate that no plot should be generated.

    Arguments:
    cohort_samples (pandas.Series or list): A collection of sample identifiers 
                                            from which a new sample is selected.

    Returns:
    None: Updates the Streamlit session state by setting the `sample_name` with a 
          new sample or by setting `no_plot` to True if no samples remain.
    """
    st.session_state['no_plot'] = False
    sample_options = cohort_samples[~cohort_samples.isin(
        st.session_state[samples_seen])]

    if len(sample_options) == 0:
        st.error("No more samples to run through with your current selected options and thresholds! Adjust in the drop-down menu above.")
        st.session_state['no_plot'] = True
    else:
        st.session_state['sample_name'] = random.choice(list(sample_options))


def plot_sample(exon_file=None, highlight_manual=False):
    """
    Plots various scatter plots for LRR and BAF visualizations
    for a selected sample.

    Arguments:
    exon_file: Path to exon intervals CSV file or uploaded file object.
    highlight_manual: Whether to enable manual region highlighting.

    Returns:
    Displays interactive Plotly charts in a Streamlit app.
    """
    sample_df_path = f"testing/app_ready/{st.session_state['cohort_choice']}/{st.session_state['model_choice']}/{st.session_state['gene_choice']}/pred_cnvs/{st.session_state['sample_name']}_full_interval.csv"
    sample_df_interval = pd.read_csv(sample_df_path)
    sample_df_interval['ALT_pred'] = sample_df_interval['ALT_pred'].fillna('<None>')
    pred_cnv = sample_df_interval[sample_df_interval['CNV_call'] == 1]

    # Get position range for highlighting controls
    pos_min, pos_max = int(sample_df_interval['POS'].min()), int(sample_df_interval['POS'].max())

    # Get gene interval from reference file for default highlighting
    genes_df = pd.read_csv('ref_files/glist_hg38_intervals.csv')
    gene_interval = genes_df[genes_df['NAME'] == st.session_state['gene_choice']]

    if len(gene_interval) > 0:
        default_start = int(gene_interval.iloc[0]['START'])
        default_stop = int(gene_interval.iloc[0]['STOP'])
        # Ensure defaults are within the data range
        default_start = max(pos_min, default_start)
        default_stop = min(pos_max, default_stop)
    else:
        default_start = pos_min
        default_stop = pos_max

    # Add position inputs to sidebar if manual highlighting is enabled
    if highlight_manual and exon_file is None:
        st.sidebar.info(f"Defaulting to {st.session_state['gene_choice']} gene boundaries")
        start_pos = st.sidebar.number_input(
            "Start position (bp)",
            min_value=pos_min,
            max_value=pos_max,
            value=default_start,
            step=1000,
            key='highlight_start'
        )
        stop_pos = st.sidebar.number_input(
            "Stop position (bp)",
            min_value=pos_min,
            max_value=pos_max,
            value=default_stop,
            step=1000,
            key='highlight_stop'
        )
        st.sidebar.markdown('---')

    # Plot variants in CNV ranges only (LRR vs. POS in base pairs)
    fig_lrr = plot_variants(
            pred_cnv[pred_cnv['ALT_pred'] != '<INS>'], 
            x_col='POS',
            y_col='LRR',
            gtype_col='ALT_pred',
            title=f'{st.session_state["gene_choice"]} Interval CNV Predictions Only'
    )
    xmin, xmax = pred_cnv['POS'].min(), pred_cnv['POS'].max()
    st.plotly_chart(fig_lrr)

    # Plot variants in CNV ranges only (BAF vs. POS in base pairs)
    fig_baf = plot_variants(
            pred_cnv[pred_cnv['ALT_pred'] == '<INS>'],
            x_col='POS',
            y_col='BAF',
            gtype_col='ALT_pred',
            title=f'{st.session_state["gene_choice"]} Interval CNV Predictions Only',
            xmin=xmin, xmax=xmax
    )
    st.plotly_chart(fig_baf)

    # Plot all variants color-coded by CNV Type (LRR)
    fig_lrr_full = plot_variants(
                sample_df_interval[sample_df_interval['ALT_pred'] != '<INS>'],
                x_col='POS',
                y_col='LRR',
                gtype_col='ALT_pred',
                title=f'{st.session_state["gene_choice"]} Interval Colored by CNV Type',
                xmin=xmin, xmax=xmax
    )
    st.plotly_chart(fig_lrr_full)

    # Plot all variants color-coded by CNV Type (BAF)
    fig_baf_full = plot_variants(
                sample_df_interval[(sample_df_interval['ALT_pred'] != '<DEL>') & (sample_df_interval['ALT_pred'] != '<DUP>')],
                x_col='POS',
                y_col='BAF',
                gtype_col='ALT_pred',
                title=f'{st.session_state["gene_choice"]} Interval Colored by CNV Type',
                xmin=xmin, xmax=xmax
    )
    st.plotly_chart(fig_baf_full)

    # Plot all variants in black & white with average line (LRR)
    bw_lrr_full = plot_variants(
                sample_df_interval,
                x_col='POS',
                y_col='LRR',
                gtype_col=None,
                midline=True,
                title=f'All Variants in {st.session_state["gene_choice"]} Interval with Average Line',
                opacity=0.3,
                xmin=xmin, xmax=xmax
    )
    st.plotly_chart(bw_lrr_full)

    # Plot all variants in black & white (BAF)
    bw_baf_full = plot_variants(
                sample_df_interval,
                x_col='POS',
                y_col='BAF',
                gtype_col=None,
                title=f'All Variants in {st.session_state["gene_choice"]} Interval',
                opacity=0.3,
                xmin=xmin, xmax=xmax
    )
    st.plotly_chart(bw_baf_full)

    # Plot with highlighted exon regions if enabled
    if exon_file is not None:
        try:
            # Read exon intervals
            if isinstance(exon_file, str):
                exons_df = pd.read_csv(exon_file)
            else:
                exons_df = pd.read_csv(exon_file)

            # Create boolean mask for variants within any exon region
            in_exon = pd.Series(False, index=sample_df_interval.index)
            for _, exon in exons_df.iterrows():
                in_exon |= (
                    (sample_df_interval['POS'] >= exon['Start']) &
                    (sample_df_interval['POS'] <= exon['Stop'])
                )

            highlighted_exons = sample_df_interval[in_exon]
            remaining_variants = sample_df_interval[~in_exon]

            if len(highlighted_exons) > 0:
                st.markdown("---")
                st.info(f"Highlighting {len(highlighted_exons)} variants across {len(exons_df)} exon regions")

                # Plot LRR with highlighted exons
                highlight_lrr = plot_variants(
                    remaining_variants,
                    x_col='POS',
                    y_col='LRR',
                    gtype_col=None,
                    title=f'{st.session_state["gene_choice"]} Exons Highlighted',
                    opacity=0.3,
                    cnvs=highlighted_exons,
                    xmin=xmin, xmax=xmax
                )
                st.plotly_chart(highlight_lrr)

                # Plot BAF with highlighted exons
                highlight_baf = plot_variants(
                    remaining_variants,
                    x_col='POS',
                    y_col='BAF',
                    gtype_col=None,
                    title=f'{st.session_state["gene_choice"]} Exons Highlighted',
                    opacity=0.3,
                    cnvs=highlighted_exons,
                    xmin=xmin, xmax=xmax
                )
                st.plotly_chart(highlight_baf)
            else:
                st.warning("No variants found within the specified exon regions")

        except Exception as e:
            st.error(f"Error reading exon intervals: {str(e)}")

    # Plot with manual highlighted region if enabled (and no exon file)
    elif highlight_manual:
        if start_pos >= stop_pos:
            st.warning("Start position must be less than stop position")
        else:
            # Split data into highlighted region and remaining variants
            highlighted_region = sample_df_interval[(sample_df_interval['POS'] >= start_pos) & (sample_df_interval['POS'] <= stop_pos)]
            remaining_variants = sample_df_interval[~((sample_df_interval['POS'] >= start_pos) & (sample_df_interval['POS'] <= stop_pos))]

            st.markdown("---")
            st.info(f"Highlighting {len(highlighted_region)} variants in region {start_pos:,} - {stop_pos:,} bp")

            # Plot LRR with highlighted region
            highlight_lrr = plot_variants(
                remaining_variants,
                x_col='POS',
                y_col='LRR',
                gtype_col=None,
                title=f'All Variants in {st.session_state["gene_choice"]} Interval with Highlighted Region',
                opacity=0.3,
                cnvs=highlighted_region,
                xmin=xmin, xmax=xmax
            )
            st.plotly_chart(highlight_lrr)

            # Plot BAF with highlighted region
            highlight_baf = plot_variants(
                remaining_variants,
                x_col='POS',
                y_col='BAF',
                gtype_col=None,
                title=f'All Variants in {st.session_state["gene_choice"]} Interval with Highlighted Region',
                opacity=0.3,
                cnvs=highlighted_region,
                xmin=xmin, xmax=xmax
            )
            st.plotly_chart(highlight_baf)


# Creates sidebar for model selections
st.sidebar.markdown('### Choose a model:')
models_dict = {'Preliminary Deletion Model': 'prelim_del_model', 'Preliminary Duplication Model': 'prelim_dup_model',
               'Updated Deletion Model': 'updated_del_model', 'Updated Duplication Model': 'updated_dup_model',
               'Final Deletion Model': 'final_del_model', 'Final Duplication Model': 'final_dup_model'}
model_name = st.sidebar.selectbox(
    label='Model Selection', label_visibility='collapsed', options=list(models_dict.keys()), index=4)

# Cohort selection
st.sidebar.markdown('### Choose a cohort:')
cohorts = next(os.walk('testing/app_ready/'))[1]
cohort_name = st.sidebar.selectbox(
    label='Cohort Selection', label_visibility='collapsed', options=sorted(cohorts))

# Gene of interest selection (load from reference file)
st.sidebar.markdown('### Choose an NDD-related gene:')
genes_df = pd.read_csv('ref_files/glist_hg38_intervals.csv')
genes = sorted(genes_df['NAME'].unique())
default_index = genes.index('PARK2') if 'PARK2' in genes else 0
gene_name = st.sidebar.selectbox(
    label='NDD-Related Gene Selection', label_visibility='collapsed', options=genes, index=default_index)

# Region highlighting controls
st.sidebar.markdown('---')
st.sidebar.markdown('### Highlight Specific Region:')

# Exon intervals file upload
exon_file_default = f'ref_files/exons/{gene_name}_exons.csv'
use_exons = st.sidebar.checkbox("Highlight exon intervals", value=False)

if use_exons:
    if os.path.isfile(exon_file_default):
        st.sidebar.info(f"Using default: {gene_name}_exons.csv")
        exon_file = exon_file_default
    else:
        st.sidebar.warning(f"No default exon file found for {gene_name}")
        uploaded_file = st.sidebar.file_uploader("Upload exon intervals CSV", type=['csv'])
        if uploaded_file is not None:
            exon_file = uploaded_file
        else:
            exon_file = None
    highlight_manual = False
else:
    exon_file = None
    # Manual region highlighting
    highlight_manual = st.sidebar.checkbox("Highlight custom region", value=False)

# Adjusts session state variables based on sidebar selection changes
# Track configuration to prevent unnecessary sample regeneration
current_config = f"{gene_name}_{models_dict[model_name]}_{cohort_name}"
if 'last_config' not in st.session_state:
    st.session_state['last_config'] = current_config
    option_change = False
else:
    option_change = (current_config != st.session_state['last_config'])
    if option_change:
        st.session_state['last_config'] = current_config

if 'threshold_submit' not in st.session_state:
    st.session_state['threshold_submit'] = False

# Initializes session state variables
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

# On-change function for adjusting predicted value threshold
def threshold_true():
    st.session_state['threshold_submit'] = True

# Defines variables necessary for pulling results to plot
st.session_state['cohort_choice'] = cohort_name
st.session_state['model_choice'] = models_dict[model_name]
st.session_state['gene_choice'] = gene_name

# Populates Main Page
st.title('Evaluation of CNV Predictions')
model_path = f'testing/app_ready/{cohort_name}/{models_dict[model_name]}/{gene_name}/{cohort_name}_{gene_name}_app_ready.csv'

# If app-ready file exists, display sample plots
if not os.path.isfile(model_path):
    st.error('No CNVs to display!')
else:
    model_results = pd.read_csv(model_path)

    if len(model_results) == 0:
        st.error('No CNVs to display!')
    else:
        with st.expander("Filter Displayed Samples"):
            # Determine model type for probability intervals
            model_type = st.session_state['model_choice'].split('_')[1] if '_' in st.session_state['model_choice'] else st.session_state['model_choice']

            if model_type == 'del':
                prob_interval_options = {
                    '≥ 0.6 and < 0.7': (0.6, 0.7),
                    '≥ 0.7 and < 0.8': (0.7, 0.8),
                    '≥ 0.8 and < 0.9': (0.8, 0.9),
                    '≥ 0.9 and ≤ 1.0': (0.9, 1.0),
                }
            else:
                prob_interval_options = {
                    '≥ 0.5 and < 0.6': (0.5, 0.6),
                    '≥ 0.6 and < 0.7': (0.6, 0.7),
                    '≥ 0.7 and < 0.8': (0.7, 0.8),
                    '≥ 0.8 and < 0.9': (0.8, 0.9),
                    '≥ 0.9 and ≤ 1.0': (0.9, 1.0),
                }

            # Calculate and display counts for each interval
            interval_counts = []
            counts_by_index = []
            for label, (low, high) in prob_interval_options.items():
                upper_mask = (model_results['Pred Values'] <= high) if high == 1.0 else (
                    model_results['Pred Values'] < high)
                count = len(
                    model_results[(model_results['Pred Values'] >= low) & upper_mask])
                interval_counts.append(
                    {'Interval': label, 'Predictions': count})
                counts_by_index.append(count)
            st.dataframe(pd.DataFrame(interval_counts),
                         hide_index=True, use_container_width=True)

            # Select default interval or closest non-empty one
            prob_interval_labels = list(prob_interval_options.keys())
            default_prob_index = len(prob_interval_labels) - 1
            selected_prob_index = default_prob_index
            if counts_by_index[default_prob_index] == 0:
                candidates = [
                    (abs(idx - default_prob_index), idx)
                    for idx, cnt in enumerate(counts_by_index)
                    if cnt > 0
                ]
                if candidates:
                    _, selected_prob_index = min(candidates)
                    st.info(
                        f"No samples found in {prob_interval_labels[default_prob_index]}; using {prob_interval_labels[selected_prob_index]} instead."
                    )

            prob_interval_label = st.selectbox(
                'Select prediction probability interval:',
                options=prob_interval_labels,
                index=selected_prob_index,
                on_change=threshold_true
            )
            prob_low, prob_high = prob_interval_options[prob_interval_label]

            # Adjustable slider for IQR filtering
            iqr_range = np.linspace(min(abs(model_results['abs_iqr_lrr'])), max(
                abs(model_results['abs_iqr_lrr'])), num=50)
            iqr_threshold = st.select_slider('Maximum Absolute Value LRR range threshold:', options=iqr_range, value=max(
                iqr_range), format_func=lambda x: "{:.2f}".format(x), on_change=threshold_true())

            # Adjustable slider for filtering based on CNV candidate counts
            min_cnv_count = min(model_results['cnv_range_count'])
            max_cnv_count = max(model_results['cnv_range_count'])
            lower_range = np.linspace(
                min_cnv_count, max_cnv_count - (0.5 * max_cnv_count), num=50, dtype=int)
            upper_range = np.linspace(
                max_cnv_count - (0.5 * max_cnv_count), max_cnv_count, num=50, dtype=int)
            selected_value_lower = min(lower_range)
            selected_value_higher = max(upper_range)

            lower_range_threshold = st.select_slider(
                'Minimum count of variants in CNV range:', options=lower_range, value=selected_value_lower, on_change=threshold_true())
            upper_range_threshold = st.select_slider(
                'Maximum count of variants in CNV range:', options=upper_range, value=selected_value_higher, on_change=threshold_true())

        # Always use the selected interval filters
        st.session_state['threshold_submit'] = True
        prob_high_inclusive = prob_high == 1.0
        threshold_results = model_results[
            (model_results['abs_iqr_lrr'] <= iqr_threshold)
            & (model_results['cnv_range_count'] >= lower_range_threshold)
            & (model_results['cnv_range_count'] <= upper_range_threshold)
            & (model_results['Pred Values'] >= prob_low)
            & (
                (model_results['Pred Values'] < prob_high)
                if not prob_high_inclusive
                else (model_results['Pred Values'] <= prob_high)
            )
        ]

        cohort_samples = threshold_results['IID']

        # Track threshold/filter changes to regenerate sample when filters change
        current_threshold_config = f"{prob_interval_label}_{iqr_threshold}_{lower_range_threshold}_{upper_range_threshold}"
        if 'last_threshold_config' not in st.session_state:
            st.session_state['last_threshold_config'] = current_threshold_config
            threshold_change = False
        else:
            threshold_change = (current_threshold_config != st.session_state['last_threshold_config'])
            if threshold_change:
                st.session_state['last_threshold_config'] = current_threshold_config

        # Variable to hold samples seen within the specified gene, cohort, and model selection
        samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'

        # Pull a new sample to display when the app starts and when a change in gene, cohort, and model is made
        if samples_seen not in st.session_state:
            st.session_state[samples_seen] = []

        # Initialize flag for button press detection
        if 'button_pressed' not in st.session_state:
            st.session_state['button_pressed'] = False

        # Generate new sample on first load
        if 'sample_name' not in st.session_state:
            generate_sample(cohort_samples)
        # Generate new sample when configuration changes
        elif option_change or threshold_change:
            generate_sample(cohort_samples)
        # Generate new sample if button was pressed in previous run
        elif st.session_state['button_pressed']:
            generate_sample(cohort_samples)
            st.session_state['button_pressed'] = False

        col1, col2 = st.columns([3, 0.7])
        btn1, btn2, btn3, btn4 = st.columns([0.5, 0.5, 0.5, 0.5])

        # Display any initial warnings
        if 'Artifact Warning' in model_results.columns and model_results['Artifact Warning'].iloc[0] == 1:
            col1.error('Please note that a high number of samples in this gene were predicted to have CNVs, which may indicate that an artifact or other array-based issue is displayed.')

        col1.markdown(
            f'##### _Would you consider Sample {st.session_state["sample_name"]} a structural variant?_')

        display_probab = model_results.loc[model_results.IID == st.session_state["sample_name"], "Pred Values"]
        if len(display_probab) > 0:
            col2.markdown(f'Prediction probability of {str(round(display_probab.iloc[0], 2))}')

        # Add first sample when app is started or when predicted value threshold is changed to samples_seen
        if len(st.session_state[samples_seen]) == 0:
            st.session_state[samples_seen].append(
                st.session_state['sample_name'])
        if st.session_state['threshold_submit']:
            st.session_state[samples_seen].append(
                st.session_state['sample_name'])

        # User can keep interacting with the buttons until all samples have been seen
        if not st.session_state['no_plot']:
            yes = btn1.button('Yes', use_container_width=True)
            maybe = btn2.button('Maybe', use_container_width=True)
            no_btn = btn3.button('No', use_container_width=True)

            # Currently only works for deletions and duplications, no custom entry
            other_cnv = btn4.button('Other CNV', use_container_width=True)
            plot_sample(exon_file=exon_file, highlight_manual=highlight_manual)
        else:
            yes = btn1.button('Yes', disabled=True, use_container_width=True)
            maybe = btn2.button('Maybe', disabled=True,
                                use_container_width=True)
            no_btn = btn3.button('No', disabled=True, use_container_width=True)
            other_cnv = btn4.button(
                'Other CNV', disabled=True, use_container_width=True)

        # Create report for sample selections that can be exported
        if yes:
            st.session_state['button_pressed'] = True
            samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'
            st.session_state[samples_seen].append(
                st.session_state['sample_name'])

            if st.session_state['threshold_submit']:
                st.session_state['yes_choices'].append(
                    st.session_state[samples_seen][-3])
            else:
                st.session_state['yes_choices'].append(
                    st.session_state[samples_seen][-2])

            st.session_state['yes_gene'].append(
                st.session_state["gene_choice"])
            cnv_type = st.session_state["model_choice"].split('_')[1]
            st.session_state['yes_type'].append(cnv_type)
        elif maybe:
            st.session_state['button_pressed'] = True
            samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'
            st.session_state[samples_seen].append(
                st.session_state['sample_name'])

            if st.session_state['threshold_submit']:
                st.session_state['maybe_choices'].append(
                    st.session_state[samples_seen][-3])
            else:
                st.session_state['maybe_choices'].append(
                    st.session_state[samples_seen][-2])

            st.session_state['maybe_gene'].append(
                st.session_state["gene_choice"])
            cnv_type = st.session_state["model_choice"].split('_')[1]
            st.session_state['maybe_type'].append(cnv_type)
        elif no_btn:
            st.session_state['button_pressed'] = True
            samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'
            st.session_state[samples_seen].append(
                st.session_state['sample_name'])

            if st.session_state['threshold_submit']:
                st.session_state['no_choices'].append(
                    st.session_state[samples_seen][-3])
            else:
                st.session_state['no_choices'].append(
                    st.session_state[samples_seen][-2])

            st.session_state['no_gene'].append(st.session_state["gene_choice"])
            cnv_type = st.session_state["model_choice"].split('_')[1]
            st.session_state['no_type'].append(cnv_type)
        elif other_cnv:
            st.session_state['button_pressed'] = True
            samples_seen = f'{st.session_state["gene_choice"]}_{st.session_state["cohort_choice"]}_{st.session_state["model_choice"]}sample_seen'
            st.session_state[samples_seen].append(
                st.session_state['sample_name'])

            if st.session_state['threshold_submit']:
                st.session_state['yes_choices'].append(
                    st.session_state[samples_seen][-3])
            else:
                st.session_state['yes_choices'].append(
                    st.session_state[samples_seen][-2])

            st.session_state['yes_gene'].append(
                st.session_state["gene_choice"])

            cnv_type = st.session_state["model_choice"].split('_')[1]
            if cnv_type == 'del':
                st.session_state['yes_type'].append('dup')
            elif cnv_type == 'dup':
                st.session_state['yes_type'].append('del')
            else:
                st.session_state['yes_type'].append('unknown')

        if yes or maybe or no_btn or other_cnv:
            st.rerun()

side_btn1, side_btn2, side_btn3 = st.sidebar.columns([0.5, 1, 0.5])

yes_report = pd.DataFrame(
    {'Yes Samples': st.session_state['yes_choices'], 'Interval': st.session_state['yes_gene'], 'Type': st.session_state['yes_type']})
maybe_report = pd.DataFrame(
    {'Maybe Samples': st.session_state['maybe_choices'], 'Interval': st.session_state['maybe_gene'], 'Type': st.session_state['maybe_type']})
no_report = pd.DataFrame(
    {'No Samples': st.session_state['no_choices'], 'Interval': st.session_state['no_gene'], 'Type': st.session_state['no_type']})

# Built-in CSV download button in Streamlit 
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
