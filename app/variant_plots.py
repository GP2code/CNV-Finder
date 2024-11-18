import numpy as np
import pandas as pd
import plotly.express as px


def plot_variants(df, x_col='BAlleleFreq', y_col='LogRRatio', gtype_col='GT', title='snp plot', opacity=1, midline=False, cnvs=None, xmin=None, xmax=None):
    """
    Plots an interactive scatter plot of genetic variants with customizable axes, colors, and features.

    Arguments:
    df (pandas.DataFrame): The DataFrame containing the data to be plotted.
    x_col (str, optional): The column name for the x-axis. Defaults to 'BAlleleFreq'.
    y_col (str, optional): The column name for the y-axis. Defaults to 'LogRRatio'.
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

    lmap = {'BAlleleFreq': 'BAF', 'LogRRatio': 'LRR'}
    smap = {'Control': 'circle', 'PD': 'diamond-open-dot'}

    # Choose color map based on genotype column
    if gtype_col == 'ALT_pred' or gtype_col == 'ALT':
        cmap_choice = cmap_ALT
    else:
        cmap_choice = cmap

    if isinstance(cnvs, pd.DataFrame):
        fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map=cmap_choice,
                         color_continuous_scale=px.colors.sequential.matter, width=650, height=497, labels=lmap, symbol_map=smap, hover_data=[gtype_col])
        fig.update_traces(opacity=opacity, marker_color='grey')

        # Overlay CNV data with a specific color (in paper: #549cdc: long read, #B371BE: short read)
        fig.add_traces(px.scatter(cnvs, x=x_col, y=y_col, hover_data=[
                       gtype_col]).update_traces(marker_color="#B371BE").data)
    else:
        if gtype_col == None:
            fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_sequence=[
                             'grey'], width=650, height=497, labels=lmap, symbol_map=smap, opacity=opacity)
        else:
            fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map=cmap_choice,
                             width=650, height=497, labels=lmap, symbol_map=smap, opacity=opacity)

    if midline:
        # Calculate the average y-value for each unique x-value
        unique_x = np.linspace(min(df[x_col]), max(df[x_col]), num=50)

        # Create bins and calculate average y within each bin
        df['x_bin'] = pd.cut(df[x_col], bins=unique_x)
        grouped_df = df[[x_col, 'x_bin', y_col]].groupby(
            'x_bin').mean().reset_index()

        # Plot the midline
        fig.add_traces(px.line(grouped_df, x=x_col, y=y_col).update_traces(
            line=dict(color='red', width=3), name='Average Line').data)

    fig.update_xaxes(range=xlim, nticks=10, zeroline=False)
    fig.update_yaxes(range=ylim, nticks=10, zeroline=False)
    fig.update_layout(margin=dict(r=76, t=63, b=75))
    fig.update_layout(legend_title_text='CNV Range Class')
    fig.update_layout(title_text=f'<b>{title}<b>')

    return fig
