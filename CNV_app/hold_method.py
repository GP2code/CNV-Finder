import os
import sys
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.io as pio
import seaborn as sns
from PIL import Image
import datetime
from io import StringIO


def plot_variants(df, x_col='BAlleleFreq', y_col='LogRRatio', gtype_col='GT', title='snp plot', opacity = 1, midline = False, cnvs = None, xmin= None, xmax = None):
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
        if gtype_col == None:
            fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_sequence=['grey'], width=650, height=497, labels=lmap, symbol_map=smap, opacity=opacity)
        else:
            fig = px.scatter(df, x=x_col, y=y_col, color=gtype_col, color_discrete_map=cmap_choice, width=650, height=497, labels=lmap, symbol_map=smap, opacity=opacity)

    if midline:
        # Calculate the average y-value for each unique x-value
        unique_x = np.linspace(min(df[x_col]), max(df[x_col]), num=50)

        # Use cut to create bins and calculate average y within each bin
        df['x_bin'] = pd.cut(df[x_col], bins=unique_x)
        grouped_df = df[[x_col, 'x_bin', y_col]].groupby('x_bin').mean().reset_index()

        fig.add_traces(px.line(grouped_df, x=x_col, y=y_col).update_traces(line=dict(color='red', width=3), name='Average Line').data)

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

