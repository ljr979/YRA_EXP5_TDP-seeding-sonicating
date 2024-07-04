import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger


#now want to group based on the variables and plot the lag time in a column plot with scatter over the top.
def bar(fit_params, output_folder, monomer_only, group='t1', hue_order= None, order=None,hue_var='t2', pal='BuGn', x_var='t3', y_var='tlag', savefig=False, ylabel='Lag Phase (h)', xlabel='Seed added (%)', ylim=None):
    fig, axes = plt.subplots(len(fit_params[group].unique()), 1, sharex=False, sharey=True, figsize=(5,15))
    for i, (t1, df) in enumerate(fit_params.groupby(group)):
        t1
        df
        i = i-1
        df = df[df['t3']!='seeds']
        if pal == p_dict:
          new_pal = p_dict[t1]

        if pal!= p_dict:
            new_pal = pal

        
        control_df = monomer_only
        control_df['t1'] = t1
        df1 = pd.concat([df, control_df])



        sns.barplot(df1, x=x_var, y=y_var, hue=hue_var, hue_order=hue_order, order=order, palette=new_pal, dodge=True, edgecolor='darkgray', linewidth=0.5,capsize=0.1, ax=axes[i])
        axes[i].annotate(f'{t1}', xy=(-0.4, 14.8))
        sns.stripplot(data=df1, x=x_var, y=y_var, hue=hue_var, hue_order=hue_order, order=order, dodge=True, jitter=False, palette=new_pal, edgecolor='black', alpha=0.8, linewidth=1, ax=axes[i])
        axes[i].legend_.remove()
        if ylim: # remove the second legend
            axes[i].set_ylim(0,ylim)
        axes[i].set_ylabel(ylabel)
        axes[i].set_xlabel(xlabel)
        plt.tight_layout()
    if savefig==True:
        plt.savefig(f'{output_folder}{y_var}_plot_all.svg')
        plt.savefig(f'{output_folder}{y_var}_plot_all.png')



if __name__ == "__main__":    
    input_folder = 'results/combined/'
    output_folder = 'results/combined/1_summaries/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    fit_params = pd.read_csv(f'{input_folder}collated_fit_params.csv')

    fit_params.drop([col for col in fit_params.columns.tolist() if col in ['Unnamed: 0.1', 'Unnamed: 0']], axis=1, inplace=True)

    fit_params[['t1', 't2', 't3', 'replicate']] = fit_params['sample'].str.split('-', expand=True)
    fit_params = fit_params[fit_params['rquare']>0.98]

    #to plot multiple at the same time
    params_to_plot = ['tlag', 'slope', 'V50', 'top', 'elongation_rate']
    p_dict = {
        'TDPLCD':'Purples',
        'TDPmRUBY':'Reds',
        'TDPt25':'Blues',
        'TDPt35':'Greens',
        'monomer':'YlOrRd'

    }

    #filter for some whacky values?
    fit_params=fit_params[fit_params['tlag']>0]
    fit_params=fit_params[fit_params['top']<700000]
    fit_params=fit_params[fit_params['slope']>0]
    fit_params=fit_params[fit_params['elongation_rate']<60000]

    for param in params_to_plot:
        #to plot the monomer alone on every plot
        param
        monomer_only = fit_params[fit_params['t2']=='monomer']
        monomer_only['t3'] = '0'
        monomer_only['t2'] = 'control'
        bar(fit_params, output_folder, monomer_only, hue_order = ['unsonicated', 'sonicated'], order=['0', '0.01', '0.03', '0.1', '0.3', '1'],group='t1', hue_var='t2', pal=p_dict, x_var='t3', y_var=param, savefig=True, ylabel=param, xlabel='Seed added (%)')



    #to play around with on its own

    bar(fit_params, output_folder=output_folder, monomer_only=monomer_only, group='t1', hue_var='t2', pal='BuGn', x_var='t3', y_var='elongation_rate', savefig=False, ylabel='Elongation rate (change in ThT/h)', xlabel='Seed added (%)', ylim=40000)