import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from lifelines import KaplanMeierFitter 

input_folder = 'results/rep_1/1_fitting_curves/'
output_folder = 'results/rep_1/2_kaplan-meier/'


if not os.path.exists(output_folder):
   os.makedirs(output_folder)

fits = pd.read_csv(f'{input_folder}fit_parameters.csv')

no_fit_treatments = pd.read_csv(f'{input_folder}no_fit.csv')

timeseries = pd.read_csv(f'{input_folder}/clean_for_fitting.csv')

endpoint = list(timeseries['timepoint'].unique())[-1]

#need to turn these into time and event data. The event in this case is the LAG phase ending. so the time will be the 'tlag' and the event will be 1. Then, the no_fit will have the time as the end point, and the 'no event' will be a 0.
#so need a new dataframe with sample, and tlag columns to start with.
#want to now look at the ones that don't fit so well
not_ideal = fits[fits['rquare']<0.95]
events = fits[['sample', 'tlag']][fits['rquare']>0.9]

#fill the events column of  this dataframe with a 1, as these all were fit to the sigmoid
#make a more sophisticated way in future? where I combine everything and say if rsquare > threshold and tlag < endpoint, assign 1, else then where tlag = endpoint, assign event = 0
events = fits[['sample', 'tlag']][fits['rquare']>0.95]
#fill the events column of  this dataframe with a 1, as these all were fit to the sigmoid
#make a more sophisticated way in future? where I combine everything and say if rsquare > threshold and tlag < endpoint, assign 1, else then where tlag = endpoint, assign event = 0
events['event'] = int(1)

no_fit.drop([col for col in no_fit.columns.tolist() if 'Unnamed: 0' in col], axis=1, inplace=True)
bad_fit_treatments = fits[['sample', 'tlag']][fits['rquare']<0.2]
no_fit_treatments = pd.concat([no_fit, bad_fit_treatments])
no_fit_treatments['tlag'] = endpoint
no_fit_treatments['event'] = int(0)

#mash the two dfs together
events_combined = pd.concat([events, no_fit_treatments])

events_combined[['t1', 't2', 't3', 'replicate']] = events_combined['sample'].str.split('-', expand=True)
#TRYING OUT needing to perform the analysis on the separated replicates, so the COMMON treatment info is t1, t2 and t3
survival = []
for (prot, t2, t3), df in events_combined.groupby(['t1','t2','t3']):
   df
   t2
   t3
   prot
   arr = df[['tlag', 'event']].to_numpy()
   easy = pd.DataFrame({'time': arr[:, 0], 'event': arr[:, 1]})
   T = easy['time']
   E = easy['event']
   kmf = KaplanMeierFitter()
   kmf.fit(T, E)
   s = kmf.survival_function_
   s['sample'] = f'{prot}_{t2}_{t3}'
   survival.append(s)
survival = pd.concat(survival).reset_index()
survival.to_csv(f'{output_folder}survival_function.csv')


#plots everything on the same graph

fig, ax = plt.subplots()
sns.lineplot(
    x = 'timeline',
    y = 'KM_estimate',
    hue = 'sample',
    data = survival,
    ax=ax
)
plt.savefig(f'{output_folder}survival_all.png')
#the below sectin is for if the plto wtih everything on the same graph is too much (i.e. ifthere are a lot of different concentrations AND multiple proteins etc.)


p_dict = {
     'TDPLCD':'Purples',
     'TDPmRUBY':'Reds',
     'TDPt25':'Blues',
     'TDPt35':'Greens'

}


#now plot the different types of TDP (variable t1) in separate graphs with their own subplots for sonicated vs non-sonicated

#need to turn this into a function obviously and have a keyword for truncating the x axis
survival[['t1', 't2', 't3']] = survival['sample'].str.split('_', expand=True)

for t1, xyz in survival.groupby(['t1']):
    palette = p_dict[t1[0]]
    t1
    xyz
    sub = len(xyz['t2'].unique().tolist())
    if sub == 1:
        row = 1
        col = 2
    if sub > 1: 
        row = 1
        col = sub
    #palette = p_dict[t1[0]]
    fig, axes = plt.subplots(row,col, figsize=(16,5))
    for i, (t2, abc) in enumerate(xyz.groupby('t2')):
        i
        t2
        abc
        sns.lineplot(
        x = 'timeline',
        y = 'KM_estimate',
        hue = 't3',
        data = abc, 
        palette=palette, 
        ax=axes[i]
        )


        axes[i].set_xlim(0,int(endpoint))
        axes[i].set_xlim(0,10)
        axes[i].set_xlabel('')
        axes[i].set_ylabel('Likelihood')
        axes[i].annotate(f'{t2}', xy=(0.8, 1))
        axes[i].legend()
    fig.suptitle (f'{t1[0]}')
    plt.savefig(f'{output_folder}{t1}-{t2}-kaplan_plots.png')


#------------------------
#NOW WE WANT TO COMPARE BETWEEN OTHER THINGS- TURN THIS INTO A FUNCTION EVENTUALLY (USING A VARIABLE AS THE THING THAT CHANGES INSTEAD OF COPYING AND CHANGING THE STRING LOL)