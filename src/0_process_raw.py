import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

input_folder = 'raw_data/rep_1/'
output_folder = 'Results/rep_1/'

if not os.path.exists(output_folder):
   os.makedirs(output_folder)

#this is the number of minutes that each time point represents
interval_mins = 
#name of the experiment you are analysing
Exp_num = 'Exp3_1-1'
#name of the file you want to be read in from this repository
filename = 'YRA_EXP3-1_TDP_fragments'

#define the palettes you want to use for the processed plotting
p_dict = {
    'normalised_to_max_change':'Purples',
    'change_in_intensity':'Reds',
}

#make sure you go into your file beforehand and remove weird columns/metadata. Save as a csv.
def read_add_time(input_folder, filename, interval_mins):
    """cleans up your raw aggregation assay data and returns it as needed.
    also adds time to the dataframe correctly 

    Args:
        input_folder (str): the location of the csv file containing your data
        filename (str): the actual name of the file you want to read in
        interval_mins (int): the number of minutes between time points


    Returns:
        _type_: df
    """
    # read in file
    df = pd.read_csv(f'{input_folder}{filename}.csv' )
    df.drop(df.index[:1], inplace=True)
    df.insert(0, 'time(min)', range(0, 0 + len(df)*interval_mins,interval_mins))
    df['time(h)'] = df['time(min)']/60
    return df

def rename_cols(df, col_name_dict):
    """renames columns according to smaple rather than well number as a column called 'ID' which is what you have defined as being in each well.

    Args:
        df (df): df you just read in
        col_name_dict (dict): the dictionary defining which samples are in which well.

    Returns:
        _type_: _description_
    """
    df.rename(columns=col_name_dict, inplace=True)
    df.drop([col for col in df.columns.tolist() if col in ['Well', 'Unnamed: 1', 'time(min)']], axis=1, inplace=True)
    df = df.T.reset_index().rename(columns={'index':'ID'})
    return df

def average_tech_replicates(df):
    new=[]
    #little loop to make some columns to identify each treatment
    for ids, d in df.groupby('ID'):    
        if 'time(h)' not in ids:
            d['treat2'] = ids.split('_')[1]
            d['treat1'] = ids.split('_')[0]
            d['replicate'] = ids.split('_')[-1]
            new.append(d)
    new=pd.concat(new)
    
    #now average the numeric data between technical replicates
    averages = []
    for treat, x in new.groupby('treat2'):
        x
        to_av=[col for col in x if col not in ['ID', 'replicate', 'treat1', 'treat2']]
        numeric = x[[c for c in x.columns if c in to_av]].astype(int)
        #averge the values between replicates
        average = pd.DataFrame(numeric.mean(axis=0)).T
        #add back in the treatment column to identify the data
        average['treat'] = treat
        averages.append(average)
    averages=pd.concat(averages)
    return averages

def subtract_background(averages):
    #identify which is blank to subtract
    blanks = averages[averages['treat']=='blank']
    #transpose and remove the time col
    blanks = blanks.T[:-1]
    #save the treatments to add back on after doing the subtraction
    treats = averages[averages['treat']!= 'blank']
    treats_list = treats['treat'].tolist()
    #make this the same arrangement/orientation at the blank 
    treats_transp=treats.T[:-1]
    #actually subtract the blank from each of the other treatments
    subtracted = (treats_transp - blanks).set_axis(treats_list, axis=1).T.reset_index().rename(columns={'index':'ID'})
    subtracted['type'] = 'blank_subtracted'
    return subtracted, treats_list, treats_transp

def plot_fluo_over_time(d, interval_mins, output_folder, Exp_num, palette, title, version, length_exp, ylim=False, savefig=True):
    """plot signal over time

    Args:
        d (df): dataframe to plot from
        interval_mins (int): how many mins per measurement
        output_folder (str): where to save
        Exp_num (str): experiment annotation
        palette (str): colour to plot
        title (str): name of the plot
        version (str): type of data (what you've done to treat it/process it)
        length_exp (float): how long you want the x axis to be (to match other experiments or show full experiment)
        savefig (bool, optional): save the figure to output foldr or not. Defaults to True.
    """

    test_Df = pd.melt(d, id_vars= ['ID', 'type'], var_name='timepoint', value_name='intensity' )

    test_Df['time (h)'] = (test_Df['timepoint']*interval_mins)/60 
    if ylim == False:
        ylim = (min(test_Df['intensity']), max(test_Df['intensity']))

    Fig, ax = plt.subplots()
    sns.lineplot(data=test_Df, x='time (h)',y='intensity', palette=palette, hue='ID', ax=ax)
    ax.set_xlabel('Time (h)')
    ax.set_xlim(0,length_exp)
    ax.set_ylim(ylim)
    ax.set_ylabel('ThT intensity (A.U.)')
    plt.legend()
    plt.title(f'{title}')
    if savefig == True:
        plt.savefig(f'{output_folder}{Exp_num}_timeseries_{version}.png')
    plt.show()

def change_over_time_normalise(treats_transp):
    change_in_intensity = []
    prop_of_total = []
    for col, vals in treats_transp.items():
        col
        #loop over each column (i.e. average intensity minus blank, for each treatment)
        z = pd.DataFrame(vals)
        start = z[f'{col}'].iloc[0]
        #for each value, subtract the INITIAL intensity of that treatment from every other value in the column. This gives the change over time in terms of fluorescence intensity values
        change = (z-start)
        #calculate the maximum change (in fluorescence values)
        x = change.max()
        #now normalise this to the maximum change in fluorescence (now between 0 & 1 rather than fluorescence values)
        if x.item() > 0:
            norm_to_max = (change/x)
            norm_to_max = norm_to_max.T.reset_index().rename(columns = {'index':'ID'})
            prop_of_total.append(norm_to_max)
        #restructure for later
        change = change.T.reset_index().rename(columns = {'index':'ID'})
        change_in_intensity.append(change)

    prop_of_total = pd.concat(prop_of_total)
    #add a column that says what stage of processing these data are 
    prop_of_total['type'] = 'normalised_to_max_change'
    change_in_intensity = pd.concat(change_in_intensity)
    change_in_intensity['type'] = 'change_in_intensity'
    #now combine them
    alls=pd.concat([change_in_intensity, prop_of_total])
    return alls

df = read_add_time(input_folder, filename, interval_mins)

#stop- here you can run this next line to get a list of the wells, and then you can add in what they're called to the definition of the dictionary.
df.columns.tolist()
#fill this dictionary to rename the columns with the 'well' number from the plate, renamed to your sample ID (this should be in the format: for blank, just id_replicate, for others, have treatment_identifiers(separated by '-')_replicate in plate)
col_name_dict = { 
 'B07':'blank_blank_1',
 'B08':'blank_blank_2',
 'B09':'blank_blank_3',
 'J07':"fragments_TDP-25-6M_1",
 'J08':"fragments_TDP-25-6M_2",
 'J09':"fragments_TDP-25-6M_3",
 'K07':"fragments_TDP-25-2M_1",
 'K08':"fragments_TDP-25-2M_2",
 'K09':"fragments_TDP-25-2M_3",
 'L07':"fragments_TDP-25-0.6M_1",
 'L08':"fragments_TDP-25-0.6M_2",
 'L09':"fragments_TDP-25-0.6M_3",
 'M07':"fragments_TDP-25-0.2M_1",
 'M08':"fragments_TDP-25-0.2M_2",
 'M09':"fragments_TDP-25-0.2M_3",
 'N07':"fragments_TDP-35-6M_1",
 'N08':"fragments_TDP-35-6M_2",
 'N09':"fragments_TDP-35-6M_3",
 'O07':"fragments_TDP-35-2M_1",
 'O08':"fragments_TDP-35-2M_2",
 'O09':"fragments_TDP-35-2M_3",
 'B11':"fragments_TDP-35-0.6M_1",
 'B12':"fragments_TDP-35-0.6M_2",
 'B13':"fragments_TDP-35-0.6M_3",
 'C11':"fragments_TDP-35-0.2M_1",
 'C12':"fragments_TDP-35-0.2M_2",
 'C13':"fragments_TDP-35-0.2M_3",
 'D11':"fragments_TDPLCD-0M_1",
 'D12':"fragments_TDPLCD-0M_2",
 'D13':"fragments_TDPLCD-0M_3",
 'E11':"fragments_GUAN-cont-2M_1",
 'E12':"fragments_GUAN-cont-2M_2",
 'E13':"fragments_GUAN-cont-2M_3",
 'F11':"fragments_TDPLCD-0.2M_1",
 'F12':"fragments_TDPLCD-0.2M_2",
 'F13':"fragments_TDPLCD-0.2M_3"}
df = rename_cols(df, col_name_dict)


time = df[df['ID']=='time(h)']

averages = average_tech_replicates(df)

subtracted, treats_list, treats_transp = subtract_background(averages)

plot_fluo_over_time(d=subtracted, interval_mins=interval_mins, output_folder=output_folder, Exp_num=Exp_num, palette='viridis', title='raw_data', version='raw', length_exp=50, savefig=True)

#then, want to calculate CHANGE over time for each treatment, and plot again.
#make sure the columns are labelled according to the treatment
treats_transp.columns = treats_list
alls = change_over_time_normalise(treats_transp)

#define the palette for plotting
for types, a in alls.groupby('type'):
    palette = p_dict[f'{types}']
    plot_fluo_over_time(d=a, interval_mins=interval_mins, output_folder=output_folder, Exp_num=Exp_num, palette=palette, title=f'{types}', version=types, length_exp=50, ylim=False, savefig=True)


#now plot a zoomed axis
for types, a in alls.groupby('type'):
    palette = p_dict[f'{types}']
    plot_fluo_over_time(d=a, interval_mins=interval_mins, output_folder=output_folder, Exp_num=Exp_num, palette=palette, title=f'{types}', version=types, length_exp=50, ylim=(0,50),savefig=False)
