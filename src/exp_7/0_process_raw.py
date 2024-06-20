import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

input_folder = 'raw_data/rep_1/'
output_folder = 'Results/rep_1/0_processing/'

if not os.path.exists(output_folder):
   os.makedirs(output_folder)

#this is the number of minutes that each time point represents
interval_mins = 4
#name of the experiment you are analysing
Exp_num = 'Exp5_1'
#name of the file you want to be read in from this repository
filename = 'raw_data'

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
    df.replace(r'\N', np.nan, inplace = True)
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
def add_info (df):
    new=[]
    #little loop to make some columns to identify each treatment
    for ids, d in df.groupby('ID'):    
        ids
        d
        if 'time(h)' not in ids:
            d['treat2'] = ids.split('_')[1]
            d['treat1'] = ids.split('_')[0]
            d['replicate'] = d.groupby('ID').cumcount()+ 1
            new.append(d)
    new = pd.concat(new)
    return new

def average_tech_replicates(new):
    #now average the numeric data between technical replicates
    averages = []
    for treat, x in new.groupby('treat2'):
        x
        to_av = [col for col in x if col not in ['ID', 'replicate', 'treat1', 'treat2']]
        numeric = x[[c for c in x.columns if c in to_av]].astype(int)
        #averge the values between replicates
        average = pd.DataFrame(numeric.mean(axis=0)).T
        #add back in the treatment column to identify the data
        average['treat'] = treat
        averages.append(average)
    averages = pd.concat(averages)
    return averages

def subtract_background(averages):
    #identify which is blank to subtract
    blanks = averages[averages['treat']=='blank-x-x-x']
    #blanks = blanks[blanks['replicate'] == '1'] 

    ifo_cols = ['ID', 'replicate', 'treat1', 'treat2', 'treat']
    to_av = blanks[[col for col in blanks if col not in ifo_cols]]
    numeric = blanks[[c for c in blanks.columns if c in to_av]].astype(int)

    #averge the values between replicates
    blank_av = pd.DataFrame(numeric.mean(axis=0), columns = ['blank']).T.reset_index().rename(columns={'index':'ID'})

    #save the treatments to add back on after doing the subtraction
    treats = averages[averages['treat']!= 'blank-x-x-x']
    treats['replicate'] = treats['replicate'].astype(str)
    #now keep the replicate information
    treats['ID'] = treats[['ID', 'replicate']].agg('-'.join, axis=1)
    #save the treatment names
    treats_list = treats['treat'].tolist()
    #now make sure each row is a treatment and each column is a timepoint, but remove the column containing the strings (i.e. the treatment names)
    treats = treats[[col for col in treats if col not in ['replicate', 'treat1', 'treat2', 'treat']]]

    treats = pd.concat([treats, blank_av])

    df = treats.T
    header = df.iloc[0]
    df = df[1:].rename(columns=header)
    #everything is a bloody string so need it to be a number!!!!!!!!!!!!!!!!!!!!!!!
    df = df.apply(pd.to_numeric)
    treats_transp = df
    #subtract the blanks from all other values
    subtracted = df.sub(df['blank'], axis=0).T.reset_index().rename(columns = {'index':'ID'})

    subtracted['type'] = 'blank_subtracted'


    return subtracted, treats_list, treats_transp

def plot_fluo_over_time(d, interval_mins, output_folder, palette, Exp_num, length_exp, types=['change_in_intensity'], ylim=False, savefig=True):
    """plot signal over time

    Args:
        d (df): dataframe to plot from
        interval_mins (int): how many mins per measurement
        output_folder (str): where to save
        Exp_num (str): experiment annotation
        palette (str): colour to plot

        length_exp (float): how long you want the x axis to be (to match other experiments or show full experiment)
        savefig (bool, optional): save the figure to output foldr or not. Defaults to True.
    """

    test_Df = pd.melt(d, id_vars= ['ID', 'type', 't1',  't2', 't3', 'replicate'], var_name='timepoint', value_name='intensity' )

    test_Df['time (h)'] = (test_Df['timepoint']*interval_mins)/60 
    if ylim == False:
        ylim = (min(test_Df['intensity']), max(test_Df['intensity']))

    for t1, xyz in test_Df[test_Df['type'].isin(types)].groupby(['t1']):
        print(t1)
        t1 = t1[0]
        sub = len(xyz['t2'].unique().tolist())
        if sub == 1:
            row = 1
            col = 2
        if sub > 1: 
            row = 1
            col = sub
        Fig, axes = plt.subplots(row, col, figsize=(16,5))
        if palette == p_dict:
            new_pal = p_dict[t1]

        if palette!= p_dict:
            new_pal = palette
        for i, (t2, abc) in enumerate(xyz.groupby('t2')):
            i
            t2
            abc
            sns.lineplot(data=abc, x='time (h)', y='intensity', palette=new_pal, hue='t3', ax=axes[i])
            axes[i].set_xlim(0,length_exp)
            axes[i].set_ylim(ylim)
            axes[i].annotate(f'{t2}', xy=(0.8, 1))
            axes[i].legend()
        Fig.suptitle(f'{t1}')        
        if savefig == True:
            plt.savefig(f'{output_folder}{Exp_num}_{t1}timeseries.png')
            plt.show()
    return test_Df

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

def get_more_info(df, col_to_split, what_to_split):
    df[['t1', 't2', 't3', 'replicate']] = df[f'{col_to_split}'].str.split(f'{what_to_split}',expand=True)
    return df

    
df = read_add_time(input_folder, filename, interval_mins)

#stop- here you can run this next line to get a list of the wells, and then you can add in what they're called to the definition of the dictionary.
df.columns.tolist()
#fill this dictionary to rename the columns with the 'well' number from the plate, renamed to your sample ID (this should be in the format: for blank, just id_replicate, for others, have treatment_identifiers(separated by '-')_replicate in plate)
#fill this dictionary to rename the columns with the 'well' number from the plate, renamed to your sample ID (this should be in the format: for blank, just id_replicate, for others, have treatment_identifiers(separated by '-')_replicate in plate)
col_name_dict = { 'F02':'blank_blank-x-x-x_1',
 'F03':'blank_blank-x-x-x_1',
 'F04':'blank_blank-x-x-x_1',
 'F06':'seeding_TDPLCD-sonicated-0.3_1',
 'F07':'seeding_TDPLCD-sonicated-0.3_1',
 'F08':'seeding_TDPLCD-sonicated-0.3_1',
 'F10':'seeding_TDPt25-unsonicated-0.01_1',
 'F11':'seeding_TDPt25-unsonicated-0.01_1',
 'F12':'seeding_TDPt25-unsonicated-0.01_1',
 'F14':'seeding_TDPt35-unsonicated-0.3_1',
 'F15':'seeding_TDPt35-unsonicated-0.3_1',
 'F16':'seeding_TDPt35-unsonicated-0.3_1',
 'F18':'control_TDPmRUBY-sonicated-seeds_1',
 'F19':'control_TDPmRUBY-sonicated-seeds_1',
 'F20':'control_TDPmRUBY-sonicated-seeds_1',
 'G02':'control_TDPLCD-monomer_1',
 'G03':'control_TDPLCD-monomer_1',
 'G04':'control_TDPLCD-monomer_1',
 'G06':'seeding_TDPLCD-sonicated-0.1_1',
 'G07':'seeding_TDPLCD-sonicated-0.1_1',
 'G08':'seeding_TDPLCD-sonicated-0.1_1',
 'G10':'control_TDPt25-unsonicated-seeds_1',
 'G11':'control_TDPt25-sonicated-seeds_1',
 'G12':'control_TDPt25-sonicated-seeds_1',
 'G14':'seeding_TDPt35-unsonicated-0.1_1',
 'G15':'seeding_TDPt35-unsonicated-0.1_1',
 'G16':'seeding_TDPt35-unsonicated-0.1_1',
 'G18':'seeding_TDPmRUBY-sonicated-1_1',
 'G19':'seeding_TDPmRUBY-sonicated-1_1',
 'G20':'seeding_TDPmRUBY-sonicated-1_1',
 'H02':'control_TDPLCD-unsonicated-seeds_1',
 'H03':'control_TDPLCD-unsonicated-seeds_1',
 'H04':'control_TDPLCD-unsonicated-seeds_1',
 'H06':'seeding_TDPLCD-sonicated-0.03_1',
 'H07':'seeding_TDPLCD-sonicated-0.03_1',
 'H08':'seeding_TDPLCD-sonicated-0.03_1',
 'H10':'seeding_TDPt25-sonicated-1_1',
 'H11':'seeding_TDPt25-sonicated-1_1',
 'H12':'seeding_TDPt25-sonicated-1_1',
 'H14':'seeding_TDPt35-unsonicated-0.03_1',
 'H15':'seeding_TDPt35-unsonicated-0.03_1',
 'H16':'seeding_TDPt35-unsonicated-0.03_1',
 'H18':'seeding_TDPmRUBY-sonicated-0.3_1',
 'H19':'seeding_TDPmRUBY-sonicated-0.3_1',
 'H20':'seeding_TDPmRUBY-sonicated-0.3_1',
 'I02':'seeding_TDPLCD-unsonicated-1_1',
 'I03':'seeding_TDPLCD-unsonicated-1_1',
 'I04':'seeding_TDPLCD-unsonicated-1_1',
 'I06':'seeding_TDPLCD-sonicated-0.01_1',
 'I07':'seeding_TDPLCD-sonicated-0.01_1',
 'I08':'seeding_TDPLCD-sonicated-0.01_1',
 'I10':'seeding_TDPt25-sonicated-0.3_1',
 'I11':'seeding_TDPt25-sonicated-0.3_1',
 'I12':'seeding_TDPt25-sonicated-0.3_1',
 'I14':'seeding_TDPt35-unsonicated-0.01_1',
 'I15':'seeding_TDPt35-unsonicated-0.01_1',
 'I16':'seeding_TDPt35-unsonicated-0.01_1',
 'I18':'seeding_TDPmRUBY-sonicated-0.1_1',
 'I19':'seeding_TDPmRUBY-sonicated-0.1_1',
 'I20':'seeding_TDPmRUBY-sonicated-0.1_1',
 'J02':'seeding_TDPLCD-unsonicated-0.3_1',
 'J03':'seeding_TDPLCD-unsonicated-0.3_1',
 'J04':'seeding_TDPLCD-unsonicated-0.3_1',
 'J10':'seeding_TDPt25-sonicated-0.1_1',
 'J11':'seeding_TDPt25-sonicated-0.1_1',
 'J12':'seeding_TDPt25-sonicated-0.1_1',
 'J14':'control_TDPt35-sonicated-seeds_1',
 'J15':'control_TDPt35-sonicated-seeds_1',
 'J16':'control_TDPt35-sonicated-seeds_1',
 'J18':'seeding_TDPmRUBY-sonicated-0.03_1',
 'J19':'seeding_TDPmRUBY-sonicated-0.03_1',
 'J20':'seeding_TDPmRUBY-sonicated-0.03_1',
 'K02':'seeding_TDPLCD-unsonicated-0.1_1',
 'K03':'seeding_TDPLCD-unsonicated-0.1_1',
 'K04':'seeding_TDPLCD-unsonicated-0.1_1',
 'K06':'control_TDPt25-unsonicated-seeds_1',
 'K07':'control_TDPt25-unsonicated-seeds_1',
 'K08':'control_TDPt25-unsonicated-seeds_1',
 'K10':'seeding_TDPt25-sonicated-0.03_1',
 'K11':'seeding_TDPt25-sonicated-0.03_1',
 'K12':'seeding_TDPt25-sonicated-0.03_1',
 'K14':'seeding_TDPt35-sonicated-1_1',
 'K15':'seeding_TDPt35-sonicated-1_1',
 'K16':'seeding_TDPt35-sonicated-1_1',
 'K18':'seeding_TDPmRUBY-sonicated-0.01_1',
 'K19':'seeding_TDPmRUBY-sonicated-0.01_1',
 'K20':'seeding_TDPmRUBY-sonicated-0.01_1',
 'L02':'seeding_TDPLCD-unsonicated-0.03_1',
 'L03':'seeding_TDPLCD-unsonicated-0.03_1',
 'L04':'seeding_TDPLCD-unsonicated-0.03_1',
 'L06':'seeding_TDPt25-unsonicated-1_1',
 'L07':'seeding_TDPt25-unsonicated-1_1',
 'L08':'seeding_TDPt25-unsonicated-1_1',
 'L10':'seeding_TDPt25-sonicated-0.01_1',
 'L11':'seeding_TDPt25-sonicated-0.01_1',
 'L12':'seeding_TDPt25-sonicated-0.01_1',
 'L14':'seeding_TDPt35-sonicated-0.3_1',
 'L15':'seeding_TDPt35-sonicated-0.3_1',
 'L16':'seeding_TDPt35-sonicated-0.3_1',
 'M02':'seeding_TDPLCD-unsonicated-0.01_1',
 'M03':'seeding_TDPLCD-unsonicated-0.01_1',
 'M04':'seeding_TDPLCD-unsonicated-0.01_1',
 'M06':'seeding_TDPt25-unsonicated-0.3_1',
 'M07':'seeding_TDPt25-unsonicated-0.3_1',
 'M08':'seeding_TDPt25-unsonicated-0.3_1',
 'M14':'seeding_TDPt35-sonicated-0.1_1',
 'M15':'seeding_TDPt35-sonicated-0.1_1',
 'N02':'control_TDPLCD-sonicated-seeds_1',
 'N03':'control_TDPLCD-sonicated-seeds_1',
 'N04':'control_TDPLCD-sonicated-seeds_1',
 'N06':'seeding_TDPt25-unsonicated-0.1_1',
 'N07':'seeding_TDPt25-unsonicated-0.1_1',
 'N08':'seeding_TDPt25-unsonicated-0.1_1',
 'N10':'control_TDPt35-unsonicated-seeds_1',
 'N11':'control_TDPt35-unsonicated-seeds_1',
 'N12':'control_TDPt35-unsonicated-seeds_1',
 'N14':'seeding_TDPt35-sonicated-0.03_1',
 'N15':'seeding_TDPt35-sonicated-0.03_1',
 'N16':'seeding_TDPt35-sonicated-0.03_1',
 'O02':'seeding_TDPLCD-sonicated-1_1',
 'O03':'seeding_TDPLCD-sonicated-1_1',
 'O04':'seeding_TDPLCD-sonicated-1_1',
 'O06':'seeding_TDPt25-unsonicated-0.03_1',
 'O07':'seeding_TDPt25-unsonicated-0.03_1',
 'O08':'seeding_TDPt25-unsonicated-0.03_1',
 'O10':'seeding_TDPt35-unsonicated-1_1',
 'O11':'seeding_TDPt35-unsonicated-1_1',
 'O12':'seeding_TDPt35-unsonicated-1_1',
 'O14':'seeding_TDPt35-sonicated-0.01_1',
 'O15':'seeding_TDPt35-sonicated-0.01_1',
 'O16':'seeding_TDPt35-sonicated-0.01_1',}


df = rename_cols(df, col_name_dict)
df = df.fillna(0)
new = add_info(df)

average = False
if average:
    averages = average_tech_replicates(new)
if not average:
    averages = new
    averages['treat'] = averages['treat2']
    averages['ID'] = averages ['treat2']

subtracted, treats_list, treats_transp = subtract_background(averages)
x=get_more_info(df=subtracted, col_to_split='ID', what_to_split='-')
# time = df[df['ID']=='time(h)']
# x = pd.concat([x, time])
x.to_csv(f'{output_folder}bg_subtracted.csv')
#-----------------------------------------------
#plotting section now

p_dict = {
     'TDPLCD':'Purples',
     'TDPmRUBY':'Reds',
     'TDPt25':'Blues',
     'TDPt35':'Greens',

}

for t, z in subtracted.groupby('t1'):
    z
    if 'blank' not in t:
        d=z
        plot_fluo_over_time(d=z, interval_mins=interval_mins, output_folder=output_folder, palette='RdBu', Exp_num=Exp_num, length_exp=20, types=['blank_subtracted'], ylim=False, savefig=True)
        # for u, dfv in z.groupby('t2'):
        #     plot_fluo_over_time(d=dfv, interval_mins=interval_mins, output_folder=output_folder, palette='RdYlGn', Exp_num=Exp_num, length_exp=20, types=['blank_subtracted'], ylim=False, savefig=False)


#then, want to calculate CHANGE over time for each treatment, and plot again.
#make sure the columns are labelled according to the treatment
alls = change_over_time_normalise(treats_transp=treats_transp)
alls.to_csv(f'{output_folder}change_in_fluo.csv')
x=get_more_info(df=alls, col_to_split='ID', what_to_split='-')




#want to just split off the 'monomer' alone (LCD alone)
monomer = x[~x['t2'].isin(['sonicated', 'unsonicated'])]
alls['replicate'][alls['t2']=='monomer'] = alls['t3'][alls['t2']=='monomer']
alls['replicate'][alls['t2']=='monmon'] = alls['t3'][alls['t2']=='monmon']

alls['t3'][alls['t2']=='monomer'] = 0
alls['t3'][alls['t2']=='monmon'] = 1

#this plots all the treatments separately (with the concentrations within those treatments staying in the same graph.)
#add different 'types' of data e.g. normalised etc. to the list 'types' if you want to plot more TYPES of data for the treatments.
alls = alls[alls['ID']!='blank']
seeding = alls[alls['t3'].isin(['0.01', '0.03', '0.1', '0.3', '1'])]
alls.to_csv(f'{output_folder}change_in_intensity_plus_info.csv')


new_plotting_melted=plot_fluo_over_time(d=alls, interval_mins=interval_mins, output_folder=output_folder, palette=p_dict, Exp_num=Exp_num, length_exp=20, types=['change_in_intensity'], ylim=False, savefig=True)

















