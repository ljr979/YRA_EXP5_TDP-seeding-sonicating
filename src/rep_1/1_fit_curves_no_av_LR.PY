import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy.optimize import curve_fit
from loguru import logger
import warnings
warnings.filterwarnings("ignore")
logger.info('Import ok')
def r_squared_calculator(x_vals, y_vals, function, popt):
    residuals = y_vals - np.array([function(x, *popt) for x in x_vals])
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_vals-np.mean(y_vals))**2)
    return 1 - (ss_res / ss_tot)
# Define the fitting functions
def sigmoid(x, bottom, top, X50):
    return bottom + ((top - bottom) / (1 + np.exp((X50 - x))))

def boltzmann(x, bottom, top, V50, slope):
    return bottom + ((top - bottom) / (1 + np.exp((V50 - x)/slope)))

def boltzmann_with_exp_decay(x, bottom, top, V50, slope, decay):
    return bottom + ((top - bottom) / (1 + np.exp((V50 - x) / slope))) - (x**decay)

def boltzmann_with_linear(x, bottom, top, V50, slope, m):
    return bottom + ((top - bottom) / (1 + np.exp((V50 - x) / slope))) - (m*x)

def boltzmann_with_both(x, bottom, top, V50, slope, decay, m):
    return bottom + ((top - bottom) / (1 + np.exp((V50 - x) / slope))) - (x**decay)-(m*x)


if __name__ == '__main__':

    filter_cols = []

    input_path = 'results/rep_1/0_processing/change_in_intensity_plus_info.csv'
    output_folder = 'results/rep_1/1_fitting_curves/'

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Read in change in fluorescence intensity data
    clusters_summary = pd.read_csv(input_path)
    change_in = clusters_summary[clusters_summary['type']=='change_in_intensity']
    info_cols = ['ID','type','t1','t2','t3', 'replicate']
    change_in.drop([col for col in change_in.columns.tolist() if 'Unnamed: 0' in col], axis=1, inplace=True)
    info_bits = change_in[info_cols]


    # Pivot data to long form for easier manipulation
    clean_data = pd.melt(
        change_in,
        id_vars=info_cols,
        value_name='intensity',
        var_name='timepoint'
    )
    clean_data['timepoint'] = (clean_data['timepoint'].astype(int))*4/60
    clean_data.to_csv(f'{output_folder}clean_for_fitting.csv')

    # Test each fit version to compare
    fit_params = []
    no_fit = []
    for (sample, t2, t3), df in clean_data.groupby(['t1', 't2', 't3']):
        sample
        t2
        t3
        df
        for replicate, df1 in df.groupby('replicate'):
            replicate
            inform = f'{sample}-{t2}-{t3}-{replicate}'
            df1

            x_vals = df1['timepoint'].values
            y_vals = df1['intensity'].values
            if 'seeds' not in t3:
                 try:
                    # Test with original 'plain' boltzmann
                    popt, pcov = curve_fit(boltzmann, x_vals, y_vals)
                    
                    x_fit = np.linspace(np.min(x_vals), np.max(x_vals), 1000)
                    y_fit = np.array([boltzmann(x, *popt) for x in x_fit])
                    
                    y_vals_fit = np.array([boltzmann(x, *popt) for x in x_vals])
                    y_subtract = y_vals-y_vals_fit

                    rsq = r_squared_calculator(x_vals, y_vals, boltzmann, popt)
                    
                    fig, axes = plt.subplots(1, 4, figsize=(18, 5))
                    
                    sns.lineplot(x=x_vals, y=y_subtract, linestyle='--', ax=axes[0], color='black')
                    sns.lineplot(x=x_vals, y=y_vals, linestyle='--', ax=axes[0], color='#086788')
                    sns.lineplot(x=x_fit, y=y_fit, ax=axes[0], color='#086788')
                    axes[0].annotate(f'R$^2$={round(rsq, 3)}', (0.8*np.max(x_vals), 0.25*np.max(y_vals)))
                    axes[0].annotate(f'boltzmann', (0.8*np.max(x_vals), 0.8*np.max(y_vals)))
                    axes[0].set_title(inform)
                    bottom, top, V50, slope = popt
                    tlag = ((4*slope*bottom) / (top - bottom)) + V50 - 2*slope
                    elongation = (top-bottom)/(4*slope)
                    fit_type='boltzmann'
                    decay='none'
                    fit_params.append([inform, bottom, top, V50, slope, decay, tlag, round(rsq,3), fit_type, elongation])

          
                    #now next type of fit
                    popt, pcov = curve_fit(boltzmann_with_exp_decay, x_vals, y_vals)
                    
                    x_fit = np.linspace(np.min(x_vals), np.max(x_vals), 1000)
                    y_fit = np.array([boltzmann_with_exp_decay(x, *popt) for x in x_fit])
                    
                    rsq = r_squared_calculator(x_vals, y_vals, boltzmann_with_exp_decay, popt)
                    

                    sns.lineplot(x=x_vals, y=y_vals, linestyle='--', ax=axes[1], color='#086788')
                    sns.lineplot(x=x_fit, y=y_fit, ax=axes[1], color='#086788')
                    axes[1].annotate(f'R$^2$={round(rsq, 3)}', (0.8*np.max(x_vals), 0.25*np.max(y_vals)))
                    axes[1].annotate(f'boltzmann_with_exp_decay', (0.8*np.max(x_vals), 0.8*np.max(y_vals)))
                    axes[1].set_title(inform)
                    
                    bottom, top, V50, slope, decay = popt
                    tlag = ((4*slope*bottom) / (top - bottom)) + V50 - 2*slope
                    elongation = (top-bottom)/(4*slope)
                    fit_type='boltzmann_with_exp_decay'

                    fit_params.append([inform, bottom, top, V50, slope, decay, tlag, round(rsq,3), fit_type, elongation])



                    #third type of fit
                    # Test with modified boltzmann
                    popt, pcov = curve_fit(boltzmann_with_linear, x_vals, y_vals)
                    
                    x_fit = np.linspace(np.min(x_vals), np.max(x_vals), 1000)
                    y_fit = np.array([boltzmann_with_linear(x, *popt) for x in x_fit])
                    
                    rsq = r_squared_calculator(x_vals, y_vals, boltzmann_with_linear, popt)

                    
                    sns.lineplot(x=x_vals, y=y_vals, linestyle='--', ax=axes[2], color='#E4572E')
                    sns.lineplot(x=x_fit, y=y_fit, ax=axes[2], color='#E4572E')

                    axes[2].annotate(f'R$^2$={round(rsq, 3)}', (0.8*np.max(x_vals), 0.25*np.max(y_vals)))
                    axes[2].annotate(f'boltzmann_with_linear', (0.8*np.max(x_vals), 0.8*np.max(y_vals)))

                    
                    bottom, top, V50, slope, m = popt
                    tlag = ((4*slope*bottom) / (top - bottom)) + V50 - 2*slope
                    elongation = (top-bottom)/(4*slope)
                    fit_type='boltzmann_with_linear'
                    fit_params.append([inform, bottom, top, V50, slope, m, tlag, round(rsq,3), fit_type, elongation])
                
                    #last fit type
                    popt, pcov = curve_fit(boltzmann_with_both, x_vals, y_vals)
                    
                    x_fit = np.linspace(np.min(x_vals), np.max(x_vals), 1000)
                    y_fit = np.array([boltzmann_with_both(x, *popt) for x in x_fit])
                    
                    rsq = r_squared_calculator(x_vals, y_vals, boltzmann_with_both, popt)
                    

                    sns.lineplot(x=x_vals, y=y_vals, linestyle='--', ax=axes[3], color='#086788')
                    sns.lineplot(x=x_fit, y=y_fit, ax=axes[3], color='#086788')
                    axes[3].annotate(f'R$^2$={round(rsq, 3)}', (0.8*np.max(x_vals), 0.25*np.max(y_vals)))
                    axes[3].annotate(f'boltzmann_with_both', (0.8*np.max(x_vals), 0.8*np.max(y_vals)))
                    axes[3].set_title(inform)
                    plt.show()

                    bottom, top, V50, slope, decay, m = popt
                    tlag = ((4*slope*bottom) / (top - bottom)) + V50 - 2*slope
                    elongation = (top-bottom)/(4*slope)

                    fit_type = 'boltzmann_with_both'
                    decay = f'{decay}_{m}'
                    fit_params.append([inform, bottom, top, V50, slope, decay, tlag, round(rsq,3), fit_type, elongation])
                 except:
                     print(f'couln"t fit {inform}')
                     continue
                
                # except Exception:
                #     logger.info(f'No fit found for {inform}-')
                #     no_fit.append(inform)
 
fit_params = pd.DataFrame(fit_params, columns=['sample', 'bottom', 'top', 'V50', 'slope', 'decay', 'tlag','rquare', 'fit_type', 'elongation_rate'])


# nots = pd.DataFrame(no_fit, columns=['sample']).to_csv(f'{output_folder}no_fit.csv')
# fit_params.to_csv(f'{output_folder}fit_parameters.csv')


#checking the fitting to decide on one
#rquare_easy = fit_params[['sample','rquare','fit_type']]
melted_params = pd.melt(fit_params, id_vars=['sample',
 'fit_type','bottom',
 'top',
 'V50',
 'slope',
 'decay',
 'tlag', 'elongation_rate'
 ], value_vars=['rquare'])

sns.barplot(melted_params, x="fit_type", y="value")


#want to have a little section to look at what the most efficient fitting is and 

for fit_type, df in melted_params.groupby('fit_type'):
    df
    fit_type
    amazing=df[df['value']>0.99]
    not_amazing=df[df['value']<0.99]

    prop_good=len(amazing)/len(df)
    print(fit_type, prop_good)



#ok so now loop through each sample and add it's best fit to a dictionary with the rsquare val
#This finds the BEST fit for every sample.
good_fits = []
for sample, df in melted_params.groupby('sample'):
    sample
    df
    maximum=max(df['value'])
    best_fit=df[df['value']==maximum]
    if len(best_fit) > 1:
        best_fit=best_fit.iloc[:1]
    good_fits.append(best_fit)

good_fits=pd.concat(good_fits)
len(good_fits)
#then it subsets to show you which ones are not so good still
not_so_good_fits=good_fits[good_fits['value']<0.98]
best_filtered = good_fits[good_fits['value']>0.98]

perc = len(not_so_good_fits)/len(good_fits)*100
print(f'{perc} are not good enough fits to keep')
best_filtered.rename(columns={'value':'rquare'}, inplace=True)
best_filtered.drop([col for col in best_filtered.columns.tolist() if 'variable' in col], axis=1, inplace=True)
best_filtered.to_csv(f'{output_folder}fit_parameters.csv')
not_so_good_fits.to_csv(f'{output_folder}no_fit.csv')