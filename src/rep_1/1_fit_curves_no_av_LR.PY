import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from scipy.optimize import curve_fit
from loguru import logger

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
                    
                    rsq = r_squared_calculator(x_vals, y_vals, boltzmann, popt)
                    
                    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
                    sns.lineplot(x=x_vals, y=y_vals, linestyle='--', ax=axes[0], color='#086788')
                    sns.lineplot(x=x_fit, y=y_fit, ax=axes[0], color='#086788')
                    axes[0].annotate(f'R$^2$={round(rsq, 3)}', (0.8*np.max(x_vals), 0.25*np.max(y_vals)))
                    axes[0].annotate(f'boltzmann', (0.8*np.max(x_vals), 0.8*np.max(y_vals)))
                    axes[0].set_title(inform)
                    
                    # Test with modified boltzmann
                    popt, pcov = curve_fit(boltzmann_with_exp_decay, x_vals, y_vals)
                    
                    x_fit = np.linspace(np.min(x_vals), np.max(x_vals), 1000)
                    y_fit = np.array([boltzmann_with_exp_decay(x, *popt) for x in x_fit])
                    
                    rsq = r_squared_calculator(x_vals, y_vals, boltzmann_with_exp_decay, popt)
                    
                    sns.lineplot(x=x_vals, y=y_vals, linestyle='--', ax=axes[1], color='#E4572E')
                    sns.lineplot(x=x_fit, y=y_fit, ax=axes[1], color='#E4572E')

                    axes[1].annotate(f'R$^2$={round(rsq, 3)}', (0.8*np.max(x_vals), 0.25*np.max(y_vals)))
                    axes[1].annotate(f'new_fit', (0.8*np.max(x_vals), 0.8*np.max(y_vals)))
                    plt.show()
                    
                
                    bottom, top, V50, slope, decay = popt
                    tlag = ((4*slope*bottom) / (top - bottom)) + V50 - 2*slope

                    fit_params.append([inform, bottom, top, V50, slope, decay, tlag, round(rsq,3)])
                except Exception:
                    logger.info(f'No fit found for {inform}')
                    no_fit.append(inform)
 
fit_params = pd.DataFrame(fit_params, columns=['sample', 'bottom', 'top', 'V50', 'slope', 'decay', 'tlag','rquare'])


nots = pd.DataFrame(no_fit, columns=['sample']).to_csv(f'{output_folder}no_fit.csv')
fit_params.to_csv(f'{output_folder}fit_parameters.csv')