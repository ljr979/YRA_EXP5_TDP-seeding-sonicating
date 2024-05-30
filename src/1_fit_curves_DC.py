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

    input_path = f'raw_data/troubleshooting_boltzmann-fit.csv'
    output_folder = f'results/fitting_curves/'

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Read in change in fluorescence intensity data
    raw_data = pd.read_csv(input_path, header=None)
    raw_data.columns = ['_', 'sample', 'data-type',] + list(range(raw_data.shape[1] - 3))
    raw_data.drop('_', axis=1, inplace=True)
    
    # Pivot data to long form for easier manipulation
    clean_data = pd.melt(
        raw_data,
        id_vars=['sample', 'data-type'],
        value_name='intensity',
        value_vars=list(range(raw_data.shape[1] - 2)),
        var_name='timepoint'
    )
    
    
    # Test each fit version to compare
    fit_params = []
    
    for sample, df in clean_data[clean_data['data-type'] == 'original'].groupby(['sample']):
             
        x_vals = df['timepoint'].values
        y_vals = df['intensity'].values
    
        # Test with original 'plain' boltzmann
        popt, pcov = curve_fit(boltzmann, x_vals, y_vals)
        
        x_fit = np.linspace(np.min(x_vals), np.max(x_vals), 1000)
        y_fit = np.array([boltzmann(x, *popt) for x in x_fit])
        
        rsq = r_squared_calculator(x_vals, y_vals, boltzmann, popt)
        
        fig, axes = plt.subplots(1, 2, figsize=(11, 5))
        sns.lineplot(x=x_vals, y=y_vals, linestyle='--', ax=axes[0], color='#086788')
        sns.lineplot(x=x_fit, y=y_fit, ax=axes[0], color='#086788')
        axes[0].annotate(f'R$^2$={round(rsq, 3)}', (0.8*np.max(x_vals), 0.25*np.max(y_vals)))
        axes[0].set_title(sample)
        
        # Test with modified boltzmann
        popt, pcov = curve_fit(boltzmann_with_exp_decay, x_vals, y_vals)
        
        x_fit = np.linspace(np.min(x_vals), np.max(x_vals), 1000)
        y_fit = np.array([boltzmann_with_exp_decay(x, *popt) for x in x_fit])
        y_bfit = np.array([boltzmann(x, *popt[:-1]) for x in x_fit])
        
        rsq = r_squared_calculator(x_vals, y_vals, boltzmann_with_exp_decay, popt)
        
        sns.lineplot(x=x_vals, y=y_vals, linestyle='--', ax=axes[1], color='#E4572E')
        sns.lineplot(x=x_fit, y=y_fit, ax=axes[1], color='#E4572E')
        sns.lineplot(x=x_fit, y=y_bfit, ax=axes[1], color='black')
        axes[1].annotate(f'R$^2$={round(rsq, 3)}', (0.8*np.max(x_vals), 0.25*np.max(y_vals)))
        plt.show()
        
    
        bottom, top, V50, slope, decay = popt
        tlag = ((4*slope*bottom) / (top - bottom)) + V50 - 2*slope

        fit_params.append([sample, bottom, top, V50, slope, decay, tlag])

fit_params = pd.DataFrame(fit_params, columns=['sample', 'bottom', 'top', 'V50', 'slope', 'decay', 'tlag'])


