import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger
from lmfit import Model, Parameters, minimize, report_fit
from scipy.optimize import curve_fit
from scipy import stats
from loguru import logger

logger.info('Import ok')
def r_squared_calculator(x_vals, y_vals, function, popt):
    residuals = y_vals - function(x_vals, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((y_vals-np.mean(y_vals))**2)
    return 1 - (ss_res / ss_tot)

# Define the fitting functions
def sigmoid(x, bottom, top, X50):
    return bottom + ((top - bottom) / (1 + np.exp((X50 - x))))


def boltzmann(x, bottom, top, V50, slope):
    return bottom + ((top - bottom) / (1 + np.exp((V50 - x)/slope)))


def denaturant(urea, top, bottom, cM, m):
    # adapted from https://en.wikipedia.org/wiki/Equilibrium_unfolding, by keeping terms for bottom as in Savitski, subbing deltaG into standard equation, and reintroducing bottom term as per boltzmann

    temp_constant = 298.15
    gas_constant = 8.31446261815324
    constant = temp_constant * gas_constant

    y = bottom + ((top - bottom) / (1 + np.exp((m*(cM-urea)/constant))))

    # deltaG can then be calculated as m(cM-urea) - generally calculated at 0M urea therefore m(cM)

    return y


def denaturant_fit(compiled, info_cols, quant_cols):
    """Attempts to fit a sigmoid to each row. Returns sample_params dict where keys are sequences"""

    fit_params = {}
    no_fit = []
    proper_fit = []
    for info, quant_data in compiled.set_index(info_cols).iterrows():
        info
        logger.info(f'working on {info}')
        quant_data
        # extract x and y vals for fitting
        y_vals = np.array(list(quant_data[quant_cols]))
        x_vals = np.array([float(x) for x in quant_cols])

        # Attempt fitting
        try:
            model = Model(denaturant)
            params = model.make_params(
                bottom=-1, top=1, cM=3, m=-10000)
            result = model.fit(y_vals, params, urea=x_vals)
            r_squared = r_squared_calculator(
                x_vals, y_vals, denaturant, result.values.values())

            # Collect fitted parameters
            fit_stats = pd.DataFrame()
            for parameter, details in result.params.items():
                fit_stats[f'{parameter}_value'] = [details.value]
                fit_stats[f'{parameter}_stderr'] = [details.stderr]
                fit_stats[f'{parameter}_relerr'] = fit_stats[f'{parameter}_stderr'].values[0] / \
                    fit_stats[f'{parameter}_value'].values[0] * 100

            # add r-squared value, key info
            fit_stats['r_squared'] = r_squared
            fit_stats['key'] = [info]

            fit_params[info] = fit_stats
            proper_fit.append(info[0])

        except:
            logger.info(f'No fit found for {info}')
            no_fit.append(info[0])
    return fit_params, no_fit, proper_fit


def sigmoid_filter(summary, filter_R2=True, filter_range=True, filter_cM=True, filter_relerr=True, filter_direction=True):
    
    # apply filtering criteria
    filtered = summary.copy()

    if filter_R2:
        # Remove R2 < filter threshold
        filtered['filter_R2'] = [1 if R2 > 0.75 else 0 for R2 in filtered['r_squared']]
        logger.info(f"R2 filter: {filtered['filter_R2'].sum()}")

    if filter_range:
        # Remove top/bottom outside range - threshold = 10?
        filtered = filtered[(abs(filtered['top_value']) < 10) & (abs(filtered['bottom_value']) < 10)]
        filtered['filter_range'] = [1 if (abs(val_1) < 10) & (abs(val_2) < 10) else 0 for val_1, val_2 in filtered[['top_value', 'bottom_value']].values]
        logger.info(f"Range filter: {filtered['filter_range'].sum()}")

    if filter_cM:
        # Remove cM outside range tested
        filtered['filter_cM'] = [1 if (val < 6) & (val > 0) else 0 for val in filtered['cM_value']]
        logger.info(f"cM filter: {filtered['filter_cM'].sum()}")

    if filter_relerr:
        # Remove fits with > 50% uncertainty in cM fit
        filtered['filter_relerr'] = [1 if val < 50 else 0 for val in filtered['cM_relerr']]
        logger.info(f"Relative cM error: {filtered['filter_relerr'].sum()}")

    if filter_direction:
        # Remove sigmoids that trend upward
        filtered['filter_direction'] = [1 if val_0 > val_6 else 0 for val_0, val_6 in zip(filtered['0M_value'], filtered['0.2M_value'])]
        logger.info(f"Sigmoid direction: {filtered['filter_direction'].sum()}")

    filter_cols = [col for col in filtered.columns.tolist() if 'filter_' in str(col)]
    filtered['filter_count'] = filtered[filter_cols].sum(axis=1)
    filtered['filter_all'] = [1 if num == len(filter_cols) else 0 for num in filtered['filter_count']]
    logger.info(f"All filters: {filtered['filter_all'].sum()}")

    # add filtering info to original df
    summary['filtered'] = filtered['filter_all']

    return summary, filtered

if __name__ == '__main__':

    filter_cols = []

    input_path = f'results/rep_1/0_processing/change_in_fluo.csv'
    output_folder = f'results/rep_1/1_fitting_curves/'

    if not os.path.exists(output_folder):
        os.mkdir(output_folder)

    # Read in cluster data
    clusters_summary = pd.read_csv(input_path)
    change_in = clusters_summary[clusters_summary['type']=='change_in_intensity']
    info_cols = ['ID','type','t1','t2','t3']
    change_in.drop([col for col in change_in.columns.tolist() if 'Unnamed: 0' in col], axis=1, inplace=True)
    info_bits = change_in[info_cols]

#------------------------
    quant_cols = [col for col in change_in.columns.tolist() if col not in info_cols]
    clusters = change_in[info_cols+quant_cols]
#__________________________________HERE IS WHERE I LOSE SEVERAL TREATMENTS- HERE NEED TO COMPARE BETWEEN AND THEN FIND THE ONES THAT ARE MISSING AND WHY? PLOT THEM? I SUCK AT THIS
    # complete denaturant fit and return a list of those that couldn't be fit
    fit_params, no_fit, proper_fit = denaturant_fit(clusters, info_cols=info_cols, quant_cols=quant_cols)
    fitting_parameters = pd.concat(fit_params.values()).reset_index(drop=True)
#_____________________
    # add back useful info
    fitting_parameters[info_cols] = pd.DataFrame(fitting_parameters['key'].tolist(), index=fitting_parameters.index)
    summary = pd.merge(clusters, fitting_parameters, on=info_cols, how='inner')
    summary.to_csv(f'{output_folder}fitting_parameters.csv')
#-----------------------
    # generate "fitted" results
    sigmoid_fitted_vals = {}
    for cluster, df in summary.iterrows():
        # generate fitted values
        (bottom, top, cM, m, r_squared, cluster, protein, sequence) = tuple(df[['bottom_value', 'top_value', 'cM_value', 'm_value', 'r_squared', 'ID', 't1', 't2']])
        y_vals = denaturant(np.array([float(x) for x in quant_cols]), top, bottom, cM, m)
        sigmoid_fitted_vals[cluster] = y_vals
    sigmoid_fitted_vals = pd.DataFrame(sigmoid_fitted_vals).T.reset_index()
    sigmoid_fitted_vals.columns = [['ID']+ quant_cols]
    sigmoid_fitted_vals.to_csv(f'{output_folder}fits.csv')
    #want to include also the data that has the actual data to be fit. so find those that got fit properly
    matched = clusters[clusters['ID'].isin(proper_fit)]
    #now get just time timeseries data
    timeseries = matched[quant_cols]
    #now add the ID
    timeseries['ID'] = matched['ID']
    #now want to concatinate these two dataframes and add a column to each which says whether it is the fitted data or the original data
    timeseries ['type'] = 'original'
    test_Df = sigmoid_fitted_vals.copy()
    test_Df['type'] = 'fit'
    # insert column using insert(position,column_name, 
    # first_column) function 
    test_Df.insert(1, 'type', test_Df.pop('type')) 
    timeseries.insert(0, 'ID', timeseries.pop('ID')) 
    timeseries.insert(1, 'type', timeseries.pop('type')) 
    #make sure column names are identical before concatinating them
    test_Df.columns=timeseries.columns.tolist()
    collated = pd.concat([test_Df, timeseries], axis=0)
    collated.to_csv(f'{output_folder}combined_fit_and_original.csv')
    #when i have time, need to plot this as GROUPED by the type of protein, ad the sonicated or not sonicated or monomer. Then plot those at all concentrations on one graph.
    #also need to make sub-output folder for each of these groupings so annoying too many graphs
    for info, df in collated.groupby(['ID']):  
        ID=info[0]
        melt = pd.melt(df, id_vars= ['ID', 'type'], var_name='timepoint', value_name='intensity' )
        melt['timepoint'] = melt['timepoint'].astype(int)
        Fig, ax = plt.subplots()
        sns.lineplot(melt, x='timepoint', y='intensity', hue='type', ax=ax)
        plt.ylim(-420,60000)
        ax.set_xlabel('Timepoint')
        ax.set_ylabel('ThT intensity (A.U.)')
        plt.title(f'{ID}')
        plt.savefig(f'{output_folder}{ID}_curve_fit.png')
        plt.show()





