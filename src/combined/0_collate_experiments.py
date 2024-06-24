import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

input_folder = 'results/'
output_folder = 'results/combined/'

if not os.path.exists(output_folder):
   os.makedirs(output_folder)

#first want to collate the fit parameters from all three experiments
fit_params_filepaths = [[f'{root}/{name}' for name in files if 'fit_parameters.csv' in name]for root, dirs, files in os.walk(f'{input_folder}/')]
#flatten list
fit_params_filepaths = [item for sublist in fit_params_filepaths for item in sublist if 'exp_7' not in item]

collated = []
for files in fit_params_filepaths:
   files = files.replace('\\', '/')
   exp_num = files.split('//')[1].split('/')[0]
   params = pd.read_csv(f'{files}')
   params['exp_num'] = exp_num
   collated.append(params)

collated = pd.concat(collated)


collated.to_csv(f'{output_folder}collated_fit_params.csv')