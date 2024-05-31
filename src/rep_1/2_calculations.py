import os, re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from loguru import logger

input_folder = 'results/rep_1/1_fitting_curves/'
output_folder = 'results/rep_1/2_calculations/'


if not os.path.exists(output_folder):
   os.makedirs(output_folder)

to_read = 