# %%
import matplotlib.pyplot as plt
from Classes.struct2dict import *  # my own function that lets us easily convert matlab to python
from Classes import analyzer_module
import numpy
# %% Loading data into python
data = struct2dict(
    'D:/Code/NeurotarCode/SampleData/PlaceCells/TSeries-07242019-1414-001_registered_data.mat')
floating = struct2dict(
    'D:/Code/NeurotarCode/SampleData/PlaceCells/floating_data_WTR010_07_24_session1.mat')

# %% Loading real data to try...
# %% Initialize the analyzer
a = analyzer_module.Analyzer(data, floating)
# %% Let's figure out which samples are when the mouse is actually moving times
a.resample_neurotar_data(avg_window=0)

a.find_moving_samples()

# %% Binning the data
a.calculate_head_direction(n_sections=4)  # Important: There is a very small discrepancy between python and MATLAB when
# they calculate correlation coefficient. It might be due to higher precision numbers in python?

is_head_direction = a.mean_quadrant_correlation > 0.2

print('Num tuned cells: {}'.format(numpy.sum(is_head_direction)))
[head_idx] = numpy.where(is_head_direction)
print('Cell IDs: {}'.format(head_idx))