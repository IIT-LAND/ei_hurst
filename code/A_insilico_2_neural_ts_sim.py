# A_insilico_2_neural_ts_sim.py
"""
Simulate neural time-series with a pre-specified 1/f slope.
"""

# libraries
import numpy as np
from neurodsp import sim, spectral
from neurodsp.utils import create_times, create_samples
from neurodsp.plts.spectral import plot_power_spectra
from neurodsp.plts.time_series import plot_time_series

# Set the random seed, for consistency simulating data
sim.set_random_seed(0)

# Set some general settings, to be used across all simulations
# sampling rate
fs = 500
# n_seconds will be 60 seconds * n minutes of data
n_minutes = 35
n_seconds = 60*n_minutes
times = create_times(n_seconds, fs)

# make 1/f slope range between 0 and -2, with steps of -0.1
oof_range = np.arange(start = 0,stop = 2.1,step = 0.1, dtype= float)*-1
oof_range[0] = 0

# some settings for where to save the data, filename stem, and peaks for oscillations
resultdir = "/Users/mlombardo/Dropbox/Manuscripts/AIMS_Hurst_Sex/github_repo/data/gao_model"

# loop over range of 1/f slopes
for oof in oof_range:

    # simulate a time-series with only 1/f slope and no oscillations
    components = {'sim_powerlaw' : {'exponent': oof}}

    # Simulate a combined signal with multiple oscillations
    sim_ts = sim.sim_combined(n_seconds, fs, components)

    # Save simulated time-series to a file
    np.savetxt(fname = "%s/sim_neuralts_sig_slope%0.2f_oof.txt" % (resultdir,oof), X = sim_ts)
