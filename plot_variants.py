
'''

'''

import numpy as np
import sciris as sc
import pylab as pl

# Load data
sc.fonts(add='C:\\Users\\jamieco\\PycharmProjects\\avenir')

# Load data
resdir = './results'
fn = f'{str(sc.path(resdir))}/data_for_run_plot.obj'
fn2 = f'{str(sc.path(resdir))}/vx_tx_rollout_data.obj'
df = sc.load(fn)
df2 = sc.load(fn2)
figdir = './figs'

variants = df2['next_variant']

# Now plot no variant
ncolors = len(variants)
from matplotlib import cm

colors = cm.rainbow(np.linspace(0, 1, ncolors))
# Rejig data
d = sc.objdict()
o_ind_sa = 3  # Index for Omicron peak
v_ind_sa = 4  # Index for new variant peak
o2_day_ind = 123  # Index for day after Omicron peak

# %% Plotting

for i, inf in enumerate(df['new_infections_by_variant']):
    sc.options(dpi=150)
    sc.options(font='Avenir')
    fig, ax = pl.subplots()
    var = df['next_variant'][i]
    x = df['datevec_sim']
    ax.plot(x, inf.values[o_ind_sa,:], c='black', label='Omicron')
    ax.fill_between(x, (inf.low[o_ind_sa, :]), (inf.high[o_ind_sa, :]), color='black', alpha=0.3)
    ax.plot(x, inf.values[v_ind_sa,:], c='red', label='Next variant')
    ax.fill_between(x, (inf.low[v_ind_sa, :]), (inf.high[v_ind_sa, :]), color='red', alpha=0.3)
    ax.set_title(var)
    ax.legend()
    fig.show()

print('Done')

