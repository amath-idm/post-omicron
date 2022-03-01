'''
Plot post omicron runs
'''

import numpy as np
import sciris as sc
import pylab as pl

# Load data
sc.fonts(add='C:\\Users\\jamieco\\PycharmProjects\\avenir')

# Load data
resdir = './results'
fn = f'{str(sc.path(resdir))}/data_for_run_plot.obj'
df = sc.load(fn)

figdir = './figs'

ncolors = 8
from matplotlib import cm
colors = cm.rainbow(np.linspace(0, 1, ncolors))
# Rejig data
d = sc.objdict()
d_inf = 2 # Index for Delta
o_ind = 3  # Index for Omicron
v_ind = 4  # Index for next variant
o2_day_ind = 123  # Index for day after Omicron peak

sc.options(dpi=150)
sc.options(font='Avenir')

new_variants = [
    'Emerged from Omicron', 'Emerged from WT', 'New cluster',
    'Emerged from Omicron, more severe', 'Emerged from WT, more severe',
    'New cluster, more severe'
]

variant_timing = ['2022-02-25', '2022-04-25', '2022-08-25']

fig, axes = pl.subplots(6, 3, figsize=(8, 10), sharex=True)

x = df['datevec_sim']

for i, inf in enumerate(df['new_infections_by_variant']):
    if df['vaccine'][i] == 'Status quo' and df['next_variant'] != 'None':
        var = df['next_variant'][i]
        var_ax = new_variants.index(var)
        var_date = df['new_variant_day'][i]
        date_ax = variant_timing.index(var_date)
        ax = axes[var_ax, date_ax]

        factor_inf = 100 / inf.values[o_ind, 0:o2_day_ind].max()
        factor_sev = 100 / df['new_severe_by_variant'][i].values[o_ind, 0:o2_day_ind].max()

        ax.plot(x, inf.values[o_ind, :]*factor_inf, c=colors[0], label='Omicron')
        ax.plot(x, inf.values[v_ind, :]*factor_inf, c=colors[1], label=var)

        ax.plot(x, df['new_severe_by_variant'][i].values[o_ind, :]*factor_sev, c=colors[0], linestyle=':')
        ax.plot(x, df['new_severe_by_variant'][i].values[v_ind, :]*factor_sev, c=colors[1], linestyle=':')

axes[0,0].set_title(variant_timing[0])
axes[0,1].set_title(variant_timing[1])
axes[0,2].set_title(variant_timing[2])
