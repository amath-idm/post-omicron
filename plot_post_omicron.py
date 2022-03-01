'''
Plot post omicron runs
'''

import numpy as np
import sciris as sc
import pylab as pl
import matplotlib.dates as mdates
import matplotlib.ticker as mtick

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
    'Emerged from Omicron, more severe', 'Emerged from WT, more severe',
    'New cluster, more severe'
]

variant_timing = ['2022-02-25', '2022-04-25', '2022-08-25']

fig, axes = pl.subplots(3, 3, figsize=(10, 10), sharex=True, sharey=False)

x = df['datevec_sim'][o2_day_ind:]

for i, inf in enumerate(df['new_infections_by_variant']):
    if df['vaccine'][i] == 'Status quo' and df['next_variant'][i] in new_variants:
        var = df['next_variant'][i]
        var_ax = new_variants.index(var)
        var_date = df['new_variant_day'][i]
        date_ax = variant_timing.index(var_date)
        ax = axes[var_ax, date_ax]

        factor_inf = 100 / inf.values[o_ind, 0:o2_day_ind].max()
        factor_sev = 100 / df['new_severe_by_variant'][i].values[o_ind, 0:o2_day_ind].max()

        ax.plot(x, inf.values[o_ind, o2_day_ind:]*factor_inf, c=colors[0], label='Omicron')
        ax.plot(x, inf.values[v_ind, o2_day_ind:]*factor_inf, c=colors[var_ax+1], label=var)

        ax.plot(x, df['new_severe_by_variant'][i].values[o_ind, o2_day_ind:]*factor_sev, c=colors[0], linestyle=':')
        ax.plot(x, df['new_severe_by_variant'][i].values[v_ind, o2_day_ind:]*factor_sev, c=colors[var_ax+1], linestyle=':')
        ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))

lgd = axes[2,1].legend(loc=9, bbox_to_anchor=(0.5,-0.2))
axes[0,0].set_title(variant_timing[0])
axes[0,1].set_title(variant_timing[1])
axes[0,2].set_title(variant_timing[2])
axes[0,0].set_ylabel('% of Omicron peak')
axes[1,0].set_ylabel('% of Omicron peak')
axes[2,0].set_ylabel('% of Omicron peak')
text1 = axes[0,2].text(1.05, 0.5 ,new_variants[0], transform=axes[0,2].transAxes)
text2 = axes[1,2].text(1.05, 0.5, new_variants[1], transform=axes[1,2].transAxes)
text3 = axes[2,2].text(1.05, 0.5, new_variants[2], transform=axes[2,2].transAxes)
fig.subplots_adjust(bottom=0.15)
fig.subplots_adjust(right=0.78)
fig.savefig(str(sc.path(figdir) / 'post_omicron.png'), bbox_inches='tight')
print('Done.')

