
'''
Script to plot vaccine impact by second- vs. third-dose coverage
'''

import numpy as np
import sciris as sc
import pandas as pd
import matplotlib.pyplot as pl
import seaborn as sns
import matplotlib.dates as mdates
import matplotlib.ticker as mtick
from matplotlib import cm


# Plot options
pad = dict(
    right  = 0.80,
    left   = 0.15,
    bottom = 0.15,
    top    = 0.95,
    hspace = 0.1,
    wspace = 0.2,
)
cax_pos = [0.82, 0.545, 0.03, 0.405]
figsize = (16,12)

bad_cmap     = 'magma_r'
good_cmap    = 'BuGn'
neutral_cmap = 'parula'
figdir = './figs'

sc.options(dpi=100)
sc.fonts(add=sc.thisdir(aspath=True) / 'avenir')
sc.options(font='Avenir')
do_show = True
do_save = True

# Create figures
figs = sc.objdict()
axes = sc.objdict()


def heatmap(data, zlabel, suptitle, filename, cmap, threshold=0.5):
    ''' Plot a matrix of results as a heat map '''
    fig,ax = pl.subplots(figsize=figsize)

    # Plot main axes
    im = ax.imshow(data, cmap=cmap, origin='lower')

    # Reconcile color limits
    cmin, cmax = im.get_clim()
    crossover = (cmax - cmin)*threshold + cmin

    # Handle labels
    for i in x:
        for j in y:

            # Label text
            k = data[j, i]
            if not np.isfinite(k):
                label = ''
            else:
                label = f'{k:0,.0f}'
            if 'Cost' in zlabel and np.isfinite(k):
                label = '$' + label
            # if 'Deaths averted' in suptitle or 'death averted' in suptitle or 'deaths averted' in suptitle:
            # if j == 0:
            #     label = 'ref'
            if 'Percent' in suptitle:
                if (j+i) > 0:
                    label += '%'
                else:
                    label = 'ref'

            # Label color
            if cmap in ['parula']:
                color = 'w' if k < crossover else 'k'
            elif cmap in ['magma_r', 'BuGn']:
                color = 'k' if (np.isnan(k) or (k < crossover)) else 'w'
            else:
                errormsg = f'Please define dark-light ordering for {cmap}'
                raise ValueError(errormsg)
            ax.text(i, j, label, ha='center', va='center', c=color)

    ax.set_xticks(x, xvals, rotation=90)
    ax.set_xlabel('Vaccine Boost Coverage', fontweight='bold')
    ax.set_yticks(y, yvals)
    ax.set_ylabel('Vaccine Prime Coverage', fontweight='bold')
    sc.boxoff(ax=ax)
    ax.axis('auto')

    # Tidying
    colorbar(fig, im, zlabel, 1)
    fig.subplots_adjust(**pad)
    sc.savefig(str(sc.path(figdir) / filename), fig=fig)

    return


def colorbar(fig, im, label, row):
    ''' Add a colorbar to the plot '''
    cbargs = dict(labelpad=15, rotation=270, fontweight='bold')
    cpos = sc.dcp(cax_pos)
    cpos[1] -= row
    cax = fig.add_axes(cpos)
    cb = fig.colorbar(im, ticklocation='right', orientation='vertical', cax=cax)
    cb.set_label(label, **cbargs)
    sc.commaticks(cax)
    return


# Load data
sc.fonts(add='C:\\Users\\jamieco\\PycharmProjects\\avenir')

# Load data
resfolder = './results'

df1 = sc.loadobj((str(sc.path(resfolder) / 'vax-cov-sweep_data_for_run_plot.obj')))
df2 = sc.loadobj((str(sc.path(resfolder) / 'vax-cov-sweep_vx_rollout_data.obj')))
df2 = pd.DataFrame.from_dict(df2)
figdir = './figs'
nab_decays = df2['nab_decay'].unique()

# Rejig data
d = sc.objdict()
o_ind_sa = 3  # Index for Omicron peak
o2_day_ind = 123  # Index for day after Omicron peak

#%% Plotting

for i, inf in enumerate(df1['new_infections_by_variant']):
    if df1['vaccine_boost'][i] == 0 and df1['vaccine_prime'][i] == 0:
        sc.options(dpi=150)
        sc.options(font='Avenir')
        fig, ax = pl.subplots()
        nab_decay = df1['nab_decay'][i]
        x = df1['datevec_sim']
        for j, var in enumerate(inf.values):
            ax.plot(x, var, label=df1['imm_source'][j])
            ax.fill_between(x, (inf.low[j,:]), (inf.high[j,:]), alpha=0.3)
        ax.set_title(f'Infections by variant, {nab_decay}')
        ax.legend()
        fig.show()

for nab_decay in nab_decays:
    df_to_use = df2[df2['nab_decay']==nab_decay]
    
    # Wrangle the data
    yvals = np.unique(df2['vaccine_boost'])
    xvals = np.unique(df2['vaccine_prime'])
    y = np.arange(len(yvals))
    x = np.arange(len(xvals))

    # Do the grouping
    df_to_use = df_to_use[df_to_use['deaths'] > 0]
    df_to_use = df_to_use[df_to_use['infections'] > 0]

    dfg = df_to_use.groupby(['vaccine_prime', 'vaccine_boost'])
    dfmean = dfg.mean().reset_index()
    dfmedian = dfg.median().reset_index()

    # First up, doses
    z_doses = dfmean.pivot('vaccine_prime', 'vaccine_boost', 'doses').reindex()
    z_doses = z_doses.reindex()
    z_doses = z_doses.values

    # Now deaths
    z_deaths = dfmean.pivot('vaccine_prime', 'vaccine_boost', 'deaths').reindex()
    z_deaths = z_deaths.reindex()
    z_deaths = z_deaths.values

    # Deaths averted
    z_deaths_averted = z_deaths[0,0] - z_deaths
    # Percent of deaths averted
    z_perc_deaths_averted = 100 * z_deaths_averted / z_deaths[0, 0]

    heatmap(
        data=z_perc_deaths_averted,
        cmap=good_cmap,
        zlabel='Percent of deaths averted',
        suptitle=f'Percent of deaths averted,, {nab_decay}',
        filename=f'vx_rollout_perc_deaths_averted_{nab_decay}.png',
    )
        
    # Now doses per deaths averted
    z_doses_deaths_averted = z_doses / z_deaths_averted
    heatmap(
        data=z_doses_deaths_averted,
        cmap=neutral_cmap,
        zlabel='Doses per death averted',
        suptitle=f'Doses per death averted, {nab_decay}',
        filename=f'vx_rollout_doses_per_death_averted_{nab_decay}.png',
    )

    # if do_show:
    #     pl.show()


print('Done')