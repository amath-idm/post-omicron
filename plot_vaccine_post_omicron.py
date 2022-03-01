
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
    right  = 0.88,
    left   = 0.15,
    bottom = 0.17,
    top    = 0.95,
    hspace = 0.1,
    wspace = 0.2,
)
cax_pos = [0.9, 0.705, 0.03, 0.25]
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

# Create the fig
fig,axs = pl.subplots(nrows=3, ncols=3, figsize=figsize)


def heatmap(data, zlabel, suptitle, filename, cmap, threshold=0.5, time_lab=None, row=None):
    ''' Plot a matrix of results as a heat map '''
    time_col_map = {'2022-02-25': 0, '2022-04-25': 1, '2022-08-25': 2}
    col = time_col_map[time_lab]
    print(f'Plotting: time={time_lab}, row={row}, col={col}')
    ax = axs[row, col]

    # Remove non-significant data
    if zlabel != 'Cumulative deaths':
        for i in x:
            for j in y:
                if z_deaths_ttest[j, i] > 0.05:
                    data[j, i] = np.nan  # To turn off these cells

        data[0,:] = np.nan

    # Plot main axes
    im = ax.imshow(data, cmap=cmap, origin='lower')

    # Reconcile color limits
    cmin, cmax = im.get_clim()
    if col == 2:
        other_im = axs[row, 0].get_images()[0]
        other_cmin, other_cmax = other_im.get_clim()
        cmax = max(cmax, other_cmax)
        im.set_clim((0, cmax))
        other_im.set_clim(0, cmax)
    crossover = (cmax - cmin) * threshold + cmin

    # Handle labels
    for i in x:
        for j in y:

            # Label text
            k = data[j, i]
            if zlabel == 'Cumulative deaths':
                label = f'{k:0,.0f}'
            else:
                if j == 0:
                    label = 'ref'
                else:
                    if not np.isfinite(k):
                        label = 'NS'
                    else:
                        label = f'{k:0,.0f}'

                        if 'Percent' in suptitle:
                            if j > 0:
                                label += '%'

            # Label color
            if cmap in ['parula']:
                color = 'w' if k < crossover else 'k'
            elif cmap in ['magma_r', 'BuGn']:
                color = 'k' if (np.isnan(k) or (k < crossover)) else 'w'
            else:
                errormsg = f'Please define dark-light ordering for {cmap}'
                raise ValueError(errormsg)
            ax.text(i, j, label, ha='center', va='center', c=color)

    if row == 2:
        ax.set_xticks(x, xlabels, rotation=90)
        ax.set_xlabel('New variant', fontweight='bold')
    else:
        ax.set_xticks([])
    if row == 0:
        ax.set_title(f'New variant introduced on\n{time_lab}')
    if col == 0:
        ax.set_yticks(y, ylabels)
        ax.set_ylabel('Vaccine strategy', fontweight='bold')
    else:
        ax.set_yticks([])
    sc.boxoff(ax=ax)
    ax.axis('auto')

    # Tidying
    if col == 2:
        colorbar(fig, im, zlabel, row)
    fig.subplots_adjust(**pad)
    if do_save and row == 2 and col == 2:
        sc.savefig(str(sc.path(figdir) / filename), fig=fig)

    return


def colorbar(fig, im, label, row):
    ''' Add a colorbar to the plot '''
    cbargs = dict(labelpad=15, rotation=270, fontweight='bold')
    cpos = sc.dcp(cax_pos)
    cpos[1] -= row*0.268
    cax = fig.add_axes(cpos)
    cb = fig.colorbar(im, ticklocation='right', orientation='vertical', cax=cax)
    cb.set_label(label, **cbargs)
    sc.commaticks(cax)
    return

vaccines = {
    'Status quo': {
        'label': 'Status quo'
    },
    'Boost after 3m (WT)': {
        'label': 'WT-based boost \nafter 3 months'
    },
    'Boost after 6m (WT)': {
        'label': 'WT-based boost \nafter 6 months'
    },
    'Omicron 1-dose after 3m': {
        'label': 'Omicron-based 1-dose \nafter 3 months'
    },
    'Omicron 2-dose after 3m': {
        'label': 'Omicron-based 2-dose \nafter 3 months'
    },
    'Omicron 1-dose after 6m': {
        'label': 'Omicron-based 1-dose \nafter 6 months'
    },
    'Omicron 2-dose after 6m': {
        'label': 'Omicron-based 2-dose \nafter 6 months'
    },
}

variants = ['None', 'Emerged from Omicron', 'Emerged from WT', 'New cluster',
    'Emerged from Omicron, more severe', 'Emerged from WT, more severe',
    'New cluster, more severe'
]

variant_labels = ['No new variant', 'Emerged from Omicron', 'Emerged from WT', 'New cluster',
    'Emerged from Omicron, \nmore severe', 'Emerged from WT, \nmore severe',
    'New cluster, \nmore severe'
]

new_variant_days = ['2022-02-25', '2022-04-25', '2022-08-25']

# Load data
sc.fonts(add='C:\\Users\\jamieco\\PycharmProjects\\avenir')

# Load data
resfolder = './results'
df = sc.loadobj((str(sc.path(resfolder) / 'vx_tx_rollout_data.obj')))
df = pd.DataFrame.from_dict(df)
figdir = './figs'
nab_decays = df['nab_decay'].unique()

#%% Plotting

for var_day in new_variant_days:
    df_to_use = df[df['new_variant_day']==var_day]
    
    # Wrangle the data
    yvals = list(vaccines.keys())
    xvals = variants
    y = np.arange(len(yvals))
    x = np.arange(len(xvals))
    ylabels = [vaccines[k]['label'] for k in yvals]
    xlabels = variant_labels

    # Calculate t-tests for deaths
    from scipy.stats import ttest_ind

    z_deaths_ttest = np.full((len(yvals), len(xvals)), fill_value=0, dtype=np.float32)

    for i, vx in enumerate(yvals):
        for j, var in enumerate(xvals):
            status_quo = df_to_use[df_to_use['vaccine'] == 'Status quo']
            status_quo = status_quo[status_quo['next_variant'] == var]
            if i > 0:
                comparison = df_to_use[df_to_use['vaccine'] == vx]
                comparison = comparison[comparison['next_variant'] == var]
                res = ttest_ind(status_quo['deaths'], comparison['deaths'])
                z_deaths_ttest[i, j] = res.pvalue

    # Do the grouping
    df_to_use = df_to_use[df_to_use['deaths'] > 0]
    df_to_use = df_to_use[df_to_use['infections'] > 0]

    dfg = df_to_use.groupby(['vaccine', 'next_variant'])
    dfmean = dfg.mean().reset_index()
    dfmedian = dfg.median().reset_index()

    # First up, doses
    z_doses = dfmean.pivot('vaccine', 'next_variant', 'doses').reindex(xvals, axis=1)
    z_doses = z_doses.reindex(yvals, axis=0)
    z_doses = z_doses.values

    # Now deaths
    z_deaths = dfmean.pivot('vaccine', 'next_variant', 'deaths').reindex(xvals, axis=1)
    z_deaths = z_deaths.reindex(yvals, axis=0)
    z_deaths = z_deaths.values

    # Deaths averted
    z_deaths_averted = z_deaths[0,:] - z_deaths
    # Percent of deaths averted
    z_perc_deaths_averted = 100 * z_deaths_averted / z_deaths[0,:]
    # Now doses per deaths averted
    z_doses_deaths_averted = z_doses / z_deaths_averted

    heatmap(
        data=z_deaths,
        cmap=bad_cmap,
        zlabel='Cumulative deaths',
        suptitle=f'Cumulative deaths,',
        filename=f'vx_rollout.png',
        row=0,
        time_lab=var_day
    )

    heatmap(
        data=z_perc_deaths_averted,
        cmap=good_cmap,
        zlabel='Percent of deaths averted',
        suptitle=f'Percent of deaths averted,',
        filename=f'vx_rollout.png',
        row=1,
        time_lab=var_day
    )
        
    heatmap(
        data=z_doses_deaths_averted,
        cmap=neutral_cmap,
        zlabel='Doses per death averted',
        suptitle=f'Doses per death averted',
        filename=f'vx_rollout.png',
        row=2,
        time_lab=var_day
    )


pl.show()


print('Done')