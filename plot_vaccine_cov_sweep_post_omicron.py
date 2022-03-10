
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
figs = sc.objdict()
axes = sc.objdict()

fig2,axs2 = pl.subplots(nrows=3, ncols=2, figsize=figsize)


def heatmap(data, zlabel, suptitle, filename, cmap, threshold=0.5, time_lab=None, var_name=None, row=None):
    ''' Plot a matrix of results as a heat map '''
    time_col_map = {'2022-02-25': 0, '2022-04-25': 1, '2022-08-25': 2}
    col = time_col_map[time_lab]
    print(f'Plotting: {var_name}, time={time_lab}, row={row}, col={col}')
    if var_name in figs.keys():
        axs = axes[var_name]
        ax = axs[row, col]
        fig = figs[var_name]
    else:
        fig, axs = pl.subplots(nrows=3, ncols=3, figsize=figsize)
        figs[var_name] = fig
        axes[var_name] = axs
        ax = axs[row, col]

    # Remove non-significant data
    if zlabel != 'Cumulative deaths':
        for i in x:
            for j in y:
                if z_deaths_ttest[j, i] > 0.05:
                    data[j, i] = np.nan  # To turn off these cells

        data[0,0] = np.nan

    # Plot main axes
    im = ax.imshow(data, cmap=cmap, origin='lower')

    # Reconcile color limits
    cmin, cmax = im.get_clim()
    crossover = (cmax - cmin) * threshold + cmin
    if col > 0:
        other_im = axs[row, 0].get_images()[0]
        other_cmin, other_cmax = other_im.get_clim()
        cmax = max(cmax, other_cmax)
        im.set_clim((0, cmax))
        other_im.set_clim(0, cmax)
        crossover = (cmax - other_cmin) * threshold + other_cmin

    # Handle labels
    for i in x:
        for j in y:

            # Label text
            k = data[j, i]
            if zlabel == 'Cumulative deaths':
                label = f'{k:0,.0f}'
            else:
                if (j+i) == 0:
                    label = 'ref'
                else:
                    if not np.isfinite(k):
                        label = 'NS'
                    else:
                        label = f'{k:0,.0f}'

                        if 'Percent' in suptitle:
                            if (j+i) > 0:
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
        ax.set_xticks(x, xlabels)
        ax.set_xlabel('Primary series vaccine coverage', fontweight='bold')
    else:
        ax.set_xticks([])
    if row == 0:
        ax.set_title(f'{var_name} on\n{time_lab}')
    if col == 0:
        ax.set_yticks(y, ylabels)
        ax.set_ylabel('Booster dose vaccine coverage', fontweight='bold')
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


def heatmap_by_var(data, zlabel, suptitle, filename, cmap, threshold=0.5, var_lab=None, time_lab=None, var_name=None, row=None):
    ''' Plot a matrix of results as a heat map '''
    var_col_map = {'Emerged from WT, more severe': 0, 'New cluster, more severe': 1}
    col = var_col_map[var_lab]
    print(f'Plotting: {var_name}, row={row}, col={col}')
    ax = axs2[row, col]

    # Remove non-significant data
    if zlabel != 'Cumulative deaths':
        for i in x:
            for j in y:
                if z_deaths_ttest[j, i] > 0.05:
                    data[j, i] = np.nan  # To turn off these cells
                if z_perc_deaths_averted[j,i] == 0:
                    data[j, i] = np.nan  # To turn off these cells

        data[0,0] = np.nan

    # Plot main axes
    im = ax.imshow(data, cmap=cmap, origin='lower')

    # Reconcile color limits
    cmin, cmax = im.get_clim()
    crossover = (cmax - cmin) * threshold + cmin
    if col >0:
        other_im = axs2[row, 0].get_images()[0]
        other_cmin, other_cmax = other_im.get_clim()
        cmax = max(cmax, other_cmax)
        im.set_clim((0, cmax))
        other_im.set_clim(0, cmax)
        crossover = (cmax - other_cmin) * threshold + other_cmin


    # Handle labels
    for i in x:
        for j in y:

            # Label text
            k = data[j, i]
            if zlabel == 'Cumulative deaths':
                label = f'{k:0,.0f}'
            else:
                if (j+i) == 0:
                    label = 'ref'
                else:
                    if not np.isfinite(k):
                        label = 'NS'
                    else:
                        label = f'{k:0,.0f}'

                        if 'Percent' in suptitle:
                            if (j+i) > 0:
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
        ax.set_xticks(x, xlabels)
        ax.set_xlabel('Primary series vaccine coverage', fontweight='bold')
    else:
        ax.set_xticks([])
    if row == 0:
        ax.set_title(f'{var_name} on\n{time_lab}')
    if col == 0:
        ax.set_yticks(y, ylabels)
        ax.set_ylabel('Booster dose vaccine coverage', fontweight='bold')
    else:
        ax.set_yticks([])
    sc.boxoff(ax=ax)
    ax.axis('auto')

    # Tidying
    if col == 1:
        colorbar(fig2, im, zlabel, row)
    fig2.subplots_adjust(**pad)
    if do_save and row == 2 and col == 1:
        sc.savefig(str(sc.path(figdir) / filename), fig=fig2)

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

new_variant_days = ['2022-02-25', '2022-04-25', '2022-08-25']
vaccine_prime = [0, 0.2, 0.4, 0.6, 0.8, 1]
vaccine_boost = [0, 0.2, 0.4, 0.6, 0.8, 1]

# Load data
sc.fonts(add='C:\\Users\\jamieco\\PycharmProjects\\avenir')

# Load data
resfolder = './results'
df = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-cov-sweep_vx_rollout_data.obj')))
df = pd.DataFrame.from_dict(df)
nab_decays = df['nab_decay'].unique()

#%% Plotting

variants = ['None', 'Emerged from Omicron', 'Emerged from WT', 'New cluster',
    'Emerged from Omicron, more severe', 'Emerged from WT, more severe',
    'New cluster, more severe'
]

variant_dict = {
    'None': {
        'label': 'No new variant',
        'fn': 'none',
    },
    'Emerged from Omicron': {
        'label': 'Emerged from Omicron',
        'fn': 'from_omicron',
    },
    'Emerged from WT': {
        'label': 'Emerged from WT',
        'fn': 'from_WT',
    },
    'New cluster': {
        'label': 'New antigenic cluster',
        'fn': 'new_cluster',
    },
    'Emerged from Omicron, more severe': {
        'label': 'Emerged from Omicron\n more severe',
        'fn': 'from_omicron_severe',
    },
    'Emerged from WT, more severe': {
        'label': 'Emerged from WT\n more severe',
        'fn': 'from_WT_severe',
    },
    'New cluster, more severe': {
        'label': 'New cluster\nmore severe',
        'fn': 'new_cluster_severe',
    },

}

for variant in variants:

    for var_day in new_variant_days:
        df_to_use = df[df['new_variant_day'] == var_day]
        df_to_use = df_to_use[df_to_use['next_variant'] == variant]

        variant_label = variant_dict[variant]['label']
        fn_to_add = variant_dict[variant]['fn']

        # Wrangle the data
        yvals = vaccine_boost
        xvals = vaccine_prime
        y = np.arange(len(yvals))
        x = np.arange(len(xvals))
        xlabels = [f'{int((.47 + .53*x)*100)}%' for x in xvals]
        ylabels = [f'{int(y*100)}%' for y in yvals]

        # Calculate t-tests for deaths
        from scipy.stats import ttest_ind

        z_deaths_ttest = np.full((len(yvals), len(xvals)), fill_value=0, dtype=np.float32)

        for i, vx in enumerate(yvals):
            for j, var in enumerate(xvals):
                status_quo = df_to_use[df_to_use['vaccine_boost'] == 0]
                status_quo = status_quo[status_quo['vaccine_prime'] == 0]
                if i > 0:
                    comparison = df_to_use[df_to_use['vaccine_boost'] == vx]
                    comparison = comparison[comparison['vaccine_prime'] == var]
                    res = ttest_ind(status_quo['deaths'], comparison['deaths'])
                    z_deaths_ttest[i, j] = res.pvalue

        # Do the grouping
        df_to_use = df_to_use[df_to_use['deaths'] > 0]
        df_to_use = df_to_use[df_to_use['infections'] > 0]

        dfg = df_to_use.groupby(['vaccine_prime', 'vaccine_boost'])
        dfmean = dfg.mean().reset_index()
        dfmedian = dfg.median().reset_index()

        # First up, doses
        z_doses = dfmean.pivot('vaccine_boost', 'vaccine_prime', 'doses').reindex()
        z_doses = z_doses.reindex()
        z_doses = z_doses.values
        z_additional_doses = z_doses - z_doses[0,0]

        # Now deaths
        z_deaths = dfmean.pivot('vaccine_boost', 'vaccine_prime', 'deaths').reindex()
        z_deaths = z_deaths.reindex()
        z_deaths = z_deaths.values

        # Deaths averted
        z_deaths_averted = z_deaths[0, 0] - z_deaths
        # Percent of deaths averted
        z_perc_deaths_averted = 100 * z_deaths_averted / z_deaths[0, 0]
        # Now doses per deaths averted
        z_doses_deaths_averted = z_additional_doses / z_deaths_averted

        heatmap(
            data=z_deaths,
            cmap=bad_cmap,
            zlabel='Cumulative deaths',
            suptitle=f'Cumulative deaths,',
            filename=f'vx_rollout_cov_sweep_{fn_to_add}.png',
            row=0,
            time_lab=var_day,
            var_name=variant_label
        )

        heatmap(
            data=z_perc_deaths_averted,
            cmap=good_cmap,
            zlabel='Percent of deaths averted',
            suptitle=f'Percent of deaths averted,',
            filename=f'vx_rollout_cov_sweep_{fn_to_add}.png',
            row=1,
            time_lab=var_day,
            var_name=variant_label
        )

        heatmap(
            data=z_doses_deaths_averted,
            cmap=neutral_cmap,
            zlabel='Doses per death averted',
            suptitle=f'Doses per death averted',
            filename=f'vx_rollout_cov_sweep_{fn_to_add}.png',
            row=2,
            time_lab=var_day,
            var_name=variant_label
        )


variants_to_plot = ['Emerged from WT, more severe', 'New cluster, more severe'
]

for variant in variants_to_plot:
    var_day = '2022-04-25'

    df_to_use = df[df['new_variant_day'] == var_day]
    df_to_use = df_to_use[df_to_use['next_variant'] == variant]

    variant_label = variant_dict[variant]['label']

    # Wrangle the data
    yvals = vaccine_boost
    xvals = vaccine_prime
    y = np.arange(len(yvals))
    x = np.arange(len(xvals))
    xlabels = [f'{int((.47 + .53 * x) * 100)}%' for x in xvals]
    ylabels = [f'{int(y * 100)}%' for y in yvals]

    # Calculate t-tests for deaths
    from scipy.stats import ttest_ind

    z_deaths_ttest = np.full((len(yvals), len(xvals)), fill_value=0, dtype=np.float32)

    for i, vx in enumerate(yvals):
        for j, var in enumerate(xvals):
            status_quo = df_to_use[df_to_use['vaccine_boost'] == 0]
            status_quo = status_quo[status_quo['vaccine_prime'] == 0]
            if i > 0:
                comparison = df_to_use[df_to_use['vaccine_boost'] == vx]
                comparison = comparison[comparison['vaccine_prime'] == var]
                res = ttest_ind(status_quo['deaths'], comparison['deaths'])
                z_deaths_ttest[i, j] = res.pvalue

    # Do the grouping
    df_to_use = df_to_use[df_to_use['deaths'] > 0]
    df_to_use = df_to_use[df_to_use['infections'] > 0]

    dfg = df_to_use.groupby(['vaccine_prime', 'vaccine_boost'])
    dfmean = dfg.mean().reset_index()
    dfmedian = dfg.median().reset_index()

    # First up, doses
    z_doses = dfmean.pivot('vaccine_boost', 'vaccine_prime', 'doses').reindex()
    z_doses = z_doses.reindex()
    z_doses = z_doses.values
    z_additional_doses = z_doses - z_doses[0, 0]

    # Now deaths
    z_deaths = dfmean.pivot('vaccine_boost', 'vaccine_prime', 'deaths').reindex()
    z_deaths = z_deaths.reindex()
    z_deaths = z_deaths.values

    # Deaths averted
    z_deaths_averted = z_deaths[0, 0] - z_deaths
    below0 = (z_deaths_averted < 0).sum()
    if below0:
        print(f'Warning: {below0} entries for deaths averted were below zero')
        z_deaths_averted = np.maximum(0, z_deaths_averted)
    # Percent of deaths averted
    z_perc_deaths_averted = 100 * z_deaths_averted / z_deaths[0, 0]
    # Now doses per deaths averted
    z_doses_deaths_averted = z_additional_doses / z_deaths_averted

    heatmap_by_var(
        data=z_deaths,
        cmap=bad_cmap,
        zlabel='Cumulative deaths',
        suptitle=f'Cumulative deaths,',
        filename=f'vx_rollout_cov_sweep.png',
        row=0,
        time_lab=var_day,
        var_lab=variant,
        var_name=variant_label
    )

    heatmap_by_var(
        data=z_perc_deaths_averted,
        cmap=good_cmap,
        zlabel='Percent of deaths averted',
        suptitle=f'Percent of deaths averted,',
        filename=f'vx_rollout_cov_sweep.png',
        row=1,
        time_lab=var_day,
        var_lab=variant,
        var_name=variant_label
    )

    heatmap_by_var(
        data=z_doses_deaths_averted,
        cmap=neutral_cmap,
        zlabel='Doses per death averted',
        suptitle=f'Doses per death averted',
        filename=f'vx_rollout_cov_sweep.png',
        row=2,
        time_lab=var_day,
        var_name=variant_label,
        var_lab=variant,
    )

print('Done')