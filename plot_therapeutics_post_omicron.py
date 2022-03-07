
'''
Script to plot combo of therapeutics and next generation vaccines
'''

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl


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

def heatmap(data, zlabel, suptitle, filename, cmap, threshold=0.5, prime_lab=None, row=None):
    ''' Plot a matrix of results as a heat map '''

    prime_col_map = {'Boost 100% vaccinated, prime 10% unvaccinated': 0,
                     'Boost 100% vaccinated, prime 50% unvaccinated': 1,
                     'Boost 100% vaccinated, prime 100% unvaccinated': 2}
    col = prime_col_map[prime_lab]
    print(f'Plotting: {prime_lab}, row={row}, col={col}')
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
        ax.set_xticks(x, xvals, rotation=90)
        ax.set_xlabel('Therapeutic strategy', fontweight='bold')
    else:
        ax.set_xticks([])
    if row == 0:
        ax.set_title(f'{prime_lab}')
    if col == 0:
        ax.set_yticks(y, yvals)
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
    cbargs = dict(labelpad=25, rotation=270, fontweight='bold')
    cpos = sc.dcp(cax_pos)
    cpos[1] -= row*0.268
    cax = fig.add_axes(cpos)
    cb = fig.colorbar(im, ticklocation='right', orientation='vertical', cax=cax)
    cb.set_label(label, **cbargs)
    sc.commaticks(cax)
    return

# Load data
sc.options(dpi=100)
sc.fonts(add=sc.thisdir(aspath=True) / 'avenir')
sc.options(font='Avenir')

# Load data
resfolder = './results'

df = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-therapeutics_vx_rollout_data.obj')))
df = pd.DataFrame.from_dict(df)
figdir = './figs'

vx_prime = list(set(df['vaccine_prime']))
treatments = ['None', '30% of 60+ symptomatic', '30% of 18+ symptomatic']
vx_prime_labels = ['Boost 100% vaccinated, prime 10% unvaccinated', 'Boost 100% vaccinated, prime 50% unvaccinated', 'Boost 100% vaccinated, prime 100% unvaccinated']
vx_labels = ['Status quo', 'Durable', 'Broadly neutralizing', 'Broadly neutralizing and durable']

for i, prime in enumerate(vx_prime):
    df_to_use = df[df['vaccine_prime'] == prime]
    prime_label = vx_prime_labels[i]

    yvals = vx_labels
    xvals = treatments
    y = np.arange(len(yvals))
    x = np.arange(len(xvals))

    from scipy.stats import ttest_ind

    z_deaths_ttest = np.full((len(yvals), len(xvals)), fill_value=0, dtype=np.float32)

    for i, vx in enumerate(yvals):
        for j, tx in enumerate(xvals):
            status_quo = df_to_use[df_to_use['vaccine'] == 'Status quo']
            status_quo = status_quo[status_quo['treatment'] == 'None']
            if i > 0:
                comparison = df_to_use[df_to_use['vaccine'] == vx]
                comparison = comparison[comparison['treatment'] == tx]
                res = ttest_ind(status_quo['deaths'], comparison['deaths'])
                z_deaths_ttest[i, j] = res.pvalue

    # Do the grouping
    df_to_use = df_to_use[df_to_use['deaths'] > 0]
    df_to_use = df_to_use[df_to_use['infections'] > 0]

    dfg = df_to_use.groupby(['vaccine', 'treatment'])
    dfmean = dfg.mean().reset_index()
    dfmedian = dfg.median().reset_index()

    # First up, doses
    z_doses = dfmean.pivot('vaccine', 'treatment', 'doses').reindex(xvals, axis=1)
    z_doses = z_doses.reindex(yvals, axis=0)
    z_doses = z_doses.values
    z_additional_doses = z_doses - z_doses[0, 0]

    z_vaccinated = z_additional_doses/2

    # Now treated
    z_treated = dfmean.pivot('vaccine', 'treatment', 'n_treated').reindex(xvals, axis=1)
    z_treated = z_treated.reindex(yvals, axis=0)
    z_treated = z_treated.values

    z_vaxed_and_treated = z_treated + z_vaccinated

    # Now deaths
    z_deaths = dfmean.pivot('vaccine', 'treatment', 'deaths').reindex(xvals, axis=1)
    z_deaths = z_deaths.reindex(yvals, axis=0)
    z_deaths = z_deaths.values

    # Deaths averted
    z_deaths_averted = z_deaths[0, 0] - z_deaths
    # Percent of deaths averted
    z_perc_deaths_averted = 100 * z_deaths_averted / z_deaths[0, 0]
    # Now doses per deaths averted
    z_doses_deaths_averted = z_doses / z_deaths_averted

    # Now NNT
    z_nnt = z_vaxed_and_treated / z_deaths_averted

    heatmap(
        data=z_deaths,
        cmap=bad_cmap,
        zlabel='Cumulative deaths',
        suptitle=f'Cumulative deaths,',
        filename=f'vx_tx_rollout.png',
        row=0,
        prime_lab=prime_label
    )

    heatmap(
        data=z_perc_deaths_averted,
        cmap=good_cmap,
        zlabel='Percent of deaths averted',
        suptitle=f'Percent of deaths averted,',
        filename=f'vx_tx_rollout.png',
        row=1,
        prime_lab=prime_label
    )

    heatmap(
        data=z_nnt,
        cmap=neutral_cmap,
        zlabel='Individuals treated and/or\nvaccinated per death averted',
        suptitle=f'Individuals treated and/or\nvaccinated per death averted',
        filename=f'vx_tx_rollout.png',
        row=2,
        prime_lab=prime_label
    )

print('Done')