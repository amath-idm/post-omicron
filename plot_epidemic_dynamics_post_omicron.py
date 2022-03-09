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
fn = f'{str(sc.path(resdir))}/post-omicron-cov-sweep_data_for_run_plot.obj'
df = sc.load(fn)

figdir = './figs'

ncolors = 8
from matplotlib import cm
# Rejig data
d = sc.objdict()
d_inf = 2 # Index for Delta
o_ind = 3  # Index for Omicron
v_ind = 4  # Index for next variant
o2_day_ind = 123  # Index for day after Omicron peak

variants = ['None', 'Emerged from Omicron', 'Emerged from WT', 'New cluster',
    'Emerged from Omicron, more severe', 'Emerged from WT, more severe',
    'New cluster, more severe'
]

variants_inf = ['None', 'Emerged from Omicron', 'Emerged from WT', 'New cluster']

variant_labels = ['No new variant', 'Emerged from Omicron', 'Emerged from WT', 'New cluster',
    'Emerged from Omicron, more severe', 'Emerged from WT, more severe',
    'New cluster, more severe'
]

variant_inf_labels = ['No new variant', 'Emerged from Omicron', 'Emerged from WT', 'New cluster'
]

variant_time = ['2022-02-25', '2022-04-25', '2022-08-25']

keys = ['inf', 'sev', 'inf_ci']
for key in keys:
    d[key] = sc.objdict(defaultdict=dict)
    for var in variants:
        d[key][var] = sc.objdict(defaultdict=dict)
        for var_date in variant_time:
            d[key][var][var_date] = sc.objdict(defaultdict=dict)

for i, inf in enumerate(df['new_infections_by_variant']):
    var = df['next_variant'][i]
    if df['vaccine_boost'][i] == 0 and df['vaccine_prime'][i] == 0:
        var_day = df['new_variant_day'][i]
        new_infs = inf.values[:,o2_day_ind:].sum()
        new_infs_high = inf.high[:,o2_day_ind:].sum()
        new_infs_low = inf.low[:,o2_day_ind:].sum()
        new_sevs = df['new_severe_by_variant'][i].values[:,o2_day_ind:].sum()
        d['inf'][var][var_day] = new_infs
        d['inf_ci'][var][var_day] = new_infs_high-new_infs_low
        d['sev'][var][var_day] = new_sevs


#%% Plotting
sc.options(dpi=150)
sc.options(font='Avenir')
colors = sc.gridcolors(ncolors)
fig = pl.figure(figsize=(7,7))
xoff = 0.3
dx   = 0.65
dy   = 0.4
ax = sc.objdict()
ax.inf = fig.add_axes([xoff, 0.55, dx, dy])
ax.sev = fig.add_axes([xoff, 0.05, dx, dy])

# Plot bar charts
n = len(variants_inf)*len(variant_time)
# Set the width of the bars
wd = .8
y_pos = n - 3*(np.arange(0, len(variants_inf), 1))

for v, var in enumerate(variants_inf):
    for j, day in enumerate(variant_time):
        val = d['inf'][var][day]
        if v == 0:
            ax['inf'].barh(y_pos[v] - (wd * j), val, color=colors[j], label=day)
        else:
            ax['inf'].barh(y_pos[v] - (wd * j), val, color=colors[j])

ax['inf'].set_yticks(y_pos-(wd/3))
ax['inf'].set_yticklabels(variant_inf_labels)

n = len(variants)*len(variant_time)
# Set the width of the bars
wd = .8
y_pos = n - 3*(np.arange(0, len(variants), 1))

for v, var in enumerate(variants):
    for j, day in enumerate(variant_time):
        val = d['sev'][var][day]
        if v == 0:
            ax['sev'].barh(y_pos[v] - (wd * j), val, color=colors[j], label=day)
        else:
            ax['sev'].barh(y_pos[v] - (wd * j), val, color=colors[j])

ax['sev'].set_yticks(y_pos-(wd/3))
ax['sev'].set_yticklabels(variant_labels)
# ax.inf.xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
# ax.sev.xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.inf.set_title('New infections')
ax.sev.set_title('Severe cases')
ax.sev.legend(frameon=False, title='Date of next variant introduction')


# Tidy up
fig.show()
fig.savefig(str(sc.path(figdir) / 'inf_by_variant.png'), bbox_inches='tight')


sc.options(dpi=150)
sc.options(font='Avenir')
colors = cm.rainbow(np.linspace(0, 1, ncolors))

new_variants_inf = [
    'Emerged from Omicron', 'Emerged from WT', 'New cluster'
]

new_variants_sev = [
    'Emerged from Omicron', 'Emerged from WT', 'New cluster',
    'Emerged from Omicron, more severe', 'Emerged from WT, more severe',
    'New cluster, more severe'
]

var_ax_code = {
    'Emerged from Omicron':0,
    'Emerged from WT':1,
    'New cluster':2,
    'Emerged from Omicron, more severe':0,
    'Emerged from WT, more severe':1,
    'New cluster, more severe': 2
}

variant_timing = ['2022-02-25', '2022-04-25', '2022-08-25']

fig, axes = pl.subplots(3, 3, figsize=(10, 10), sharex='col', sharey='row')

x = df['datevec_sim'][o2_day_ind:]

omicron_lines = []
variant_lines = []
for i, inf in enumerate(df['new_infections_by_variant']):
    if df['vaccine_boost'][i] == 0 and df['vaccine_prime'][i] == 0:
        if df['next_variant'][i] in new_variants_inf:
            var = df['next_variant'][i]
            var_ax = var_ax_code[var]
            var_date = df['new_variant_day'][i]
            date_ax = variant_timing.index(var_date)
            ax = axes[var_ax, date_ax]

            factor_inf = 100 / inf.values[o_ind, 0:o2_day_ind].max()
            l1, = ax.plot(x, inf.values[o_ind, o2_day_ind:] * factor_inf, c=colors[0], label='Omicron')
            ax.fill_between(x, (inf.low[o_ind, o2_day_ind:] * factor_inf), (inf.high[o_ind, o2_day_ind:]*factor_inf),
                            color=colors[0], alpha=0.3)
            omicron_lines.append(l1)
            l2, = ax.plot(x, inf.values[v_ind, o2_day_ind:] * factor_inf, c=colors[var_ax + 1], label=var)
            ax.fill_between(x, (inf.low[v_ind, o2_day_ind:] * factor_inf), (inf.high[v_ind, o2_day_ind:]*factor_inf),
                            color=colors[var_ax + 1], alpha=0.3)
            variant_lines.append(l2)
            ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
            ax.tick_params(axis='x', rotation=90)

pl.legend([omicron_lines[0], variant_lines[0], variant_lines[3], variant_lines[6]],
          ['Omicron', 'Emerged from Omicron', 'Emerged from WT', 'New cluster'], loc='right', bbox_to_anchor=(2.2, 2),
          title='Variant')
axes[0,0].set_title(f'Introduced on \n{variant_timing[0]}')
axes[0,1].set_title(f'Introduced on \n{variant_timing[1]}')
axes[0,2].set_title(f'Introduced on \n{variant_timing[2]}')
axes[0,0].set_ylabel('% of Omicron peak')
axes[1,0].set_ylabel('% of Omicron peak')
axes[2,0].set_ylabel('% of Omicron peak')
# fig.suptitle('New Infections')
fig.subplots_adjust(right=0.75)
fig.show()
fig.savefig(str(sc.path(figdir) / 'post_omicron_inf.png'), bbox_inches='tight')


fig, axes = pl.subplots(3, 3, figsize=(10, 10), sharex='col', sharey='row')

x = df['datevec_sim'][o2_day_ind:]

omicron_lines = []
variant_lines = []
sev_lines = []
for i, sev in enumerate(df['new_severe_by_variant']):
    if df['vaccine_boost'][i] == 0 and df['vaccine_prime'][i] == 0:
        if df['next_variant'][i] in new_variants_sev:
            var = df['next_variant'][i]
            var_ax = var_ax_code[var]
            var_date = df['new_variant_day'][i]
            date_ax = variant_timing.index(var_date)
            ax = axes[var_ax, date_ax]

            factor_sev = 100 / sev.values[o_ind, 0:o2_day_ind].max()
            l1, = ax.plot(x, sev.values[o_ind, o2_day_ind:] * factor_sev, c=colors[0], label='Omicron')
            ax.fill_between(x, (sev.low[o_ind, o2_day_ind:] * factor_sev), (sev.high[o_ind, o2_day_ind:]*factor_sev),
                            color=colors[0], alpha=0.3)
            omicron_lines.append(l1)

            if var in new_variants_inf:
                l2, = ax.plot(x, sev.values[v_ind, o2_day_ind:] * factor_sev, c=colors[var_ax + 1], label=var)
                ax.fill_between(x, (sev.low[v_ind, o2_day_ind:] * factor_sev), (sev.high[v_ind, o2_day_ind:]*factor_sev),
                            color=colors[var_ax + 1], alpha=0.3)
                variant_lines.append(l2)
            else:
                l3, = ax.plot(x, sev.values[v_ind, o2_day_ind:] * factor_sev, c=colors[var_ax + 1], linestyle=(0,(1,10)))
                ax.fill_between(x, (sev.low[v_ind, o2_day_ind:] * factor_sev),
                                (sev.high[v_ind, o2_day_ind:] * factor_sev),
                                color=colors[var_ax + 1], alpha=0.3)
                sev_lines.append(l3)
            ax.tick_params(axis='x', rotation=90)
            ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))

legend2 = pl.legend([omicron_lines[0], variant_lines[0], variant_lines[1], variant_lines[2]],
          ['Omicron', 'Emerged from Omicron', 'Emerged from WT', 'New cluster'], loc='right', bbox_to_anchor=(2.2, 2.1),
          title='Variant')
pl.legend([variant_lines[0], sev_lines[0]], ['mild', 'virulent'], loc='right', bbox_to_anchor=(2, 1.6), title='Severity of variant')

pl.gca().add_artist(legend2)
# pl.xticks(rotation = 90)
axes[0,0].set_title(f'Introduced on \n{variant_timing[0]}')
axes[0,1].set_title(f'Introduced on \n{variant_timing[1]}')
axes[0,2].set_title(f'Introduced on \n{variant_timing[2]}')
axes[0,0].set_ylabel('% of Omicron peak')
axes[1,0].set_ylabel('% of Omicron peak')
axes[2,0].set_ylabel('% of Omicron peak')
# fig.suptitle('New Severe')
fig.subplots_adjust(right=0.75)

fig.show()
fig.savefig(str(sc.path(figdir) / 'post_omicron_sev.png'), bbox_inches='tight')

#%% Plotting
sc.options(dpi=150)
sc.options(font='Avenir')
colors = sc.gridcolors(ncolors)
fig, ax = pl.subplots(figsize=(7,7))

# Plot bar charts
n = len(variants_inf)
# Set the width of the bars
y_pos = n-(np.arange(0, len(variants_inf), 1))

for v, var in enumerate(variants_inf):
    for j, day in enumerate(variant_time):
        if j == 0:
            val = d['inf'][var][day]/1e6
            err = d['inf_ci'][var][day]/1e6/2
            ax.barh(y_pos[v], val, xerr=err, capsize=7, color=colors[j])

ax.set_yticks(y_pos)
ax.set_yticklabels(variant_inf_labels)
ax.xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.set_xlabel('Cumulative new infections (millions)')
ax.set_ylabel('New variant')
# fig.subplots_adjust(bottom=0.28)
fig.subplots_adjust(left=0.3)
# Tidy up
fig.show()
fig.savefig(str(sc.path(figdir) / 'inf_by_variant_v2.png'), bbox_inches='tight')