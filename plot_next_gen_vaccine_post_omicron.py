
'''
Script to plot next generation vaccine impact
'''

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.dates as mdates
import matplotlib.ticker as mtick



vaccine_breadth = [0, 1]
vaccine_durability = [0, 1]

vaccine_breadth_label = {
    0: 'antigen-specific',
    1: 'broadly neutralizing'
}

vaccine_durability_label = {
    0: 'short-lived',
    1: 'long-lasting'
}

# Load data
sc.options(dpi=100)
sc.fonts(add=sc.thisdir(aspath=True) / 'avenir')
sc.options(font='Avenir')

# Load data
resfolder = './results'
df1 = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-next-gen-vaccines_vx_rollout_data.obj')))
df2 = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-next-gen-vaccines_data_for_run_plot.obj')))
df1 = pd.DataFrame.from_dict(df1)
figdir = './figs'

vx_breadth = list(set(df2['vaccine_breadth']))
vx_durability = list(set(df2['vaccine_durability']))
vx_prime = list(set(df2['vaccine_prime']))

z_perc_deaths_averted_dict = []
z_doses_deaths_averted_dict = []
for prime in vx_prime:

    df_to_use = df1[df1['vaccine_prime'] == prime]
    dfg = df_to_use.groupby(['vaccine_breadth', 'vaccine_durability'])
    dfmean = dfg.mean().reset_index()
    dfmedian = dfg.median().reset_index()

    # First up, doses
    z_doses = dfmean.pivot('vaccine_breadth', 'vaccine_durability', 'doses').reindex()
    z_doses = z_doses.reindex()
    z_doses = z_doses.values
    z_additional_doses = z_doses - z_doses[0, 0]

    # Now deaths
    z_deaths = dfmean.pivot('vaccine_breadth', 'vaccine_durability', 'deaths').reindex()
    z_deaths = z_deaths.reindex()
    z_deaths = z_deaths.values

    # Deaths averted
    z_deaths_averted = z_deaths[0, 0] - z_deaths
    # Percent of deaths averted
    z_perc_deaths_averted = 100 * z_deaths_averted / z_deaths[0, 0]
    z_perc_deaths_averted_dict.append(z_perc_deaths_averted)
    # Now doses per deaths averted
    z_doses_deaths_averted = z_additional_doses / z_deaths_averted
    z_doses_deaths_averted_dict.append(z_doses_deaths_averted)


#%% Plotting

# fig, axs = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
#
# for i, ax in enumerate(axs):
#     x = df2['datevecs'][i]
#     for j, nabs in enumerate(df2['pop_nabs'][i]):
#         if i == 1:
#             ax.plot(x, nabs, label=df2['imm_source'][j])
#         else:
#             ax.plot(x, nabs)
#     ax.set_title(f'{vaccine_breadth_label[vx_breadth[i]]}, {vaccine_durability_label[vx_durability[i]]} vaccine')
#
# axs[1].legend()
# fig.show()
# sc.savefig(str(sc.path(figdir) / 'pop_nabs_next_gen.png'), fig=fig)
#
# fig, axs = plt.subplots(4, 1, figsize=(8, 10), sharex=True)
#
# for i, ax in enumerate(axs):
#     x = df2['datevec_sim'][1:]
#     for j, infs in enumerate(df2['new_infections_by_variant'][i]):
#         if i == 1:
#             ax.plot(x, infs[1:], label=df2['imm_source'][j])
#         else:
#             ax.plot(x, infs[1:])
#     ax.set_title(f'{vaccine_breadth_label[vx_breadth[i]]}, {vaccine_durability_label[vx_durability[i]]} vaccine')
#     ax.set_ylabel('New Infections')
#     ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
#     ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
#
#     ax = axs[i].twinx()
#     for j, sevs in enumerate(df2['new_severe_by_variant'][i]):
#         ax.plot(x, sevs[1:], linestyle='--')
#     ax.set_ylabel('New Severe Cases')
#     ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
#
# axs[1].legend()
# fig.show()
# sc.savefig(str(sc.path(figdir) / 'infs_next_gen.png'), fig=fig)
#
#
fig, axv = plt.subplots(3, 1, figsize=(8, 10))

# FIRST up TS with no new vaccination
ncolors = 8
from matplotlib import cm
colors = cm.rainbow(np.linspace(0, 1, ncolors))
ax = axv[0]
x = df2['datevec_sim'][136:]
for j, infs in enumerate(df2['new_infections_by_variant'][0]):
    if j == 3:
        ax.plot(x, infs[136:], label=df2['imm_source'][j], color=colors[0])
        ax.fill_between(x, (df2['new_infections_by_variant'][0].low[j,136:]), (df2['new_infections_by_variant'][0].high[j,136:]), color=colors[0], alpha=0.3)
    elif j == 4:
        ax.plot(x, infs[136:], linestyle='--', label=df2['imm_source'][j], color=colors[0])
        ax.fill_between(x, (df2['new_infections_by_variant'][0].low[j, 136:]),
                        (df2['new_infections_by_variant'][0].high[j, 136:]), color=colors[0], alpha=0.3)
ax.set_ylabel('New Infections')
ax.legend(loc='upper left')
# ax.axvline(x[136], linestyle='--', color='black')
# ax.text(x[140], 2.5e6, 'Date next generation\nvaccines would begin roll-out')
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))

ax.set_title('New infections and deaths (no additional vaccination)')

ax = axv[0].twinx()
ax.plot(x, df2['new_deaths'][0].values[136:], color=colors[1])
ax.fill_between(x, (df2['new_deaths'][0].low[136:]), (df2['new_deaths'][0].high[136:]), alpha=0.3, color=colors[1])
ax.set_ylabel('New Deaths', color=colors[1])
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))

colors = sc.gridcolors(5)

ax = axv[1]
# Plot bar charts
n = 3*3
wd=0.8
y_pos = n - 3*(np.arange(0, 3, 1))
ax.barh(y_pos[0]-(wd*0), z_perc_deaths_averted_dict[0][0,1], color=colors[0], label='Durable')
ax.barh(y_pos[0]-(wd*1), z_perc_deaths_averted_dict[0][1,0], color=colors[1], label='Broadly neutralizing')
ax.barh(y_pos[0]-(wd*2), z_perc_deaths_averted_dict[0][1,1], color=colors[2], label='Durable & broadly neutralizing')
ax.barh(y_pos[1]-(wd*0), z_perc_deaths_averted_dict[1][0,1], color=colors[0])
ax.barh(y_pos[1]-(wd*1), z_perc_deaths_averted_dict[1][1,0], color=colors[1])
ax.barh(y_pos[1]-(wd*2), z_perc_deaths_averted_dict[1][1,1], color=colors[2])
ax.barh(y_pos[2]-(wd*0), z_perc_deaths_averted_dict[2][0,1], color=colors[0])
ax.barh(y_pos[2]-(wd*1), z_perc_deaths_averted_dict[2][1,0], color=colors[1])
ax.barh(y_pos[2]-(wd*2), z_perc_deaths_averted_dict[2][1,1], color=colors[2])
ax.set_yticks(y_pos-(wd/3))
ax.set_yticklabels(['Boost 100% vaccinated,\nprime 10% unvaccinated', 'Boost 100% vaccinated,\nprime 50% unvaccinated', 'Boost 100% vaccinated,\nprime 100% unvaccinated'])
ax.set_xlabel('% of deaths averted')

ax.legend(title='Next generation vaccine')

ax = axv[2]
# Plot bar charts
ax.barh(y_pos[0]-(wd*0), z_doses_deaths_averted_dict[0][0,1], color=colors[0], label='Durable')
ax.barh(y_pos[0]-(wd*1), z_doses_deaths_averted_dict[0][1,0], color=colors[1], label='Broadly neutralizing')
ax.barh(y_pos[0]-(wd*2), z_doses_deaths_averted_dict[0][1,1], color=colors[2], label='Durable & broadly neutralizing')
ax.barh(y_pos[1]-(wd*0), z_doses_deaths_averted_dict[1][0,1], color=colors[0])
ax.barh(y_pos[1]-(wd*1), z_doses_deaths_averted_dict[1][1,0], color=colors[1])
ax.barh(y_pos[1]-(wd*2), z_doses_deaths_averted_dict[1][1,1], color=colors[2])
ax.barh(y_pos[2]-(wd*0), z_doses_deaths_averted_dict[2][0,1], color=colors[0])
ax.barh(y_pos[2]-(wd*1), z_doses_deaths_averted_dict[2][1,0], color=colors[1])
ax.barh(y_pos[2]-(wd*2), z_doses_deaths_averted_dict[2][1,1], color=colors[2])
ax.set_yticks(y_pos-(wd/3))
ax.set_yticklabels(['Boost 100% vaccinated,\nprime 10% unvaccinated', 'Boost 100% vaccinated,\nprime 50% unvaccinated', 'Boost 100% vaccinated,\nprime 100% unvaccinated'])
ax.set_xlabel('Doses per death averted')
ax.xaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))

fig.subplots_adjust(left=0.25)
fig.subplots_adjust(right=0.85)
fig.tight_layout(pad=3.0)
fig.show()
sc.savefig(str(sc.path(figdir) / 'next_gen_vx.png'), fig=fig)






print('Done')