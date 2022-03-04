
'''
Script to plot variant-chasing vaccine impact
'''

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.dates as mdates
import matplotlib.ticker as mtick



days_to_start = np.arange(5, 100, 5)
vaccine_prime = [0, 0.5, 1]
vaccine_boost = [0, 0.5, 1]

# Load data
sc.options(dpi=100)
sc.fonts(add=sc.thisdir(aspath=True) / 'avenir')
sc.options(font='Avenir')

# Load data
resfolder = './results'
df1 = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-variant-chasing_vx_rollout_data.obj')))
df2 = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-variant-chasing_data_for_run_plot.obj')))
df1 = pd.DataFrame.from_dict(df1)
figdir = './figs'

#%% Plotting

n_reps = 50
n_days = len(days_to_start)
total = n_reps * n_days

deaths_no_vx = (df1['deaths'].iloc[0:950]).values
deaths_averted = [0]*total
for i in range(1,9):
    deaths_averted += list(deaths_no_vx - (df1['deaths'].iloc[total*i:total*i+total]).values)

df1['Deaths averted'] = deaths_averted
df1['% Deaths averted'] = (deaths_averted/df1['deaths'])*100

df_to_use = df1[df1['vaccine_prime'] + df1['vaccine_boost']>0]
df_to_use = df_to_use[df_to_use['vaccine_prime'] == 1]
df_to_use = df_to_use[df_to_use['Deaths averted'] > 0]
df_to_use['Doses per death averted'] = df_to_use['doses']/df_to_use['Deaths averted']
df_to_use = df_to_use[df_to_use['Doses per death averted'] > 0]
vaccine_strategy = []
for row in range(len(df_to_use)):
    vaccine_strategy.append(f'Boost {int(df_to_use.iloc[row,1]*100)}% vaccinated, prime {int(df_to_use.iloc[row,0]*100)}% unvaccinated')
df_to_use['Coverage'] = vaccine_strategy


fig, axv = plt.subplots(3, 1, figsize=(8, 10))
ax = axv[0]

new_inf = df2['new_infections_by_variant'][0]
x = df2['datevec_sim'][333:433]
ax.plot(x, new_inf.values[4,333:433], label='Infections')
ax.fill_between(x, (new_inf.low[4,333:433]), (new_inf.high[4,333:433]), alpha=0.3)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.set_title('Next variant infections and deaths (no additional vaccination)')
ax.grid(alpha=0.3)
ax.set_ylabel('New infections')

ax = axv[0].twinx()
new_deaths = df2['new_deaths'][0]
x = df2['datevec_sim'][333:433]
ax.plot(x, new_deaths.values[333:433], label='Deaths', color='green', ls='--')
ax.fill_between(x, (new_deaths.low[333:433]), (new_deaths.high[333:433]), alpha=0.3, color='green')
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.set_ylabel('New deaths', color='green')

ax = axv[1]
sns.lineplot(data=df_to_use.reset_index(), x='days_to_start', y='% Deaths averted', hue='Coverage',
ax=ax, palette='flare')
# ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.grid(alpha=0.3)
ax.legend(title='Coverage')
ax.set_title('% deaths averted')
ax.set_xticks([])
ax.set_xticklabels([])
ax.set_xlabel('')

ax = axv[-1]
l1 = sns.lineplot(data=df_to_use.reset_index(), x='days_to_start', y='Doses per death averted', hue='Coverage',
ax=ax, palette='flare')
# l1.set(yscale='log')
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.grid(alpha=0.3)
ax.get_legend().remove()
ax.set_xlabel('Days from variant introduction to vaccine delivery')
ax.set_title('Doses per death averted')

fig.show()
# ax.set_ylim(bottom=0, top=2000000)
sc.savefig(str(sc.path(figdir) / 'variant_chasing.png'), fig=fig)
print('Done')