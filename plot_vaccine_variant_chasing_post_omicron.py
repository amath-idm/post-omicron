
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



days_to_start = np.arange(-60, 100, 5)
vaccine_prime = [0, 0.5, 1]
vaccine_boost = [1]

# Load data
sc.options(dpi=100)
sc.fonts(add=sc.thisdir(aspath=True) / 'avenir')
sc.options(font='Avenir')

# Load data
resfolder = './results'
df1 = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-variant-chasing_vx_rollout_data.obj')))
df2 = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-variant-chasing_data_for_run_plot.obj')))
df = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-cov-sweep_vx_rollout_data.obj')))
df = pd.DataFrame.from_dict(df)
df1 = pd.DataFrame.from_dict(df1)

df_nextgen = sc.loadobj((str(sc.path(resfolder) / 'post-omicron-next-gen-vaccines_vx_rollout_data.obj')))
df_nextgen = pd.DataFrame.from_dict(df_nextgen)
figdir = './figs'

#%% Plotting

n_reps = 50
n_days = len(days_to_start)
total = n_reps * n_days


df_nextgen_novax = df_nextgen[df_nextgen['vaccine_breadth'] + df_nextgen['vaccine_durability'] == 0]
deaths_no_vax = np.mean(df_nextgen_novax['deaths'].values)

df = df[df['vaccine_prime'] + df['vaccine_boost']==2]
df = df[df['next_variant'] == 'New cluster, more severe']
df = df[df['new_variant_day'] == '2022-08-25']
df_deaths_pfizer = np.mean(df['deaths'])

df1 = df1[df1['vaccine_prime'] + df1['vaccine_boost']==2]
# deaths_averted = df_deaths_pfizer - df1['deaths']
# perc_deaths_averted = deaths_averted/df_deaths_pfizer
deaths_averted = deaths_no_vax - df1['deaths']
perc_deaths_averted = deaths_averted/deaths_no_vax
df1['Deaths averted'] = deaths_averted
df1['% Deaths averted'] = perc_deaths_averted

df_nextgen = df_nextgen[df_nextgen['vaccine_breadth']>0]
# deaths_averted = df_deaths_pfizer - df_nextgen['deaths']
# perc_deaths_averted = deaths_averted/df_deaths_pfizer
deaths_averted = deaths_no_vax - df_nextgen['deaths']
perc_deaths_averted = deaths_averted/deaths_no_vax
df_nextgen['Deaths averted'] = deaths_averted
df_nextgen['% Deaths averted'] = perc_deaths_averted
df_nextgen['Doses per death averted'] = df_nextgen['doses']/df_nextgen['Deaths averted']
vaccine_strategy = []
for row in range(len(df_nextgen)):
    if df_nextgen.iloc[row,0]:
        vaccine = 'Broadly neutralizing'
        if df_nextgen.iloc[row, 1]:
            vaccine += ' & durable'
    else:
        vaccine = ''
        if df_nextgen.iloc[row, 1]:
            vaccine += ' Durable'

    vaccine_strategy.append(vaccine)
df_nextgen['Coverage'] = vaccine_strategy
# df_nextgen = df_nextgen[df_nextgen['Deaths averted'] > 0]

df_to_use = df1[df1['vaccine_prime'] + df1['vaccine_boost']==2]
# df_to_use = df1[df1['vaccine_prime'] + df1['vaccine_boost']>0]
# df_to_use = df_to_use[df_to_use['Deaths averted'] > 0]
df_to_use['Doses per death averted'] = df_to_use['doses']/df_to_use['Deaths averted']
df_to_use_doses = df_to_use[df_to_use['Doses per death averted'] > 0]

vaccine_strategy = []
for row in range(len(df_to_use)):
    vaccine_strategy.append('Variant-chasing vaccine')
df_to_use['Coverage'] = vaccine_strategy

plt.rcParams['font.size'] = 20
fig, axv = plt.subplots(2, 1, figsize=(12, 16))
sns.lineplot(data=df_to_use.reset_index(), x='days_to_start', y='% Deaths averted', hue='Coverage',
ax=axv[0], palette='flare')
sns.lineplot(data=df_nextgen.reset_index(), x='days_to_start', y='% Deaths averted', hue='Coverage',
ax=axv[0], palette='crest')
axv[0].yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
axv[0].grid(alpha=0.3)
axv[0].legend(title='Vaccine')
axv[0].set_title('Deaths averted (%) compared to no additional vaccination')
axv[0].set_xlabel('')
sns.lineplot(data=df_to_use.reset_index(), x='days_to_start', y='Doses per death averted', hue='Coverage',
ax=axv[1], palette='flare')
sns.lineplot(data=df_nextgen.reset_index(), x='days_to_start', y='Doses per death averted', hue='Coverage',
ax=axv[1], palette='crest')
axv[1].grid(alpha=0.3)
axv[1].set(yscale='log')
axv[1].set_title('Doses per death averted compared to no additional vaccination')
# axv[1].set_ylim(bottom=0, top=15000)
axv[1].set_xlabel('Days from variant introduction to vaccine delivery')
axv[1].get_legend().remove()
fig.show()
# ax.set_ylim(bottom=0, top=2000000)
sc.savefig(str(sc.path(figdir) / 'variant_chasing.png'), fig=fig)
print('Done')