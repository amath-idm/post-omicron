
'''
Script to plot vaccine efficacy as a function of time
'''

import sciris as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.dates as mdates
import matplotlib.ticker as mtick

# Load data
sc.fonts(add=sc.thisdir(aspath=True) / 'avenir')

# Load data
resfolder = './results'

inf_df = sc.loadobj((str(sc.path(resfolder) / 'inf.obj')))
sev_df = sc.loadobj((str(sc.path(resfolder) / 'sev.obj')))
exp_df = sc.loadobj((str(sc.path(resfolder) / 'exp.obj')))
deaths_df = sc.loadobj((str(sc.path(resfolder) / 'deaths.obj')))
res = sc.loadobj((str(sc.path(resfolder) / 'res.obj')))
figdir = './figs'

death_df = pd.concat(deaths_df)

# Get doses per death averted ready
n_reps = len(deaths_df)
dates = np.unique(res['vx_day'].values)
n_dates = len(dates)
doses_per_death_averted = []
nnv_tot = []

for rep in range(n_reps):
    deaths_df[rep] = deaths_df[rep].reset_index()
    for i, day in enumerate(dates):
        day_ind = deaths_df[rep].index[deaths_df[rep]['Date'] == day]
        index = i + n_dates*rep
        vx_deaths = res['cum_deaths'][index]
        vx_doses = res['cum_doses'][index]
        no_vx_deaths = deaths_df[rep].iloc[-1]['Cumulative Deaths'] - deaths_df[rep].iloc[day_ind]['Cumulative Deaths'].values[0]
        deaths_averted = no_vx_deaths - vx_deaths
        doses_per = vx_doses/deaths_averted
        nnv = vx_doses/2/deaths_averted
        if np.isinf(doses_per):
            doses_per = 0
            nnv = 0
        if doses_per < 0:
            doses_per = 0
            nnv = 0
        doses_per_death_averted.append(doses_per)
        nnv_tot.append(nnv)

res['Doses per death averted'] = doses_per_death_averted
res['NNV'] = nnv_tot


fig, axv = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

# FIRST AXIS
ax = axv[0]
sns.lineplot(data=inf_df.reset_index(), x='Date', y='Infections', hue='Variant', ci='sd', 
ax=ax, palette='flare')
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.axvspan(inf_df.iloc[100]['Date'], inf_df.iloc[180]['Date'], alpha=0.5, color='goldenrod')
ax.grid(alpha=0.3)
ax.set_title('Infections by variant and percent exposed (no vaccination)')
ax.set_ylim(bottom=0, top=2000000)
ax.legend(bbox_to_anchor=(0.15, 1), title='Variant')

# TWIN FIRST AXIS
ax = axv[0].twinx()
sns.lineplot(data=exp_df.reset_index(), x='Date', y='Exposed (%)', ci='sd', color='green', ls='--', palette='tab10',
                 lw=2, ax=ax)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.PercentFormatter())
# ax.grid()
ax.set_ylabel('Exposed (%)', color='green')
ax.set_ylim(bottom=0, top=100)

# SECOND AXIS
ax = axv[1]
ax.grid(alpha=0.3)
sns.lineplot(data=sev_df.reset_index(), x='Date', y='Severe', hue='Variant', ci='sd', ax=ax, palette='flare')
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.axvspan(inf_df.iloc[100]['Date'], inf_df.iloc[180]['Date'], alpha=0.5, color='goldenrod')
ax.get_legend().remove()
# ax.text(res['vx_day'][0], 100000, 'Example efficacy \nwindow')
ax.annotate('Example efficacy window', xy=(0.083, -0.05), xytext=(0.083, -0.15), xycoords='axes fraction',
            fontsize=8, ha='center', va='bottom',
            bbox=dict(boxstyle='square', fc='white'),
            arrowprops=dict(arrowstyle='-[, widthB=0.75, lengthB=.8', lw=1.0))
ax.set_title('Severe cases by variant and cumulative deaths (no vaccination)')
ax.set_ylim(bottom=0, top=150000)

# TWIN SECOND AXIS
ax = axv[1].twinx()
sns.lineplot(data=death_df.reset_index(), x='Date', y='Cumulative Deaths', ci='sd', color='blue', ls='--', palette='tab10',
                 lw=2, ax=ax)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
# ax.grid()
ax.set_ylabel('Cumulative deaths', color='blue')
# ax.set_ylim(bottom=0, top=100)

# LAST AXIS
ax = axv[-1]
# sns.lineplot(data=res, x='vx_day', y='VE_inf', ax=ax, lw=2, label='Infection')
# sns.lineplot(data=res, x='vx_day', y='VE_symp', ax=ax, lw=2, label='Symptomatic disease')
sns.lineplot(data=res, x='vx_day', y='VE_sev', ax=ax, lw=2, label='Severe disease')
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax.set_xlabel('Date')
ax.set_ylabel('Vaccine efficacy (60 day window)')
ax.grid(alpha=0.3)
ax.set_title('Vaccine efficacy and efficiency')
ax.set_ylim(bottom=0, top=1)
ax.legend()

# TWIN LAST AXIS
ax = axv[-1].twinx()
dose_per = sns.lineplot(data=res, x='vx_day', y='Doses per death averted', color='purple', ls='--', palette='tab10',
                 lw=2, ax=ax)
dose_per.set(yscale='log')
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.set_ylabel('Doses per death averted', color='purple')
# ax.grid()
# ax.set_ylim(bottom=0, top=10000)

sc.dateformatter()
fig.subplots_adjust(right=0.88)
fig.subplots_adjust(left=0.15)
# fig.show()
sc.savefig(str(sc.path(figdir) / 'vaccine_efficacy.png'), fig=fig)
print('Done.')

fig, axv = plt.subplots(3, 1, figsize=(8, 8), sharex=True)

# FIRST AXIS
ax = axv[0]
sns.lineplot(data=inf_df.reset_index(), x='Date', y='Infections', hue='Variant', ci='sd', 
ax=ax, palette='flare')
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.axvspan(inf_df.iloc[100]['Date'], inf_df.iloc[180]['Date'], alpha=0.5, color='goldenrod')
ax.grid(alpha=0.3)
ax.set_title('Infections by variant and percent exposed (no vaccination)')
ax.set_ylim(bottom=0, top=2000000)
ax.legend(bbox_to_anchor=(0.15, 1), title='Variant')
ax.annotate('Example efficacy window', xy=(0.083, -0.05), xytext=(0.083, -0.15), xycoords='axes fraction',
            fontsize=8, ha='center', va='bottom',
            bbox=dict(boxstyle='square', fc='white'),
            arrowprops=dict(arrowstyle='-[, widthB=0.75, lengthB=.8', lw=1.0))



# TWIN FIRST AXIS
ax = axv[0].twinx()
sns.lineplot(data=exp_df.reset_index(), x='Date', y='Exposed (%)', ci='sd', color='k', ls='--', palette='tab10',
                 lw=2, ax=ax)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.PercentFormatter())
# ax.grid()
# ax.set_title('Exposed (%)')
ax.set_ylim(bottom=0, top=100)

# SECOND AXIS
ax = axv[1]
sns.lineplot(data=res, x='vx_day', y='VE_inf', ax=ax, lw=2, label='Infection')
sns.lineplot(data=res, x='vx_day', y='VE_symp', ax=ax, lw=2, label='Symptomatic disease')
sns.lineplot(data=res, x='vx_day', y='VE_sev', ax=ax, lw=2, label='Severe disease')
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax.set_xlabel('Date')
ax.set_ylabel('Vaccine efficacy\n(60 day window)')
ax.grid(alpha=0.3)
ax.set_title('Vaccine efficacy')
ax.set_ylim(bottom=0, top=1)
ax.legend()


# LAST AXIS
ax = axv[-1]
dose_per = sns.lineplot(data=res, x='vx_day', y='Doses per death averted', color='k', palette='tab10',
                 lw=2, ax=ax)
dose_per.set(yscale='log')
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.grid(alpha=0.3)
ax.set_title('Vaccine efficiency')
ax.set_ylim(bottom=0, top=10000)
ax.set_xlabel('Date')
sc.dateformatter()
fig.subplots_adjust(right=0.88)
fig.subplots_adjust(left=0.15)
# fig.show()
sc.savefig(str(sc.path(figdir) / 'vaccine_efficacy_v2.png'), fig=fig)
print('Done.')