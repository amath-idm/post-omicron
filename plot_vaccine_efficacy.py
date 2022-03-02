
'''
Script to plot vaccine efficacy as a function of time
'''

import sciris as sc
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

# Get doses per death averted ready
window = 60
dates = inf_df['Date']
doses_per_death_averted = []
for i, day in enumerate(res['vx_day']):
    deaths_df[i] = deaths_df[i].reset_index()
    day_ind = deaths_df[i].index[deaths_df[i]['Date'] == day]
    end_of_window_ind = day_ind + window
    vx_deaths = res['cum_deaths'][i]
    vx_doses = res['cum_doses'][i]
    no_vx_deaths = deaths_df[i].iloc[end_of_window_ind,1].values - deaths_df[i].iloc[day_ind,1].values
    deaths_averted = no_vx_deaths - vx_deaths
    doses_per_death_averted.append(vx_deaths/deaths_averted)

res['Doses per death averted'] = doses_per_death_averted


fig, axv = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

# FIRST AXIS
ax = axv[0]
sns.lineplot(data=inf_df.reset_index(), x='Date', y='Infections', hue='Variant', ci='sd', 
ax=ax, palette='flare')
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.axvspan(res['vx_day'][0], res['vx_day'][60], alpha=0.5, color='goldenrod')
ax.get_legend().remove()
ax.grid(alpha=0.3)
ax.set_title('Infections by variant and percent exposed')
ax.set_ylim(bottom=0, top=2000000)

# TWIN FIRST AXIS
ax = axv[0].twinx()
sns.lineplot(data=exp_df.reset_index(), x='Date', y='Exposed (%)', ci='sd', color='k', ls='--', palette='tab10',
                 lw=2, ax=ax)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
# ax.grid()
# ax.set_title('Exposed (%)')
ax.set_ylim(bottom=0, top=100)

# SECOND AXIS
ax = axv[1]
sns.lineplot(data=sev_df.reset_index(), x='Date', y='Severe', hue='Variant', ci='sd', ax=ax, palette='flare')
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.axvspan(res['vx_day'][0], res['vx_day'][60], alpha=0.5, color='goldenrod')
# ax.text(res['vx_day'][0], 100000, 'Example efficacy \nwindow')
ax.annotate('Example efficacy window', xy=(0.13, -0.05), xytext=(0.13, -0.15), xycoords='axes fraction',
            fontsize=8, ha='center', va='bottom',
            bbox=dict(boxstyle='square', fc='white'),
            arrowprops=dict(arrowstyle='-[, widthB=5.0, lengthB=1', lw=1.0))
ax.grid(alpha=0.3)
ax.set_title('Severe cases by variant')
ax.set_ylim(bottom=0, top=200000)

# LAST AXIS
ax = axv[-1]
sns.lineplot(data=res, x='vx_day', y='VE_inf', ax=ax, lw=2, label='Infection')
sns.lineplot(data=res, x='vx_day', y='VE_symp', ax=ax, lw=2, label='Symptomatic disease')
sns.lineplot(data=res, x='vx_day', y='VE_sev', ax=ax, lw=2, label='Severe disease')
ax.set_xlabel('Date')
ax.set_ylabel('Vaccine efficacy (60 day window)')
ax.grid(alpha=0.3)
ax.set_title('Vaccine efficacy')

sc.dateformatter()
fig.show()
sc.savefig(str(sc.path(figdir) / 'vaccine_efficacy.png'), fig=fig)
print('Done.')