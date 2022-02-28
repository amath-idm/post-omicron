
'''
Script to plot vaccine efficacy as a function of time
'''

import sciris as sc
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.dates as mdates
import matplotlib.ticker as mtick

# Load data
sc.fonts(add='C:\\Users\\jamieco\\PycharmProjects\\avenir')

# Load data
resfolder = './results'

inf_df = sc.loadobj((str(sc.path(resfolder) / 'inf.obj')))
sev_df = sc.loadobj((str(sc.path(resfolder) / 'sev.obj')))
exp_df = sc.loadobj((str(sc.path(resfolder) / 'exp.obj')))
res = sc.loadobj((str(sc.path(resfolder) / 'res.obj')))
figdir = './figs'

fig, axv = plt.subplots(3, 1, figsize=(8, 10), sharex=True)

# FIRST AXIS
ax = axv[0]
sns.lineplot(data=inf_df.reset_index(), x='Date', y='Infections', hue='Variant', ci='sd', 
ax=ax, palette='flare')
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
ax.axvline(x=res['vx_day'][0], color='goldenrod', ls='-.', label='Window start')
ax.axvline(x=res['vx_day'][60], color='goldenrod', ls='-.', label='Window end')
ax.text(res['vx_day'][0], 0.3, 'Example window start', rotation=90, transform=ax.get_xaxis_text1_transform(0)[0])
ax.text(res['vx_day'][60], 0.3, 'Example window end', rotation=90, transform=ax.get_xaxis_text1_transform(0)[0])
ax.get_legend().remove()
ax.grid()
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
ax.axvline(x=res['vx_day'][0], color='goldenrod', ls='-.', label='Window start')
ax.axvline(x=res['vx_day'][60], color='goldenrod', ls='-.', label='Window end')
ax.text(res['vx_day'][0], 0.3, 'Example window start', rotation=90, transform=ax.get_xaxis_text1_transform(0)[0])
ax.text(res['vx_day'][60], 0.3, 'Example window end', rotation=90, transform=ax.get_xaxis_text1_transform(0)[0])
ax.grid()
ax.set_title('Severe cases by variant')
ax.set_ylim(bottom=0, top=200000)

# LAST AXIS
ax = axv[-1]
sns.lineplot(data=res, x='vx_day', y='VE_inf', ax=ax, lw=2, label='Infection')
sns.lineplot(data=res, x='vx_day', y='VE_symp', ax=ax, lw=2, label='Symptomatic disease')
sns.lineplot(data=res, x='vx_day', y='VE_sev', ax=ax, lw=2, label='Severe disease')
ax.axvline(x=res['vx_day'][0], color='goldenrod', ls='-.', label='Window start')
ax.axvline(x=res['vx_day'][60], color='goldenrod', ls='-.', label='Window end')
ax.text(res['vx_day'][0], 0.3, 'Example window start', rotation=90, transform=ax.get_xaxis_text1_transform(0)[0])
ax.text(res['vx_day'][60], 0.3, 'Example window end', rotation=90, transform=ax.get_xaxis_text1_transform(0)[0])
ax.set_xlabel('Date')
ax.set_ylabel('Vaccine efficacy (60 day window)')
ax.grid()
ax.set_title('Vaccine efficacy (if vaccinating on this date)')

# fig.show()
sc.savefig(str(sc.path(figdir) / 'vaccine_efficacy.png'), fig=fig)
print('Done.')