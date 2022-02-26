
'''

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

df = sc.loadobj((str(sc.path(resfolder) / 'inf.obj')))
exp_df = sc.loadobj((str(sc.path(resfolder) / 'exp.obj')))
res = sc.loadobj((str(sc.path(resfolder) / 'res.obj')))
figdir = './figs'

fig, axv = plt.subplots(2, 1, figsize=(7, 6), sharex=True)
# FIRST AXIS
ax = axv[0]
sns.lineplot(data=df.reset_index(), x='Date', y='Infections', hue='Variant', ci='sd', ax=ax)
ax.xaxis.set_major_formatter(mdates.ConciseDateFormatter(ax.xaxis.get_major_locator()))
ax.yaxis.set_major_formatter(mtick.FuncFormatter(lambda x, p: format(int(x), ',')))
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

# LAST AXIS
ax = axv[-1]
sns.lineplot(data=res, x='vx_day', y='VE_inf', ax=ax, lw=2, label='Infection')
sns.lineplot(data=res, x='vx_day', y='VE_symp', ax=ax, lw=2, label='Symptomatic disease')
sns.lineplot(data=res, x='vx_day', y='VE_sev', ax=ax, lw=2, label='Severe disease')
ax.set_xlabel('Date')
ax.set_ylabel('Vaccine efficacy (60 day window)')
ax.grid()
ax.set_title('Vaccine efficacy (if vaccinating on this date)')

fig.show()
sc.savefig(str(sc.path(figdir) / 'vaccine_efficacy.png'), fig=fig)
print('Done.')