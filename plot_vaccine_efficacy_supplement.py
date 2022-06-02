
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

res = sc.loadobj((str(sc.path(resfolder) / 'res.obj')))
figdir = './figs'

fig, ax = plt.subplots(figsize=(8, 8), sharex=True)

sns.lineplot(data=res, x='vx_day', y='VE_sev', ax=ax, lw=2, label='Severe disease')
sns.lineplot(data=res, x='vx_day', y='VE_symp', ax=ax, lw=2, label='Symptomatic disease')
sns.lineplot(data=res, x='vx_day', y='VE_inf', ax=ax, lw=2, label='Infection')
ax.yaxis.set_major_formatter(mtick.PercentFormatter(1.0))
ax.set_xlabel('Date')
ax.set_ylabel('Vaccine effectiveness (60 day window)')
ax.grid(alpha=0.3)
ax.set_title('Vaccine effectiveness')
ax.set_ylim(bottom=0, top=1)
ax.legend()
sc.dateformatter()
fig.show()
sc.savefig(str(sc.path(figdir) / 'vaccine_efficacy_supp.png'), fig=fig)
print('Done.')

