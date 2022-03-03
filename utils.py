'''

Utils for the Omicron Endgame analysis
'''

import numpy as np
import sciris as sc
import pandas as pd
import covasim as cv
import covasim.parameters as cvpar
import os
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from interventions import treatment

module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

cv.check_version('>=3.1.1', die=True)
verbose = 0
debug = 0
plot_uncertainty = True
figfolder = 'figs'


def sweep_params(sample=True):
    p = sc.objdict(
        rel_beta_omicron = [2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6],
        rel_imm = [2**n/100 for n in range(7)], #[1/50, 1 / 25, 1 / 10, 1 / 5],
        nab_decay = ['both_fast',
            # 'both_slow'
        ])
    if sample:
        return sc.objdict({pname: np.random.choice(pval) for pname, pval in par_ranges.items()}), 1
    else:
        dims = [len(v) for k, v in p.items()]
        npars = np.prod(dims)
        return p, npars



def make_vx_intv(vaccine, interval=90, day=None, coverage=1):
    '''
    Return an age-targeted booster vaccine intervention
        * vaccine: the name of the vaccine, e.g. "pfizer"
        * day: day on which to vaccinate
        * coverage: vaccine coverage
    '''

    day = np.arange(cv.day('2022-02-15', start_date='2021-10-01'), cv.day('2022-12-15', start_date='2021-10-01')) if \
        day is None else day
    subtarget = {
        'inds': lambda sim: cv.true((sim.people.doses == 2) & ((sim.t - sim.people.date_vaccinated) >= interval)),
        'vals': coverage}
    if vaccine == 'Omicron 1-dose':
        pfizer = cvpar.get_vaccine_dose_pars(vaccine='pfizer')
        pfizer['nab_init']['par1'] = 2
        intv = cv.vaccinate_prob(pfizer, days=day, label='omicron-1-dose', prob=0, subtarget=subtarget, booster=True,
                                 do_plot=False)
    elif vaccine == 'Omicron 2-dose':
        pfizer = cvpar.get_vaccine_dose_pars(vaccine='pfizer')
        intv = cv.vaccinate_prob(pfizer, days=day, label='omicron-2-dose', prob=0, subtarget=subtarget, booster=False,
                                 do_plot=False)
    else:
        intv = cv.vaccinate_prob(vaccine=vaccine, days=day, prob=0, subtarget=subtarget, booster=True,
                                 do_plot=False)

    return sc.promotetolist(intv)


def make_tx_intv(start_day, rel_severe=0.11, coverage=1, age_min=None, age_max=None, label="Paxlovid"):
    '''
    Return an age-targeted treatment intervention
        * rel_severe: change in severe probability
        * start_day: day to start treatment
        * coverage: treatment coverage
        * age_min: the minimum age (inclusive), or None for no minimum
        * age_max: the maximum age (exclusive), or None for no maximum
    '''
    age_min = 0 if age_min is None else age_min
    age_max = 1000 if age_max is None else age_max
    subtarget = {'inds': lambda sim: cv.true((sim.people.age >= age_min) & (sim.people.age < age_max) & (~sim.people.severe)
                                         & (~sim.people.critical) & (sim.people.date_symptomatic == sim.t)), 'vals': coverage}
    label = f'{label}: {age_min}-{age_max} @ {coverage}'
    intv = treatment(rel_severe=rel_severe, start_day=start_day, prob=0, subtarget=subtarget, eligible='symptomatic', do_plot=False, label=label)
    return intv


def get_fit_params(index):
    params, nparams = sweep_params(sample=False)
    dims = [len(v) for k, v in params.items()]
    indices = index % nparams
    inds = np.unravel_index(indices, dims)
    f = sc.objdict()
    for ik, (k, v) in enumerate(params.items()):
        f.update({k: v[inds[ik]]})
    return f



def plotts(msims, filename='time_series.png'):
    if verbose >= 0: print('Plotting timeseries')
    # Just plot one

    labels = {'new_infections_by_variant': 'Number of new infections by variant',
              'new_severe_by_variant': 'Number of new severe cases by variant',
              }

    for which in ['new_infections_by_variant', 'new_severe_by_variant']:
        Y = dict()
        X = dict()
        for j, s in enumerate(msims):
            variant_map = s.base_sim['variant_map']
            hist_wave = s.base_sim.get_intervention(cv.historical_wave)
            X[j] = np.array(list(hist_wave.__getattribute__('results')['datevec']))
            Y[j] = hist_wave.__getattribute__('results')['variant'][which].values

        colors = sc.gridcolors(len(df.imm_source))
        fig, ax = plt.subplots()
        plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
        plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=60))
        for j, data in Y.items():
            for i, var_data in enumerate(data):
                if i != 3:
                    x = X[i]
                    if j == 0:
                        ax.plot(x, var_data, color=colors[i], label=variant_map[i])
                    else:
                        ax.plot(x, var_data, color=colors[i])
        plt.gcf().autofmt_xdate()
        ax.legend()
        ax.set_title(labels[which])
        fig.savefig(f'{figfolder}/{filename}')
    if verbose >= 0: print('loaded')


def results_after(sim, date, channel):
    '''
    Sum results for a specified channel after a specified date

    date (str): Like '2021-11-15'
    channel (str): sim.results channel like "new_deaths"
    '''
    idx = np.argmax(sim.results['date'] == sc.date(date))
    return np.sum(sim.results[channel][idx:])


def get_data_for_plots(msims, filename='data_for_plot.obj', historical_wave=True):

    df = sc.objdict()
    df.pop_nabs = []
    df.pop_nabs_high = []
    df.pop_nabs_low = []
    if historical_wave:
        df.pop_nabs_historic = []
        df.delta_peak_inf = []
        df.delta_peak_sev = []
    df.new_infections_by_variant = []
    df.new_severe_by_variant = []
    df.new_deaths = []
    df.datevecs = []
    params = msims[0].base_sim.run_info.keys()
    for param in params:
        df[param] = []
    for s in msims:  # per parameter set
        for param in params:
            df[param].append(s.base_sim.run_info[param])
        df.new_infections_by_variant.append(s.results['variant']['new_infections_by_variant'])
        df.new_severe_by_variant.append(s.results['variant']['new_severe_by_variant'])
        df.new_deaths.append(s.results['new_deaths'])
        datevec_sim = s.results['date']
        if historical_wave:
            hist_wave = s.base_sim.get_intervention(cv.historical_wave)
            hist_vax = s.base_sim.get_intervention(cv.historical_vaccinate_prob)
            df.delta_peak_inf.append(np.max(hist_wave.results['variant']['new_infections_by_variant'].values[2, :]))
            df.delta_peak_sev.append(np.max(hist_wave.results['variant']['new_severe_by_variant'].values[2, :]))
            new_date_vec = list(max(hist_vax.results['datevec'], hist_wave.results['datevec'], key=len)) + list(
            s.results['date'])
            df.datevecs.append(np.array(new_date_vec))
            additional_values = hist_wave.results['pop_nabs_by_source'].values
            x = abs(len(hist_vax.results['datevec']) - len(hist_wave.results['datevec']))
            additional_values[:, x:] += hist_vax.results['pop_nabs_by_source'].values
            df.pop_nabs_historic.append(additional_values)
            pop_nab = np.hstack((additional_values, s.results['pop_nabs_by_source'].values))
            pop_nab_high = np.hstack((additional_values, s.results['pop_nabs_by_source'].high))
            pop_nab_low = np.hstack((additional_values, s.results['pop_nabs_by_source'].low))
        else:
            pop_nab = s.results['pop_nabs_by_source'].values
            pop_nab_high = s.results['pop_nabs_by_source'].high
            pop_nab_low = s.results['pop_nabs_by_source'].low
        
        df.pop_nabs.append(pop_nab)
        df.pop_nabs_high.append(pop_nab_high)
        df.pop_nabs_low.append(pop_nab_low)
        base = s.base_sim
        imm_source = base.pars['variant_map']
        n_var = base.pars['n_variants']
        for i, vac in base.pars['vaccine_map'].items():
            imm_source[i + n_var] = vac

    df.imm_source = imm_source
    df.datevec_sim = datevec_sim
    sc.saveobj(filename, df)
    return df


def plot_the_rest(df):
    colors = sc.gridcolors(len(df.imm_source))
    fig, ax = plt.subplots()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=60))
    for i, imm in df.imm_source.items():
        for j, nabs in enumerate(df.pop_nabs):
            if i < nabs.shape[0]:
                if j == len(df.pop_nabs) - 1:
                    ax.plot(df.datevecs[j], nabs[i, :], color=colors[i], label=imm)
                    if plot_uncertainty:
                        ax.fill_between(df.datevecs[j], df.pop_nabs_low[j][i, :], df.pop_nabs_high[j][i, :],
                                        color=colors[i], alpha=0.4)
                else:
                    ax.plot(df.datevecs[j], nabs[i, :], color=colors[i])
                    if plot_uncertainty:
                        ax.fill_between(df.datevecs[j], df.pop_nabs_low[j][i, :], df.pop_nabs_high[j][i, :],
                                        color=colors[i], alpha=0.4)
    ax.set_ylabel('Population NAbs')
    plt.gcf().autofmt_xdate()
    ax.legend()
    fig.savefig(f'{figfolder}/pop_nabs.png')

    colors = sc.gridcolors(len(df.new_infections_by_variant))
    fig, ax = plt.subplots()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    interval = 15

    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=interval))
    x = df.datevec_sim
    for i, res in enumerate(df.new_infections_by_variant):
        ax.plot(x, (res.values[3:, ][0] / df.delta_peak_inf[i]) * 100, color=colors[i], label=df.labels[i]['nab_decay'])
        if plot_uncertainty:
            ax.fill_between(x, (res.low[3:, ][0] / df.delta_peak_inf[i]) * 100,
                            (res.high[3:, ][0] / df.delta_peak_inf[i]) * 100, color=colors[i], alpha=0.4)
    ax.set_title('Omicron infections (% of delta peak)')
    plt.gcf().autofmt_xdate()
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    # Put a legend to the right of the current axis
    # ax.legend(loc='center left', fontsize=8, bbox_to_anchor=(1, 0.5))
    #fig.savefig(f'{args.root}_infections.png')
    fig.savefig(f'{figfolder}/infections.png')

    fig, ax = plt.subplots()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=interval))
    x = df.datevec_sim
    for i, res in enumerate(df.new_severe_by_variant):
        ax.plot(x, (res.values[3:, ][0] / df.delta_peak_sev[i]) * 100, color=colors[i], label=df.labels[i]['nab_decay'])
        if plot_uncertainty:
            ax.fill_between(x, (res.low[3:, ][0] / df.delta_peak_sev[i]) * 100,
                            (res.high[3:, ][0] / df.delta_peak_sev[i]) * 100, color=colors[i], alpha=0.4)
    ax.set_title('Omicron severe (% of delta peak)')
    plt.gcf().autofmt_xdate()
    # box = ax.get_position()
    # ax.set_position([box.x0, box.y0, box.width * 0.6, box.height])
    # Put a legend to the right of the current axis
    # ax.legend(loc='center left', fontsize=8, bbox_to_anchor=(1, 0.5))
    #fig.savefig(f'{args.root}_severe.png')
    fig.savefig(f'{figfolder}/severe.png')

    fig, ax = plt.subplots()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    interval = 60
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=interval))
    x = df.datevec_sim[1:]
    for i, res in enumerate(df.new_infections_by_variant):
        for j, row in enumerate(res.values):
            if i == 0:
                ax.plot(x, row[1:], color=colors[j], label=df.imm_source[j])
            else:
                ax.plot(x, row[1:], color=colors[j])
            if plot_uncertainty:
                ax.fill_between(x, res.low[j, 1:], res.high[j, 1:], color=colors[j], alpha=0.4)
    ax.set_title('Infections by variant')
    plt.gcf().autofmt_xdate()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', fontsize=8, bbox_to_anchor=(1, 0.5))
    #fig.savefig(f'{args.root}_infections_by_variant.png')
    fig.savefig(f'{figfolder}/infections_by_variant.png')

    fig, ax = plt.subplots()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
    interval = 60
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=interval))
    x = df.datevec_sim[1:]
    for i, res in enumerate(df.new_severe_by_variant):
        for j, row in enumerate(res.values):
            if i == 0:
                ax.plot(x, row[1:], color=colors[j], label=df.imm_source[j])
            else:
                ax.plot(x, row[1:], color=colors[j])
        if plot_uncertainty:
            ax.fill_between(x, res.low[j, 1:], res.high[j, 1:], color=colors[j], alpha=0.4)
    ax.set_title('Severe cases by variant')
    plt.gcf().autofmt_xdate()
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.85, box.height])
    # Put a legend to the right of the current axis
    ax.legend(loc='center left', fontsize=8, bbox_to_anchor=(1, 0.5))
#    fig.savefig(f'{args.root}_severe_by_variant.png')
    fig.savefig(f'{figfolder}/severe_by_variant.png')

    colors = sc.gridcolors(len(df.imm_source))
    fig, ax = plt.subplots()
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    plt.gca().xaxis.set_major_locator(mdates.DayLocator(interval=60))
    for i, imm in df.imm_source.items():
        for j, nabs in enumerate(df.pop_nabs_historic):
            if i != 3:
                if j == 0:
                    ax.plot(df.datevec_historic, nabs[i, :], color=colors[i], label=imm)
                else:
                    ax.plot(df.datevec_historic, nabs[i, :], color=colors[i])

    ax.set_ylabel('Population NAbs')
    plt.gcf().autofmt_xdate()
    ax.legend()
#    fig.savefig(f'{args.root}_pop_nabs_historic.png')
    fig.savefig(f'{figfolder}/pop_nabs_historic.png')
    return


