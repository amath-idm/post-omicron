'''

Script to show how vaccine efficacy varies as a function of time and immunity

on some day (to be swept), pick X% of [target cohort] to vaccinate and mark non-vaccinees as placebo.
Compute infections in each arm, relative risk / OR / ..., and ultimate VE over a period of Y days.
So it will be VE as a function of day on which randomization and vaccination occur

Simulate 4 discrete waves.

'''

import numpy as np
import sciris as sc
import pandas as pd
import covasim as cv
import os
import sys



module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

cv.check_version('>=3.1.1', die=True)
verbose = 0
debug = 0
plot_uncertainty = True
resfolder = 'results'
figdir = 'figs'


# Covasim default parameters will be overridden with the following
base_pars = sc.objdict(
    scaled_pop=54_000_000,  # Population size - does not seem to work?
    n_agents=54_000,  # The number of simulated agents
    beta=0.015,  # Base transmission probability per contact, per day
    pop_infected=10,  # Number of seed infections
    start_day='2020-03-01',  # First day of simulation
    end_day='2022-03-15',  # Last day of simulation
    interventions=[],  # Interventions to be added later
    analyzers=[],  # Analyzers to be added later
    use_waning=True,  # Enable waning immunity
    verbose=0,  # Turn off outputs to avoid notebook clutter
)

fast_decay = dict(form='nab_growth_decay', growth_time=21, decay_rate1=np.log(2) / 100, decay_time1=100,
                  decay_rate2=np.log(2) / 250, decay_time2=200)
slow_decay = dict(form='nab_growth_decay', growth_time=21, decay_rate1=0.006741981, decay_time1=47,
                  decay_rate2=0.001764455, decay_time2=106)

nab_decay_params = {
    'vax_fast_nat_slow': dict(natural=slow_decay, vaccine=fast_decay),
    'both_fast': dict(natural=fast_decay, vaccine=fast_decay),
    'nat_fast_vax_slow': dict(natural=fast_decay, vaccine=slow_decay),
    'both_slow': dict(natural=slow_decay, vaccine=slow_decay),
}


window = 80 # 60 days after 2nd dose
trial_size = 8000


# construct analyzer to select placebo and vaccine arms
class vaccine_trial_arms(cv.Analyzer):
    def __init__(self, day, trial_size, **kwargs):
        super().__init__(**kwargs)
        self.day = day
        self.trial_size = trial_size
        return

    def initialize(self, sim=None):
        self.placebo_inds = []
        self.vacc_inds = []
        self.initialized = True
        return

    def apply(self, sim):
        if sim.t == sim.day(self.day):
            placebo_eligible = cv.true(~sim.people.vaccinated)
            placebo_eligible = placebo_eligible[cv.true(sim.people.age[placebo_eligible]>= 12)]
            self.placebo_inds = placebo_eligible[cv.choose(len(placebo_eligible), min(self.trial_size, len(placebo_eligible)))]
            vacc_eligible = cv.true(sim.people.vaccinated)
            vacc_eligible = vacc_eligible[cv.true(sim.people.age[vacc_eligible]>= 12)]
            self.vacc_inds = vacc_eligible[cv.choose(len(vacc_eligible), min(self.trial_size, len(vacc_eligible)))]
        return


def make_vaccine_intv(day, vaccine='pfizer', coverage=0.35):
    '''
    Return a vaccine intervention on a single day
    vaccine (str): Vaccine to implement
    day (str): Day to implement vaccine
    coverage (float): Share of population to target with vaccine
    '''
    subtarget = {
        'inds': lambda sim: cv.true(sim.people.age > 12), 'vals': coverage}
    intv = cv.vaccinate_prob(vaccine=vaccine, days=day, prob=0, subtarget=subtarget, do_plot=False)

    return intv


def make_scenarios(n_reps=5, vx_res=5):
    scenarios = []

    for rand_seed in range(n_reps):
        scenarios += [{
            # No vaccination
            'label': 'No vaccination',
            'meta': {
                'pars': dict(base_pars, **{  # Update base_pars with scenario-specific parameters
                    'rand_seed': rand_seed,
                })
            },
        }]
        for vx_day_offset in np.linspace(1, sc.daydiff(base_pars['start_day'], base_pars['end_day'])-window, vx_res):
            vx_day = sc.datedelta(base_pars['start_day'], days=vx_day_offset)
            scenario = {
                'label': f'vx_day {vx_day}',
                'meta': {
                    'vx': {
                        'vaccine': 'pfizer',
                        'day': vx_day,
                    },
                    'pars': dict(base_pars, **{  # Update base_pars with scenario-specific parameters
                        'rand_seed': rand_seed,
                    })
                },
            }
            scenarios.append(scenario)
    return scenarios


def make_sim(label, meta):
    '''
    This "builder" function creates a single simulation. It is intended to be called in parallel

    vax_day (int): Day of vaccine intervention

    '''

    p = sc.dcp(meta['pars'])
    if 'vx' in meta:
        vx = meta['vx']
        p['interventions'] += [make_vaccine_intv(**vx)]
        vax_date = vx['day']
        vax_day = cv.day(vax_date, start_date=p['start_day'])
        p['analyzers'] += [vaccine_trial_arms(vax_date, trial_size=trial_size), cv.snapshot(days=vax_day + window)]

    # Create variants
    beta = cv.variant('beta', days='2020-10-15', n_imports=4000)
    beta.p['rel_beta'] = 1.4
    delta = cv.variant('delta', days='2021-05-01', n_imports=6000)
    omicron = cv.variant('omicron', days='2021-10-25', n_imports=4000)
    variants = [beta, delta, omicron]

    # Create beta interventions
    beta_interventions = [
        cv.change_beta('2020-06-20', 0.6),  # shut down
        cv.change_beta('2020-07-15', 0.3),  # shut down
        cv.change_beta('2020-10-15', 1),  # reopen
        cv.change_beta('2020-12-15', 0.6),  # shut down
        cv.change_beta('2021-01-01', 0.3),  # shut down
        cv.change_beta('2021-05-01', 1),  # shut down
        cv.change_beta('2021-07-01', 0.6),  # shut down
        cv.change_beta('2021-07-15', 0.3),  # shut down
        cv.change_beta('2021-10-25', 0.9),  # reopen
        cv.change_beta('2021-11-15', 0.6)   # shut down
    ]

    p['interventions'] += beta_interventions

    sim = cv.Sim(label=label, pars=p, variants=variants)
    sim.meta = meta

    return sim


def create_dfs(msim):
    exp_dfs = []
    dfs = []
    ret = []
    no_vax_res = []
    for sim in msim.sims:
        if 'vx' in sim.meta:
            vx_day = sim.day(sim.meta['vx']['day'])
            placebo_inds = sim['analyzers'][0].placebo_inds
            vacc_inds = sim['analyzers'][0].vacc_inds
            snap = sim['analyzers'][1].snapshots[0]
            VE_inf = 1 - (cv.true(snap.date_exposed[vacc_inds] > vx_day).sum() / cv.true(
                snap.date_exposed[placebo_inds] > vx_day).sum())
            if VE_inf < 0:
                VE_inf = 0
            VE_symp = 1 - (cv.true(snap.date_symptomatic[vacc_inds] > vx_day).sum() / cv.true(
                snap.date_symptomatic[placebo_inds] > vx_day).sum())
            if VE_symp < 0:
                VE_symp = 0
            VE_sev = 1 - (cv.true(snap.date_severe[vacc_inds] > vx_day).sum() / cv.true(
                snap.date_severe[placebo_inds] > vx_day).sum())
            if VE_sev < 0:
                VE_sev = 0
            ret.append({
                'label': sim.label,
                'VE_inf': VE_inf,
                'VE_symp': VE_symp,
                'VE_sev': VE_sev,
                'vx_day': sim.meta['vx']['day'],
            })
        else:
            # New infections by variant
            reskey = 'new_infections_by_variant'
            dat = sc.dcp(sim.results['variant'][reskey].values)
            d = pd.DataFrame(dat.T, index=pd.DatetimeIndex(sim.results['date'], name='Date'),
                             columns=sim['variant_map'].values())
            dfs.append(d)
            no_vax_res.append(sim.results)
            # Num exposed
            reskey = 'n_naive'
            dat = sc.dcp(sim.results[reskey].values)
            vals = dat.T

            d = pd.DataFrame(vals, index=pd.DatetimeIndex(sim.results['date'], name='Date'), columns=['Exposed'])
            d['Exposed (%)'] = 100 - 100 * d[
                'Exposed'] / sim.scaled_pop_size  # sim.results['cum_deaths'][:] - sim.results['n_recovered'][:] - sim.results['n_exposed'][:]
            d['rep'] = sim['rand_seed']
            exp_dfs.append(d)

    df = pd.concat(dfs).stack().reset_index().rename(columns={'level_1': 'Variant', 0: 'Infections'})
    exp_df = pd.concat(exp_dfs)  # .stack().reset_index().rename(columns={'level_1': 'Variant', 0:'Infections'})

    res = pd.DataFrame(ret)
    res['vx_day'] = pd.to_datetime(res['vx_day'])
    return df, exp_df, res



if __name__ == '__main__':
    n_reps = 2
    vx_res = 5
    scenarios = make_scenarios(n_reps=n_reps, vx_res=vx_res)

    print('Building scenarios...')
    sims = sc.parallelize(make_sim, iterkwargs=scenarios)

    print('Running simulations...')
    msim = cv.MultiSim(sims)
    msim.run()

    print('Processing and saving results...')
    df, exp_df, res = create_dfs(msim)

    sc.saveobj(f'{resfolder}/inf.obj', df)
    sc.saveobj(f'{resfolder}/exp.obj', exp_df)
    sc.saveobj(f'{resfolder}/res.obj', res)
    print('Done.')
