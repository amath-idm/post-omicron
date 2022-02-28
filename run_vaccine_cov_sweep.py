'''

Script to show trade-off between second and third dose coverage
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
from utils import get_data_for_plots

module_path = os.path.abspath(os.path.join('..'))
if module_path not in sys.path:
    sys.path.append(module_path)

cv.check_version('>=3.1.1', die=True)
verbose = 0
plot_uncertainty = True
resfolder = 'results'

############################################################################
# Define variables to sweep over
############################################################################
fast_decay = dict(form='nab_growth_decay', growth_time=21, decay_rate1=np.log(2) / 100, decay_time1=100,
                  decay_rate2=np.log(2) / 250, decay_time2=200)
slow_decay = dict(form='nab_growth_decay', growth_time=21, decay_rate1=0.006741981, decay_time1=47,
                  decay_rate2=0.001764455, decay_time2=106)

nab_decay_params = {
    # 'vax_fast_nat_slow': dict(natural=slow_decay, vaccine=fast_decay),
    'both_fast': dict(natural=fast_decay, vaccine=fast_decay),
    # 'nat_fast_vax_slow': dict(natural=fast_decay, vaccine=slow_decay),
    # 'both_slow': dict(natural=slow_decay, vaccine=slow_decay),
}

vaccine_prime = [0, 0.2, 0.4, 0.6, 0.8, 1]
vaccine_boost = [0, 0.2, 0.4, 0.6, 0.8, 1]

def make_vx_intv(vaccine='pfizer', boost=False, coverage=1):
    '''
    Return a 90-day vaccine intervention
        * vaccine: the name of the vaccine, e.g. "pfizer"
        * booster: whether or not this is a booster intervention
        * coverage: vaccine coverage
    '''
        
    if boost:
        day = np.arange(cv.day('2021-11-01', start_date='2020-03-01'), cv.day('2022-02-01', start_date='2020-03-01'))
        interval = 180
        subtarget = {'inds': lambda sim: cv.true((sim.people.doses==2) & ((sim.t - sim.people.date_vaccinated) >= interval)),
                         'vals': coverage/90}
        intv = cv.vaccinate_prob(vaccine=vaccine, days=day, prob=0, subtarget=subtarget, booster=True,
                                     do_plot=False)
    else:
        day = np.arange(cv.day('2021-05-01', start_date='2020-03-01'), cv.day('2021-08-01', start_date='2020-03-01'))
        subtarget = {'inds': lambda sim: cv.true((sim.people.age >= 18)), 'vals': coverage/90}
        intv = cv.vaccinate_prob(vaccine=vaccine, days=day, prob=0, subtarget=subtarget, do_plot=False)
    return sc.promotetolist(intv)


def scen_params():
    p = sc.objdict()
    p['vaccine_boost'] = vaccine_boost
    p['vaccine_prime'] = vaccine_prime
    p['nab_decay'] = list(nab_decay_params.keys())
    dims = [len(v) for k, v in p.items()]
    npars = np.prod(dims)
    return p, npars


def get_run_params(index):
    params, nparams = scen_params()
    dims = [len(v) for k, v in params.items()]
    indices = index % nparams
    inds = np.unravel_index(indices, dims)
    f = sc.objdict()
    for ik, (k, v) in enumerate(params.items()):
        f.update({k: v[inds[ik]]})
    return f

############################################################################
# Make the sims
############################################################################


def make_sim(p):
    
    # Set the parameters
    pars = sc.objdict(
        n_agents=54e3,
        pop_scale=1e3,
        location='south africa',
        rand_seed=p.idx,
        nab_decay=nab_decay_params[p.nab_decay],
        beta=0.015,  # Base transmission probability per contact, per day
        pop_infected=70,  # Number of seed infections
        start_day='2020-03-01',  # First day of simulation
        end_day='2022-02-15',  # Last day of simulation
        interventions=[],  # Interventions to be added later
        analyzers=[],  # Analyzers to be added later
        use_waning=True,  # Enable waning immunity
        verbose=0,  # Turn off outputs to avoid notebook clutter
    )
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
        cv.change_beta('2021-11-25', 0.6)   # shut down
    ]

    pars['interventions'] += beta_interventions

    # Diagnostic testing interventions
    # Testing is limited, 0.625% chance per day of testing while symptomatic with 0.8% chance per day if critical
    symp_prob = 0.00625
    st = {'inds': lambda sim: cv.true(sim.people.critical), 'vals': 0.008}
    pars['interventions'] += [cv.test_prob(symp_prob=symp_prob, subtarget=st, start_day=0, test_delay=1)]

    future_vax_prime = make_vx_intv(coverage=p.vaccine_prime)
    future_vax_boost = make_vx_intv(coverage=p.vaccine_boost, boost=True)

    pars['interventions'] += future_vax_prime
    pars['interventions'] += future_vax_boost
    for intervention in pars['interventions']:
        intervention.do_plot = False

    sim = cv.Sim(pars=pars, variants=variants)

    sim.run_info = sc.objdict()
    sim.run_info.update(p)
    return sim


############################################################################
# Utils for building and running
############################################################################

def build_sim(index, idx):
    print(f'Building sim {index}, iteration {idx}')
    p = get_run_params(index) # get sim params
    print(p)
    p.idx = idx
    sim = make_sim(p) # get sim
    sim.run_info.update(p) # save info
    return sim


def run_sim(sim):
    # run sim
    sim.run()

    # shrink sim
    sim.shrink()

    return sim


def build_run_sim(paramset, seed):
    ''' Build then run the simulations '''
    sim = build_sim(paramset, seed)
    sim = run_sim(sim)
    return sim


############################################################################
# Run analyses
############################################################################

if __name__ == '__main__':

    # Define defaults for arguments
    from argparse import ArgumentParser

    parser = ArgumentParser()

    parser.add_argument('--nsamples', default=50, type=int)  # How many samples to run (stochastic uncertainty)
    parser.add_argument('--ncpus', default=None, type=int)  # CPUs
    parser.add_argument('--root', default='vax-cov-sweep', type=str)
    args = parser.parse_args()

    # Customize, initialize timer
    T = sc.tic()

    # Get parameters and number of draws
    # Run!
    params, nparams = scen_params()
    nsim_total = nparams * args.nsamples

    # Make sims
    sc.heading(f'Making {nsim_total} sims...')
    iterkwargs = sc.autolist()
    for paramset in range(nparams):
        for seed in range(args.nsamples):
            iterkwargs += dict(paramset=paramset, seed=seed)

    sims = sc.parallelize(build_run_sim, iterkwargs=iterkwargs)
    msim = cv.MultiSim(sims)
    msims = msim.split(chunks=nparams)

    for msim in msims:
        msim.reduce()

    if verbose >= 0:
        print('Done running')

    # Save summary measures for fitting
    d = sc.objdict()
    ref_sim = msim.base_sim
    params = msims[0].base_sim.run_info.keys()

    for param in params:
        d[param] = []
    d.doses = []
    d.n_vaccinated = []
    d.deaths = []
    d.infections = []
    d.severe = []
    for msim in msims:
        for sim in msim.sims:
            for param in params:
                d[param].append(sim.run_info[param])
            date_after = np.argmax(sim.results['date'] == sc.date('2021-05-01'))
            d.deaths.append(np.sum(sim.results['new_deaths'][date_after:]))
            d.doses.append(np.sum(sim.results['new_doses'][date_after:]))
            d.n_vaccinated.append(np.sum(sim.results['new_vaccinated'][date_after:]))
            d.infections.append(np.sum(sim.results['new_infections'][date_after:]))
            d.severe.append(np.sum(sim.results['new_severe'][date_after:]))
                
    sc.saveobj(f'{resfolder}/{args.root}_vx_rollout_data.obj', d)
        

    # Also create and save data for making additional plots
    plotdf = get_data_for_plots(msims, filename=f'{resfolder}/{args.root}_data_for_run_plot.obj', historical_wave=False)
    print('Done.')