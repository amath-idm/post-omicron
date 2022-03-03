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
debug = 0
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

variants = {
    'None': {
        'rel_beta_next': 3.5,
        'rel_imm_omicron_next': 0.5,
        'rel_imm_WT_next': 0.02,
        'rel_severe_next': 1,
        'n_imports': 0
    },
    'Emerged from Omicron': {
        'rel_beta_next': 3.5,
        'rel_imm_omicron_next': 0.5,
        'rel_imm_WT_next': 0.02,
        'rel_severe_next': 1,
        'n_imports': 25
    },
    'Emerged from WT': {
        'rel_beta_next': 3.5,
        'rel_imm_omicron_next': 0.02,
        'rel_imm_WT_next': 0.5,
        'rel_severe_next': 1,
        'n_imports': 25
    },
    'New cluster': {
        'rel_beta_next': 3.5,
        'rel_imm_omicron_next': 0.02,
        'rel_imm_WT_next': 0.02,
        'rel_severe_next': 1,
        'n_imports': 25
    },
    'Emerged from Omicron, more severe': {
        'rel_beta_next': 3.5,
        'rel_imm_omicron_next': 0.5,
        'rel_imm_WT_next': 0.02,
        'rel_severe_next': 3.5,
        'n_imports': 25
    },
    'Emerged from WT, more severe': {
        'rel_beta_next': 3.5,
        'rel_imm_omicron_next': 0.02,
        'rel_imm_WT_next': 0.5,
        'rel_severe_next': 3.5,
        'n_imports': 25
    },
    'New cluster, more severe': {
        'rel_beta_next': 3.5,
        'rel_imm_omicron_next': 0.02,
        'rel_imm_WT_next': 0.02,
        'rel_severe_next': 3.5,
        'n_imports': 25
    }
}

new_variant_days = ['2022-02-25', '2022-04-25', '2022-08-25']
vaccine_prime = [0, 0.2, 0.4, 0.6, 0.8, 1]
vaccine_boost = [0, 0.2, 0.4, 0.6, 0.8, 1]

def make_vx_intv(vaccine='pfizer', boost=False, day=None, coverage=1):
    '''
    Return an age-targeted booster vaccine intervention
        * vaccine: the name of the vaccine, e.g. "pfizer"
        * booster: whether or not this is a booster intervention
        * day: day on which to vaccinate
        * coverage: vaccine coverage
    '''

    day = np.arange(cv.day('2022-02-15', start_date='2021-10-01'), cv.day('2022-12-15', start_date='2021-10-01')) if \
        day is None else day

    if boost:
        interval = 180
        subtarget = {'inds': lambda sim: cv.true((sim.people.doses==2) & ((sim.t - sim.people.date_vaccinated) >= interval)),
                         'vals': coverage/90}
        intv = cv.vaccinate_prob(vaccine=vaccine, days=day, prob=0, subtarget=subtarget, booster=True,
                                     do_plot=False)
    else:
        subtarget = {'inds': lambda sim: cv.true((sim.people.age >= 12)), 'vals': coverage/90}
        intv = cv.vaccinate_prob(vaccine=vaccine, days=day, prob=0, subtarget=subtarget, do_plot=False)
    return sc.promotetolist(intv)


def scen_params():
    p = sc.objdict()
    p['vaccine_boost'] = vaccine_boost
    p['vaccine_prime'] = vaccine_prime
    p['nab_decay'] = list(nab_decay_params.keys())
    p['next_variant'] = list(variants.keys())
    p['new_variant_day'] = new_variant_days
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
    ''' RSA sim using historical waves'''

    # Variant definitions
    beta = cv.variant('beta', n_imports=0, days=0)
    delta = cv.variant('delta', n_imports=0, days=0)

    # Set the parameters
    pars = sc.objdict(
        n_agents=54e3,
        pop_scale=1e3,
        location='south africa',
        rand_seed=p.idx,
        nab_decay=nab_decay_params[p.nab_decay],
        pop_infected = 0,
        rescale = False,
        pop_type = 'hybrid',
        start_day = '2021-10-01',
        end_day = '2022-12-15',
        verbose = 0.1
    )

    # Intervention definitions
    beta_interventions = [
        cv.change_beta(pars['start_day'], 0.1),  # shut down any lingering delta
        cv.change_beta('2021-11-01', 0.8), # reopen to let omicron start to take off
        cv.change_beta('2021-11-15', 0.6), # behavior change during omicron
        cv.change_beta('2022-02-01', 1),  # behavior change after omicron
    ]

    # Add interventions
    interventions = sc.dcp(beta_interventions)

    # Diagnostic testing interventions
    # Testing is limited, 0.625% chance per day of testing while symptomatic with 0.8% chance per day if critical
    symp_prob = 0.00625
    st = {'inds': lambda sim: cv.true(sim.people.critical), 'vals': 0.008}
    interventions += [cv.test_prob(symp_prob=symp_prob, subtarget=st, start_day=0, test_delay=1)]

    # Create historical waves
    wave_dist = {'dist': 'normal', 'par1': 0, 'par2': 5 * 7 / 2.355}
    interventions += [cv.historical_wave(['2020-07-15', '2021-01-12', '2021-07-08'], [0.3, 0.3, 0.35],
                                        variant=['wild', 'beta', 'delta'], dist=wave_dist)]

    # Create historical vaccination
    vax_days_prior_pf = sc.daydiff('2021-05-18', pars['start_day'])
    subtarget_pfizer = {'inds': lambda sim: cv.true((sim.people.age >= 18)), 'vals': 0.0025}
    interventions += [cv.historical_vaccinate_prob('pfizer', days=np.arange(-vax_days_prior_pf, 0), prob=0,
                                            subtarget=subtarget_pfizer)]
    # Add current vaccination
    interventions += [cv.vaccinate_prob('pfizer', days=np.arange(0, sc.daydiff(pars['start_day'], '2022-02-14')),
                                    prob=0, subtarget=subtarget_pfizer)]

    # Create omicron
    omicron = cv.variant('omicron', days=['2021-10-25'], n_imports=25 * pars['pop_scale'])

    # Add a new variant
    var_pars = variants[p.next_variant]
    variant_pars = sc.mergedicts(cvpar.get_variant_pars('omicron'), {'rel_beta': var_pars['rel_beta_next'],
                                                                     'rel_severe_prob': var_pars['rel_severe_next']})
    next_variant = cv.variant(variant_pars, label='next_variant', days=p.new_variant_day,
                              n_imports=var_pars['n_imports'] * pars['pop_scale'])

    future_vax_prime = make_vx_intv(coverage=p.vaccine_prime)
    future_vax_boost = make_vx_intv(coverage=p.vaccine_boost, boost=True)

    interventions += future_vax_prime
    interventions += future_vax_boost
    for intervention in interventions:
        intervention.do_plot = False

    sim = cv.Sim(pars=pars, interventions=interventions, variants=[beta, delta, omicron, next_variant])

    # Now initialize the sim and update immunity
    sim.initialize()
    # Get index of omicron in immunity matrx -- WARNING, FRAGILE!
    inv_variant_map = {v: k for k, v in sim['variant_map'].items()}
    for k, v in sim['vaccine_map'].items():
        inv_variant_map[v] = k + sim['n_variants']
    oi = inv_variant_map['omicron']  # Index of omicron results
    nvi = inv_variant_map['next_variant']  # Index of next variant results

    immunity = sim['immunity']

    immunity[oi, nvi] = var_pars['rel_imm_omicron_next']  # Relative immunity of omicron to next variant
    immunity[nvi, :] = sc.dcp(immunity[oi, :])
    immunity[nvi, :oi] = var_pars['rel_imm_WT_next']  # Relative immunity of next variant to variants prior to omicron
    immunity[nvi, oi] = var_pars['rel_imm_omicron_next']  # Relative immunity of next variant to omicron
    immunity[nvi, nvi] = 1
    immunity[nvi, nvi + 1:] = var_pars['rel_imm_WT_next']

    pars.update({'immunity': immunity})

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
    parser.add_argument('--root', default='post-omicron-cov-sweep', type=str)
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

    if debug:
        print('Note, debug run!')
        sim = build_run_sim(paramset=paramset, seed=seed)

    else:
        sims = sc.parallelize(build_run_sim, iterkwargs=iterkwargs)
        msim = cv.MultiSim(sims)
        msims = msim.split(chunks=nparams)

        for msim in msims:
            msim.reduce()

    if verbose >= 0:
        print('Done running')

    # Save summary measures for fitting
    d = sc.objdict()
    if debug:
        ref_sim = sim
        params = sim.run_info.keys()
    else:
        ref_sim = msim.base_sim
        params = msims[0].base_sim.run_info.keys()

    inv_variant_map = {v: k for k, v in ref_sim['variant_map'].items()}
    oi = inv_variant_map['omicron']  # Index of omicron results
    di = inv_variant_map['delta']  # Index of delta results
    for param in params:
        d[param] = []
    d.perc_peak = []
    d.doses = []
    d.n_vaccinated = []
    d.deaths = []
    d.infections = []
    d.severe = []
    d.infections_next_variant = []
    d.severe_next_variant = []
    if not debug:
        for msim in msims:
            for sim in msim.sims:
                for param in params:
                    d[param].append(sim.run_info[param])
                d.perc_peak.append(
                    np.max(sim.results['variant']['new_infections_by_variant'][oi, :]) / np.max(
                        sim.results['variant']['new_infections_by_variant'][di, :]))
                date_after = np.argmax(sim.results['date'] == sc.date('2022-02-14'))
                d.deaths.append(np.sum(sim.results['new_deaths'][date_after:]))
                d.doses.append(np.sum(sim.results['new_doses'][date_after:]))
                d.n_vaccinated.append(np.sum(sim.results['new_vaccinated'][date_after:]))
                d.infections.append(np.sum(sim.results['new_infections'][date_after:]))
                d.severe.append(np.sum(sim.results['new_severe'][date_after:]))
                d.infections_next_variant.append(
                    np.sum(sim.results['variant']['new_infections_by_variant'][oi + 1, date_after:]))
                d.severe_next_variant.append(
                    np.sum(sim.results['variant']['new_severe_by_variant'][oi + 1, date_after:]))
                peak_day = list(sim.results['variant']['new_infections_by_variant'][oi + 1, :]).index(
                    np.max(sim.results['variant']['new_infections_by_variant'][oi + 1, :]))

        sc.saveobj(f'{resfolder}/{args.root}_vx_rollout_data.obj', d)

    # Also create and save data for making additional plots
    plotdf = get_data_for_plots(msims, filename=f'{resfolder}/{args.root}_data_for_run_plot.obj')
    print('Done.')