import covasim.utils as cvu
import covasim.parameters as cvpar
import covasim.defaults as cvd
import covasim.immunity as cvi
import covasim.interventions as cvintv
import numpy as np
from numpy.lib.arraysetops import isin
import sciris as sc

class treatment(cvintv.Intervention):
    '''
        Apply a treatment to eligible people.
        By default, the treatment has 90% efficacy against
        severe disease (equivalent to hospitalizations).
        Args:
            start_day    (int): intervention start day (default: 0, i.e. the start of the simulation)
            end_day      (int): intervention end day (default: no end)
            doses        (int): number of doses per treatment course
            prob       (float): probability of successful treatment uptake if eligible
            rel_severe (float): relative change in severe probability (default: 89% efficacy)
            rel_crit   (float): relative change in critical probability (default: 95% efficacy
            rel_death  (float): relative change in death probability (default: 100% efficacy)
            rel_trans  (float): relative change in transmissibility (default: 10% efficacy)
            eligible    (list): who is eligible: list of states or comma-separated string; by default, ``'diagnosed, symptomatic'`` (NB: false positives are ignored)
            subtarget   (dict): subtarget intervention to people with particular indices (see test_num() for details)
            kwargs      (dict): passed to Intervention()
        **Examples**::
            tr = cv.treatment(start_day=30, end_day=60, rel_severe=0.8, rel_trans=0.1, eligible='all')
        '''

    def __init__(self, start_day=0, end_day=None, doses=10, prob=1.0, eligible='diagnosed, symptomatic',
                 rel_severe=0.11, rel_crit=0.05, rel_death=0, rel_trans=1, subtarget=None, **kwargs):
        super().__init__(**kwargs)  # Initialize the Intervention object
        self.start_day = start_day
        self.end_day = end_day
        self.doses_per_course = doses
        self.prob = prob
        self.rel_severe = rel_severe
        self.rel_crit = rel_crit
        self.rel_death = rel_death
        self.rel_trans = rel_trans
        if eligible in [None, 'default']:
            eligible = 'diagnosed, symptomatic'
        if eligible == 'all':
            eligible = 'exposed'
        if isinstance(eligible, str):
            eligible = [s.lstrip().rstrip() for s in eligible.split(',')]
        self.eligible = eligible
        self.subtarget = subtarget
        return

    def initialize(self, sim):
        ''' Fix the dates and store the doses '''
        super().initialize()
        self.start_day = cvintv.preprocess_day(self.start_day, sim)
        self.end_day = cvintv.preprocess_day(self.end_day, sim)
        self.doses = np.zeros(sim.n, dtype=cvd.default_int)  # Number of doses given per person of this treatment
        self.n_treated = 0
        self.deaths_averted = 0
        self.on_treat = np.zeros(sim.n, dtype=bool)  # Person is currently being treated
        self.treat_dates = [[]] * sim.n  # Store the dates when people are treated
        return

    def apply(self, sim):
        ''' Perform treatment '''

        t = sim.t
        start_day = cvintv.get_day(self.start_day, self, sim)
        end_day = cvintv.get_day(self.end_day, self, sim)
        if t < start_day:
            return

        if end_day is None or t <= end_day:

            # Construct the testing probabilities piece by piece -- complicated, since need to do it in the right order
            on_treat_inds = sc.findinds(self.on_treat)  # People currently on treatment, do not want to double count
            treat_probs = np.full(sim.n, self.prob, dtype=np.float)  # Begin by assigning equal treatment probability to everyone
            if self.subtarget is not None:
                subtarget_inds, subtarget_vals = cvintv.get_subtargets(self.subtarget, sim)
                treat_probs[subtarget_inds] = subtarget_vals  # People being explicitly subtargeted
            if self.eligible:
                for cond in self.eligible:
                    treat_probs *= sim.people[cond]

            treat_probs[on_treat_inds] = 0
            treat_inds = cvu.true(cvu.binomial_arr(treat_probs))  # Calculate who actually gets treated

            self.on_treat[treat_inds] = True
            sim.people['rel_trans'][treat_inds] *= self.rel_trans

            # Find those who become severe/die, apply efficacy and determine if they no longer become severe/die
            sev_inds = treat_inds[np.isfinite(sim.people.date_severe[treat_inds])]
            crit_inds = sev_inds[np.isfinite(sim.people.date_critical[sev_inds])]
            dead_inds = crit_inds[np.isfinite(sim.people.date_dead[crit_inds])]
            sev_probs = np.full(len(sev_inds), 1 - self.rel_severe)
            crit_probs = np.full(len(crit_inds), 1 - self.rel_crit)
            dead_probs = np.full(len(dead_inds), 1 - self.rel_death)
            no_longer_sev_inds = sev_inds[cvu.true(cvu.binomial_arr(sev_probs))]  # Calculate who doesn't get severe disease
            no_longer_dead_inds = np.intersect1d(no_longer_sev_inds, dead_inds)
            # no_longer_crit_inds = crit_inds[cvu.true(cvu.binomial_arr(crit_probs))]  # Calculate who doesn't become critical (overlap with above)
            # no_longer_dead_inds = dead_inds[cvu.true(cvu.binomial_arr(dead_probs))]  # Calculate who doesn't become critical (overlap with above)

            # New durations and dates of recovery for these folks
            durpars = sc.dcp(sim['dur'])

            # First those who don't become severe (they therefore never become critical or die)
            sim.people.date_severe[no_longer_sev_inds] = np.nan
            sim.people.date_critical[no_longer_sev_inds] = np.nan
            sim.people.date_dead[no_longer_sev_inds] = np.nan
            dur_mild2rec = cvu.sample(**durpars['mild2rec'], size=len(no_longer_sev_inds))
            sim.people.date_recovered[no_longer_sev_inds] = sim.people.date_symptomatic[
                                                                no_longer_sev_inds] + dur_mild2rec  # Date they recover
            sim.people.dur_disease[no_longer_sev_inds] = sim.people.dur_exp2inf[no_longer_sev_inds] + \
                                                         sim.people.dur_inf2sym[
                                                             no_longer_sev_inds] + dur_mild2rec  # Store how long this person had COVID-19

            # # Then those who don't become critical (they therefore never die)
            # sim.people.date_critical[no_longer_crit_inds] = np.nan
            # sim.people.date_dead[no_longer_crit_inds] = np.nan
            # dur_sev2rec = cvu.sample(**durpars['sev2rec'], size=len(no_longer_crit_inds))
            # sim.people.date_recovered[no_longer_crit_inds] = sim.people.date_severe[
            #                                                      no_longer_crit_inds] + dur_sev2rec  # Date they recover
            # sim.people.dur_disease[no_longer_crit_inds] = sim.people.dur_exp2inf[no_longer_crit_inds] + \
            #                                               sim.people.dur_inf2sym[no_longer_crit_inds] + \
            #                                               sim.people.dur_sym2sev[
            #                                                   no_longer_crit_inds] + dur_sev2rec  # Store how long this person had COVID-19
            #
            # # Then those who don't die (they don't die)
            # sim.people.date_dead[no_longer_dead_inds] = np.nan
            # dur_crit2rec = cvu.sample(**durpars['crit2rec'], size=len(no_longer_dead_inds))
            # sim.people.date_recovered[no_longer_dead_inds] = sim.people.date_critical[
            #                                                      no_longer_dead_inds] + dur_crit2rec  # Date they recover
            # sim.people.dur_disease[no_longer_dead_inds] = sim.people.dur_exp2inf[no_longer_dead_inds] + \
            #                                               sim.people.dur_inf2sym[no_longer_dead_inds] + \
            #                                               sim.people.dur_sym2sev[no_longer_dead_inds] + \
            #                                               sim.people.dur_sev2crit[
            #                                                   no_longer_dead_inds] + dur_crit2rec  # Store how long this person had COVID-19

            # For people who recovered today, take off treatment and reset rel trans
            finish_inds = on_treat_inds[sim.people.date_recovered[on_treat_inds] == t]
            self.on_treat[finish_inds] = False
            sim.people['rel_trans'][finish_inds] /= self.rel_trans

            # Update counters
            factor = sim.rescale_vec[t]
            self.doses[treat_inds] += int(self.doses_per_course*factor)  # twice daily
            self.n_treated += len(treat_inds)*factor
            self.deaths_averted += len(no_longer_dead_inds)*factor
            for t_ind in treat_inds:
                self.treat_dates[t_ind].append(sim.t)

        return
