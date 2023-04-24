import argparse
import numpy as np
from msb_functions import *

parser = argparse.ArgumentParser(description='')
parser.add_argument('--abc', dest='abc', type=str)
parser.add_argument('--spectra', dest='spectra', type=str)
parser.add_argument('--sites', dest='sites', type=str)
parser.add_argument('--ntotal', dest='ntotal', type=str)
parser.add_argument('--cancer', dest='cancer', type=str)
parser.add_argument('--gene', dest='gene', type=str)
parser.add_argument('--samples', dest='samples', type=int)

args = parser.parse_args()

# load data

cancer = args.cancer
gene = args.gene
[intercept, slope, rate0, sites, mut, samples, obs] = read_data(args.spectra,args.sites,args.ntotal,None,cancer,gene)

abc = pd.read_csv(args.abc, engine='python', delim_whitespace=True)
abc = abc[(abc['cancer'] == cancer) & (abc['gene'] == gene)]
maxstep = max(abc['step'])
abc = abc[abc['step'] == maxstep]
n_samples = abc.shape[0]
samples = list(samples)

result = []
sim_samples = args.samples
if not sim_samples:
    sim_samples = len(samples)

result = []

# simulate

print('cancer gene log_N log_alpha log_s data_samples sim_samples n_driver n_all n_driver_uniq n_all_uniq n_driver_max n_all_max sites LL0 N99')
for k in range(n_samples):
    samples1 = list(np.random.choice(samples,sim_samples,replace=True))
    if sim_samples == len(samples):
        samples1 = samples
    slope1 = np.log(np.array([samples1]).T) * np.array(slope)
    rate1 = np.exp(intercept + slope1).T

    # assign driver mutations

    [mutations, mutation_counts, n_eff, adapt, mtype] = sim_adapt(abc['log_N'].iloc[k], abc['log_alpha'].iloc[k], abc['log_s'].iloc[k], sites, rate1, 0)
    adapt1 = sorted(adapt / sum(adapt) )[::-1]
    adapt1 = np.cumsum(adapt1)
    N99 = sum(adapt1 <= 0.99) + 1

    # simulate accumulation of driver mutations

    res = sim_stats(abc['log_N'].iloc[k], abc['log_alpha'].iloc[k], abc['log_s'].iloc[k], samples1, sites, intercept, slope, rate0, 0, 1)
    [n, n_uniq, nmax, ll] = [np.exp(x) for x in res[0]]
    n_uniq=n/n_uniq

    # simulate accumulation of driver and passenger mutations

    res = sim_stats(abc['log_N'].iloc[k], abc['log_alpha'].iloc[k], abc['log_s'].iloc[k], samples1, sites, intercept, slope, rate0, 1, 1)
    [n_neut, n_uniq_neut, nmax_neut, ll_neut] = [np.exp(x) for x in res[0]]
    n_uniq_neut=n_neut/n_uniq_neut

    print(cancer, gene, abc['log_N'].iloc[k], abc['log_alpha'].iloc[k], abc['log_s'].iloc[k], len(samples), sim_samples, n,n_neut,n_uniq, n_uniq_neut,nmax,nmax_neut,sum(sites),ll_neut,N99)

