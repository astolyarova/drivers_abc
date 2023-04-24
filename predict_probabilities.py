import sys
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
[intercept, slope, rate0, sites, muttypes, samples, obs] = read_data(args.spectra,args.sites,args.ntotal,None,cancer,gene)

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

prob_dr = {}
prob = {}
expected = {}

# simulate

for k in range(n_samples):

    if not k % 10:
        sys.stderr.write('simulation #' + str(k) + '/5000 done\n')
    # assign driver mutations

    [mutations, mutation_counts, n_eff, adapt, mtype] = sim_adapt(abc['log_N'].iloc[k], abc['log_alpha'].iloc[k], abc['log_s'].iloc[k], sites, np.array([np.exp(rate0)]).T, 1)

    # calculate conditional probabilities

    mutations = np.arange(len(mtype))
    for j in mutations:
        s=adapt[j]
        rate=np.exp(rate0[mtype[j]])/len(samples)/sites[mtype[j]]
        prop=min(1,(s+1)*rate)
        mut = muttypes.iloc[mtype[j]]
        expected[mut] = expected.get(mut, 0) + prop
        for nsamples in [sim_samples]:
            prop_not=pow((1-prop),(nsamples-1))
            weight=prop*prop_not
            prob[mut] = prob.get(mut, 0) + weight
            if s > 0:
                prob_dr[mut] = prob_dr.get(mut, 0) + weight

print('cancer gene samples mutation expected_mutations driver_probability')
for mut in prob:
    print(cancer, gene, nsamples, mut, expected[mut] / n_samples, prob_dr.get(mut, 0) / prob[mut])
