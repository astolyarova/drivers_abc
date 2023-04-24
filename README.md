
This page containts simulation code from the manuscript *"Mutation bias in driver genes reveals distribution of effects of oncogenic mutations"* (AV Stolyarova, GA Bazykin 2023).  

We developed an ABC-based approach to reconstruct the distribution of fitness effects (DFE) of mutations in driver genes based on several characteristics of patterns of mutations observed in these genes, including the mutational spectrum bias (MBS). The strength of MBS describes how the distribution of mutations across 96 three-nucleotide contexts differs from the same distribution of neutral mutations in non-driver genes in the same cancer type. 

Here, we performed ABC for 30 driver genes in 22 cancer type datasets (total 74 cancer-gene pairs). DFE in a driver gene is described by parameters *N* (number of potential driver mutations in a gene), *alpha* (parameter of the DFE, 1/*alpha* proportional to variance of fitness effects) and *s* (scaling coefficient). ABC uses simulations to generate the posterior distribution of the parameters giving best fit to the data. File `abc_results.out` contains parameters of 5000 simulations accepted in the final step of simulated annealing ABC (SABC), preformed separately for each of 74 cancer-gene pairs. For more details, see manuscript.

For any cancer-gene pair from `gene_cancer.out`, provided scripts can be used to simulate accumulation of mutations in this gene in a given sample size, e.g. the expected number of observed mutations, number of unique mutations, the strength of mutational spectrum bias and the probability of a mutation in a specific context to be driver.

**Predict mutations in driver genes**

The following command will use DFE parameters estimated for *KRAS* gene in PAAD (pancreatic adenocarcinoma) cancer type dataset to simulate accumulation of mutations in this gene in a given sample size (here, sample size = 802 = the size of PAAD cohort). 

    python3 predict_mutations.py --abc abc_results.out \
    --spectra spectra.csv --ntotal samples.out \
    --sites sites_by_genes.out --cancer PAAD --gene KRAS --samples 802

It will use the same simulation script used for ABC to perform 5000 simulations under parameters accepted in the final iteration of ABC and produce summary statistics for each simulation. The calculated statistics include:

Summary statistics used in ABC:
* the number of unique observed mutations in the gene
* the total number of observed mutations
* log-likelihood of the observed spectrum of mutations under neutral null-model (a measure of MSB)
* recurrence of the most frequent mutation

Other:
* the number of unique observed driver mutations 
* the total number of observed driver mutations
* N99 - the number of potential driver mutations expected to carry 99% of DFE (N99 <= N)

Similarly, the script can be run for any driver gene with estimated DFE in the corresponding cancer type (`gene_cancer.out` has the list of all available gene-cancer pairs). To generate predictions under custom DFE, instead of `abc_results.out` use your file of the same format containing the values for *N*, *alpha* and *s*.

**Predict probability of the observed mutation to be driver**

For each three-nucleotide context, the following command will calculate the probability of the mutation in this context observed in the given sample size to be driver, under condition that it wasn't sampled before (i.e. in the smaller sample sizes):

    python3 predict_probabilities.py --abc abc_results.out \
    --spectra spectra.csv --ntotal samples.out \
    --sites sites_by_genes.out --cancer PAAD --gene KRAS --samples 802

This script again uses parameters of 5000 simulations accepted in the final iteration of ABC to generate DFEs for given driver gene. For each possible mutation in the gene, by knowing its neutral mutation rate $\mu$ and the probability of fixation $x$, we can calculate the probability to sample this mutation in a single cancer sample $p = \mu* x$. The  probability to sample it for the first time in the *k*-th sample is:

$$ p(k) = p * (1 - p ^{k-1}) $$

For each mutation type, we can calculate the expected fraction of driver mutations in newly sampled mutations: 

$$ prob_{driver} = \sum_{drivers} p(k) / \sum_{all} p(k) $$

The results for 5000 simulations are pooled together to calculate the output probabilities. This command was used to produce Figures 5cd of the manuscript. Again, in can be run for any driver gene-cancer pair from `gene_cancer.out` file and for any sample size.
 
**Input files**

* `samples.out` - number of cancer samples within cancer type datasets and number of misense mutations per sample
* `spectra.csv` - neutral mutation rates in 96 three-nucleotide contexts, estimated in universal non-driver genes within cancer type datasets;
* `sites_by_genes.out` - occurrences of three-nucleotide contexts in driver genes

* `stat_genes.csv` - statistics for individual driver genes with no less than 50 observed mutations
* `stat_roles.csv` - statistics for driver genes pooled by role (oncogenes and TSGs)

* `abc_results.out` - accepted simulations

The scripts require `numpy` and `scipy` packages.
