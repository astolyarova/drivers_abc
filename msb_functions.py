import numpy as np
import pandas as pd
import scipy
import statsmodels.formula.api as smf
import statsmodels.api as sm


def sim_adapt(N, sel, sel_scale, sites, rate, neut):
    N = int(round(pow(10, N)))
    sel_scale = pow(10, sel_scale)
    sel = pow(10, sel)
    mtypes = np.arange(0, len(sites), dtype=int)
    sites = np.array(sites)
    rate = np.array(rate)
    Nsites = sum(sites)
    Npass = Nsites - N
    Npass = max(0, Npass)
    mtype = np.repeat(mtypes, sites)
    adapt = []
    while not sum(adapt):
        adapt = np.random.dirichlet([sel]*N,size=1)[0]
        adapt = np.nan_to_num(adapt,nan=0,posinf=1e20,neginf=0)
    adapt = adapt / sum(adapt)

    adapt = adapt * sel_scale
    adapt = np.concatenate([adapt,np.zeros(Npass)])
    np.random.shuffle(adapt)

    if neut and Npass > 0:
        mtype0 = mtype[adapt == 0]
        rate0 = rate[mtype0]
        rate0 = np.sum(rate0, axis=1)
        mutation_counts0 = np.random.poisson(rate0)
        is_mutation0 = np.repeat(mtype0, mutation_counts0)
    else:
        is_mutation0 = np.array([])
        mutation_counts0 = np.array([])
        mtype0 = []

    mtype = mtype[adapt > 0]
    adapt1 = adapt[adapt > 0]
    rate = rate[mtype]
    rate = np.minimum(rate * np.array([adapt1 + 1]).T, 1)

    mutation_counts = np.random.binomial(n=1,p=rate)
    mutation_counts = np.sum(mutation_counts, axis=1)
    is_mutation = np.concatenate([np.repeat(mtype, mutation_counts),is_mutation0])
    unique, counts = np.unique(is_mutation, return_counts=True)
    is_mutation = dict(zip(unique, counts))
    mutations = [is_mutation.get(x,0) for x in mtypes]

    return(mutations, sum(mutation_counts>0)+sum(mutation_counts0>0), list(mutation_counts0)+list(mutation_counts), adapt, list(mtype0)+list(mtype))


def sim_stats(N, sel, sel_scale, samples, sites, intercept, slope, rate0, neut, k):
    results = []
    slope = np.log(np.array([samples]).T) * np.array(slope)
    rate1 = np.exp(intercept + slope).T
    for i in range(k):
        [mutations, mutation_counts,  n_eff, adapt, mtype] = sim_adapt(N, sel, sel_scale, sites, rate1, neut)
        if sum(mutations) == 0:
            results.append([-1,-1,-1,-1])
        else:
            mutations = np.array(mutations)
            n_eff = n_eff
            if not mutation_counts:
                results.append([-1,-1,-1,-1])
            else:
                ds=pd.DataFrame({'n':mutations,'rate0':rate0}) #,'predictor':predictor}) #'slope':slope1,'intercept':intercept1,'ntotal':ntotal1})
                m0_pois = smf.glm(formula = "n ~ 1", data=ds, offset=rate0,family=sm.families.Poisson()).fit()
                llf=m0_pois.llf
                results.append([np.log(sum(mutations)),np.log(mutation_counts),np.log(max(n_eff)),np.log(-llf)])#,m2_pois.llf])
    results = [np.array(x) for x in results]
    return(results)


def read_cvsq(file):
    rm_quote = lambda x: x.replace('"', '')
    df = pd.read_csv(file, engine='python', delim_whitespace=True).replace('"','', regex=True)
    df = df.rename(columns=rm_quote)
    df = df.fillna(0)
    return(df)


def read_data(file_spectra, file_sites, file_ntotal, file_stats, cancer, gene):
    df_spectra = read_cvsq(file_spectra)
    df_sites = read_cvsq(file_sites)
    df_samples = read_cvsq(file_ntotal)
    obs=[0,0,0,0]
    if file_stats:
        df_data = read_cvsq(file_stats)
        df_data = df_data.loc[(df_data['cancer'] == cancer) & (df_data['gene'] == gene)]
        ds = df_data.iloc[0]
        obs = [np.array([np.log(ds['n']),np.log(ds['n_uniq']),np.log(ds['n_max']),np.log(-ds['LL0_ABC'])])]

    ds_spectra = df_spectra.loc[(df_spectra['cancer'] == cancer)]
    ds_samples = df_samples.loc[(df_samples['cancer'] == cancer)]
    df_sites.loc[:,'mutation'] = df_sites.loc[:,'context'] + '_' + df_sites.loc[:,'alt']
    ds_sites = df_sites.loc[(df_sites['gene'] == gene)]
    ds_params = pd.merge(ds_spectra,ds_sites,on=['mutation'])
    sites = ds_params['n_sites']
    slope = list(ds_params['mu1'][sites > 0])
    intercept = list(ds_params['mu0'][sites > 0])
    samples = list(ds_samples['ntotal'])
    rate0 = list(ds_params['mu'][sites > 0])
    sites = list(sites[sites > 0])

    rate0 = np.log(np.exp(rate0)*np.array(sites)*len(samples))

    ksites = len(sites)
    predictor = [str(x) for x in range(ksites)]

    return(intercept, slope, rate0, sites, ds_params['mutation'], samples, obs)
