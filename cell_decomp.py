
import pymc as pm
import numpy as np

#==============================================================
def model_stats(idata, prop_vec, n_clusts):
#==============================================================
    mean_post = np.mean(idata.posterior['beta'][0],axis=0)

    # #Calculate relative proportions
    Nd = np.sum(mean_post, axis=1) 
    Nd = np.asarray([Nd for i in range(n_clusts)]).T 
    from scipy.stats import linregress
    mean_post = np.divide(np.mean(idata.posterior['beta'][0],axis=0),Nd)
    line_fit=linregress(np.ravel(prop_vec), np.ravel(mean_post))
    return(mean_post, line_fit.rvalue**2)


##########################################################
##########################################################
#MY MODELS
def basic_pymc(n_clusts, n_spots, n_genes, ref_exp, spots, noise_type, likelihood):
    #Simple Linear regression
    with pm.Model(coords={"celltypes": np.arange(n_clusts),
                        "spots": np.arange(n_spots),
                        "genes": np.arange(n_genes) }) as basic_model:
        #Declare data 
        mean_exp = pm.Data('mean_exp', ref_exp, mutable=False, dims=['celltypes','genes'])
        # Priors for unknown model parameters
        beta=pm.HalfNormal("beta", sigma=1, dims=['spots','celltypes'])
        lmd= pm.Deterministic('lmd', pm.math.dot(beta, mean_exp), dims=['spots','genes'])
        #Convert from proportions to counts
        N_g = pm.Data('N_g', np.sum(spots, 1).reshape(n_spots,1), mutable=False)
        #Likelihood of observed data given Poisson rates
        if likelihood == 'poisson': y=pm.Poisson("y", mu=lmd*N_g, observed=spots)
        if likelihood == 'negbin': y=pm.NegativeBinomial("y", mu=lmd*N_g,alpha=1, observed=spots)
    return(basic_model)


def epsilon_pymc(n_clusts, n_spots, n_genes, ref_exp, spots, noise_type, likelihood):
    with pm.Model(coords={"celltypes": np.arange(n_clusts),
                        "spots": np.arange(n_spots),
                        "genes": np.arange(n_genes),
                        "1": np.arange(1) }) as epsilon_model:
        #Declare data 
        mean_exp = pm.Data('mean_exp', ref_exp, mutable=False, dims=['celltypes','genes'])
        # Priors for unknown model parameters
        beta=pm.HalfNormal("beta", sigma=1, dims=['spots','celltypes']) # celltype proportions
        if noise_type == 'gaussian': eps=pm.Normal("eps", mu=0, sigma=1, dims=['spots', 'genes'])
        if noise_type == 'exponential': eps=pm.Exponential("eps", lam=1, dims=['spots', 'genes']) # random noise at each spot and gene

        lmd= pm.Deterministic('lmd', np.exp(eps)*pm.math.dot(beta, mean_exp), dims=['spots','genes'])
        #Convert from proportions to counts
        N_g = pm.Data('N_g', np.sum(spots, 1).reshape(n_spots,1), mutable=False)

        #Likelihood of observed data given Poisson rates
        if likelihood == 'poisson': y=pm.Poisson("y", mu=lmd*N_g, observed=spots)
        if likelihood == 'negbin': y=pm.NegativeBinomial("y", mu=lmd*N_g,alpha=1, observed=spots)
    return(epsilon_model)


def gamma_pymc(n_clusts, n_spots, n_genes, ref_exp, spots, noise_type, likelihood):
    #Poisson noise
    with pm.Model(coords={"celltypes": np.arange(n_clusts),
                        "spots": np.arange(n_spots),
                        "genes": np.arange(n_genes),
                        "1": np.arange(1) }) as gamma_model:
        #Declare data 
        mean_exp = pm.Data('mean_exp', ref_exp, mutable=False, dims=['celltypes','genes'])
        # Priors for unknown model parameters
        beta=pm.HalfNormal("beta", sigma=1, dims=['spots','celltypes']) # celltype proportions
        if noise_type == 'gaussian': gamma = pm.Normal("gamma", mu=0, sigma=1, dims=['1','genes'])
        if noise_type == 'exponential': gamma = pm.Exponential("gamma", lam=1, dims=['1','genes']) # random noise at each spot and gene

        lmd= pm.Deterministic('lmd', np.exp(gamma)*pm.math.dot(beta, mean_exp), dims=['spots','genes'])
        #Convert from proportions to counts
        N_g = pm.Data('N_g', np.sum(spots, 1).reshape(n_spots,1), mutable=False)

        #Likelihood of observed data given Poisson rates
        if likelihood == 'poisson': y=pm.Poisson("y", mu=lmd*N_g, observed=spots)
        if likelihood == 'negbin': y=pm.NegativeBinomial("y", mu=lmd*N_g,alpha=1, observed=spots)
    return(gamma_model)


def alpha_pymc(n_clusts, n_spots, n_genes, ref_exp, spots, noise_type, likelihood):
    #Poisson noise
    with pm.Model(coords={"celltypes": np.arange(n_clusts),
                        "spots": np.arange(n_spots),
                        "genes": np.arange(n_genes),
                        "1": np.arange(1) }) as alpha_model:
        #Declare data 
        mean_exp = pm.Data('mean_exp', ref_exp, mutable=False, dims=['celltypes','genes'])
        # Priors for unknown model parameters
        beta=pm.HalfNormal("beta", sigma=1, dims=['spots','celltypes']) # celltype proportions

        if noise_type == 'gaussian': alpha = pm.Normal("alpha",mu=0, sigma=1, dims=['spots','1'])
        if noise_type == 'exponential': alpha = pm.Exponential("alpha", lam=1, dims=['spots','1']) # random noise at each spot and gene

        lmd= pm.Deterministic('lmd', np.exp(alpha)*pm.math.dot(beta, mean_exp), dims=['spots','genes'])
        #Convert from proportions to counts
        N_g = pm.Data('N_g', np.sum(spots, 1).reshape(n_spots,1), mutable=False)

        #Likelihood of observed data given Poisson rates
        if likelihood == 'poisson': y=pm.Poisson("y", mu=lmd*N_g, observed=spots)
        if likelihood == 'negbin': y=pm.NegativeBinomial("y", mu=lmd*N_g,alpha=1, observed=spots)
    return(alpha_model)


def noise_pymc(n_clusts, n_spots, n_genes, ref_exp, spots, noise_type, likelihood):
    #Poisson noise
    with pm.Model(coords={"celltypes": np.arange(n_clusts),
                        "spots": np.arange(n_spots),
                        "genes": np.arange(n_genes),
                        "1": np.arange(1) }) as noise_model:
        #Declare data 
        mean_exp = pm.Data('mean_exp', ref_exp, mutable=False, dims=['celltypes','genes'])
        # Priors for unknown model parameters
        beta=pm.HalfNormal("beta", sigma=1, dims=['spots','celltypes']) # celltype proportions

        if noise_type == 'gaussian': 
            alpha = pm.Normal("alpha",mu=0, sigma=1, dims=['spots','1'])
            gamma = pm.Normal("gamma", mu=0, sigma=1, dims=['1','genes'])
            eps=pm.Normal("eps", mu=0, sigma=1, dims=['spots', 'genes'])

        if noise_type == 'exponential': 
            alpha = pm.Exponential("alpha", lam=1, dims=['spots','1']) # random noise at each spot and gene
            gamma = pm.Exponential("gamma", lam=1, dims=['1','genes'])
            eps=pm.Exponential("eps", lam=1, dims=['spots', 'genes'])

        lmd= pm.Deterministic('lmd', np.exp(alpha)*np.exp(gamma)*np.exp(eps)*pm.math.dot(beta, mean_exp), dims=['spots','genes'])
        #Convert from proportions to counts
        N_g = pm.Data('N_g', np.sum(spots, 1).reshape(n_spots,1), mutable=False)

        #Likelihood of observed data given Poisson rates
        if likelihood == 'poisson': y=pm.Poisson("y", mu=lmd*N_g, observed=spots)
        if likelihood == 'negbin': y=pm.NegativeBinomial("y", mu=lmd*N_g,alpha=1, observed=spots)
    return(noise_model)