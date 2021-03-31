# -*- coding: utf-8 -*-
"""
Carriage duration parameter fitting for pneumo IBM
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotnine import *
from plotnine.data import mpg
from collections import Counter
import scipy.stats as stats
from scipy.stats import gamma
import emcee
import corner
import time
import math

carr_dat = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/carriage_episodes_AL.csv")
paramsdf = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/carrdurMCMCparams.csv")
seroparams = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/seroparams.csv")

# remove any rows with non positive length of carr
carr_dat = carr_dat.loc[np.where(carr_dat.length > 0)].reset_index(drop = True) 
carr_dat.t_f = carr_dat.t_i + carr_dat.length
# remove dates columns`
carr_dat = carr_dat.iloc[:,2:len(carr_dat.columns)+1]
# create ID instead of bday
carr_dat['ID'] = (carr_dat.sort_index().groupby(['bday'], sort=False).ngroup()+1)
carr_dat = carr_dat.drop(['bday','lanes'], axis = 1)
carr_dat = carr_dat.sort_values('t_i')
carr_dat['coln_num'] = carr_dat.groupby('ID').cumcount()
carr_dat.iloc[carr_dat.serotype == '15B', 4] = '15B/C'
carr_dat.iloc[carr_dat.serotype == '15C', 4] = '15B/C'
t = np.arange(0,max(carr_dat.t_f)+1)
max_t = max(t)
unique_sero = np.unique(carr_dat.serotype)

def logprior(theta, ndim):

    if ndim == 2:
        A, B = theta
        
    if ndim == 3:
        A, B, alpha = theta

    # uniform prior on alpha, A and B
    alphamin = 0 # lower range of prior
    alphamax = 5.  # upper range of prior
    Amin = 100. # lower range of prior
    Amax = 300.  # upper range of prior
    Bmin = -2. # lower range of prior
    Bmax = 2.  # upper range of prior

    # set prior to 1 (log prior to 0) if in the range and zero (-inf) outside the range 
    if ndim == 3 and Amin < A < Amax and Bmin < B < Bmax and alphamin < alpha < alphamax:
        return 0.0
        
    if ndim == 2 and Amin < A < Amax and Bmin < B < Bmax:
        return 0.0
    
    return -np.inf

### log likelihood function
def loglikelihood(theta, prevcol, duration, ndim):

    if ndim == 2:
        alpha = 1
        A, B = theta
        
    if ndim == 3:
        A, B, alpha = theta

    mu = A*np.exp(-prevcol*B) # mean of gamma distribution
    beta = alpha/mu
    LL = np.sum(gamma.logpdf(x = duration, a = alpha, loc = 0, scale = 1/beta))
    return LL

def logposterior(theta, prevcol, duration, ndim):
    if not np.isfinite(logprior(theta, ndim)):
        return -np.inf
    return logprior(theta, ndim) + loglikelihood(theta, prevcol, duration, ndim)

def carrdur_mcmc(serotype, alpha_tf): 
    # function that runs MCMC to find best fitting parameters for a serotype
    # alpha_tf is 0/False if alpha will be fixed (alpha = 1) or 1/True if alpha will be estimated
    
    ds = carr_dat[carr_dat['serotype'] == serotype]
    
    if len(ds['coln_num']) == 1:
        prevcol_obs = np.array(ds['coln_num'])
        duration_obs = np.array(ds['length'])
    else:
        prevcol_obs = np.array([ds['coln_num']])
        duration_obs = np.array([ds['length']])

    ndim = 2+alpha_tf  # number of parameters in the model = 2 if alpha fixed and = 3 if alpha is estimated
    nwalkers = 50  # number of MCMC walkers
    nburn = 1000  # "burn-in" period to let chains stabilize
    nsteps = 5000  # number of MCMC steps to take
    
    alphaini = np.random.uniform(0, 5, nwalkers)
    Aini = np.random.uniform(50, 300, nwalkers)
    Bini = np.random.uniform(-2, 2, nwalkers)
    
    if ndim == 2: 
        inisamples = np.array([Aini, Bini]).T # initial samples
    else:
        inisamples = np.array([Aini, Bini, alphaini]).T # initial samples
    
    np.random.seed(0)
    
    sampler = emcee.EnsembleSampler(nwalkers, ndim, logposterior, args=[prevcol_obs, duration_obs, ndim])
    #%time 
    sampler.run_mcmc(inisamples, nsteps)
    print("done")
    
    samples_emcee = sampler.chain[:, nburn:, :].reshape((-1, ndim))
    
    resdict = {}
    resdict['serotype'] = serotype
    resdict['A_emcee'] = np.mean(samples_emcee[:,0])
    resdict['A_stdev'] = np.std(samples_emcee[:,0])
    resdict['B_emcee'] = np.mean(samples_emcee[:,1])
    resdict['B_stdev'] = np.std(samples_emcee[:,1])
    
    if ndim == 3:
        resdict['alpha_emcee'] = np.mean(samples_emcee[:,2])
        resdict['alpha_stdev'] = np.std(samples_emcee[:,2])
        fig = corner.corner(samples_emcee, labels=[r"$A$", r"$B$", r"$alpha$"], 
                            truths=[resdict['A_emcee'], resdict['B_emcee'], resdict['alpha_emcee']])
    else:
        fig = corner.corner(samples_emcee, labels=[r"$A$", r"$B$"], 
                            truths=[resdict['A_emcee'], resdict['B_emcee']])
    
    return resdict

def qualfit(serotype, prevcol): #resdict
    ds = carr_dat[carr_dat['serotype'] == serotype] #dataset
    prevcol_obs = np.array([ds['coln_num']])
    duration_obs = np.array([ds['length']])
    resdict = paramsdf[paramsdf['serotype'] == serotype] # comment out if providing resdict as param
    A = resdict['A_emcee']
    B = resdict['B_emcee']
    alpha = resdict['alpha_emcee']
    string = 'alpha = ' + str(round(alpha,2))
    string2 = 'alphaest'
    # if len(resdict) == 4:
    #     alpha = 1
    #     string = 'alpha = 1'
    #     string2 = 'alphafixed'
    # else:
    #     alpha = resdict['alpha_emcee']
    #     string = 'alpha = ' + str(round(alpha,2))
    #     string2 = 'alphaest'    
    # model
    x = np.sort(duration_obs[prevcol_obs == prevcol])
    y1 = stats.gamma.pdf(x, a=alpha, scale=(A*np.exp(-0*B))/alpha)
    #plt.plot(x, y1, "y-")
    # observed
    #plt.hist(duration_obs[prevcol_obs == prevcol], bins = 15)
    title = serotype + ', ' + string
    filenam = serotype + '_' + string2
    fig, axs = plt.subplots(2)
    fig.suptitle(title)
    axs[0].plot(x, y1, "y-")
    axs[1].hist(duration_obs[prevcol_obs == prevcol], bins = 15)
    fig.savefig(filenam+'.png', bbox_inches='tight')
    return 

# test
# sero19F_mcres = carrdur_mcmc('19F', True)
# sero19F_mcres_alphafix = carrdur_mcmc('19F', False)
# qualfit(serotype = '19F', resdict = sero19F_mcres, prevcol = 4)
# qualfit(serotype = '19F', resdict = sero19F_mcres_alphafix, prevcol = 4)

# params for all serotypes
paramslist = []
for i in unique_sero:
    resdict = carrdur_mcmc(i, True)
    paramslist.append(resdict)

paramsdf = pd.DataFrame(paramslist, 
                        columns= ['serotype', 'A_emcee','A_stdev', 'B_emcee', 'B_stdev', 
                                  'alpha_emcee', 'alpha_stdev'])
newfiledest = "//qdrive/homes/al817/Technical/Python/carrdurMCMCparams.csv"
paramsdf.to_csv(path_or_buf= newfiledest, index = False)

# initialise fake baby population
indivs = len(np.unique(carr_dat['ID']))
num_sero = len(unique_sero)
max_carr_cap = 1
children_data = np.empty([indivs, num_sero])
for i in range(num_sero):
    children_data[:,i] = np.array([np.random.choice([int(0),int(1)], size = indivs, replace = True, p = [0.99, 0.01])])

carriers_children = pd.DataFrame(data = children_data, columns = unique_sero)
carriers_children = carriers_children.astype(int) # change floating 0s and 1s to integers
prev_col_children = carriers_children.copy(deep = True) # serotype-specific previous colonizations
susc_children = max_carr_cap - carriers_children.sum(axis = 1)
carriers_children['susc'] = susc_children
carriers_children['agegroup'] = int(0)
carriers_children['age'] = carr_dat['age_d']

# initial FOIs
sero_events = dict(Counter(carr_dat['serotype']))
sero_events_df = pd.DataFrame(sero_events.items(), columns=['serotype', 'count'])
sero_events_df['initFOI'] = sero_events_df['count']/(len(carr_dat.ID.unique())*max_t)
sero_events_df = sero_events_df.drop(['count'], axis = 1)
# initial beta eff contact rates = initial FOI/average carriers existing at any given time
# avg carr existing at any given time = sum carr lengths/total time
sumcarrlength = carr_dat.groupby(['serotype'])['length'].sum().reset_index()
avgcarrpertime = sumcarrlength['length']/max_t
avgcarrpertime = pd.concat([sumcarrlength['serotype'].reset_index(drop=True), avgcarrpertime], axis=1)
FOImergedf = sero_events_df.merge(avgcarrpertime, left_on='serotype', right_on='serotype')
FOImergedf['beta'] = FOImergedf['initFOI']/FOImergedf['length']
beta_rates = FOImergedf.drop(['initFOI', 'length'], axis = 1)

carrdurparams = paramsdf.drop(['A_stdev', 'B_stdev', 'alpha_stdev'], axis = 1)
seroparams = carrdurparams.merge(beta_rates, left_on = 'serotype', right_on = 'serotype')
seroparams['maxcarrdur']=np.random.gamma(shape = seroparams['alpha_emcee'], scale = (seroparams['A_emcee']*np.exp(-prevcol0*seroparams['B_emcee']))/seroparams['alpha_emcee'])
seroparams.to_csv(path_or_buf= r"//qdrive/homes/al817/Technical/Python/seroparams.csv", index = False)
     
# beta (transmission or FOI) rates ## these will be different for children and adults but are the same for all serotypes
# values taken from Melegaro 2004 longitudinal carriage study in the UK
beta_rates = [0.012/indivs]*num_sero # per day
seroparams['beta'] = beta_rates
   
sigma = 0.4 # specific immunity
epsilon = 0.1 # non-specific immunity
theta = 0.01 # competition

def new_lambdas(carrier_person, ind_lambdas, prev_col, sigma): # function that computes new lambdas based on serotype switching event
    carrier_person = carrier_person.to_numpy()
    susc = carrier_person[num_sero]
    ind_lambdas = ind_lambdas.to_numpy()
    prev_col = prev_col.to_numpy()
    index = np.where(carrier_person[0:num_sero] == 0) # index of serotypes not carried
    non_index = np.bool_(carrier_person[0:num_sero]) # index of serotypes carried
    # draw from distribution of parameters with mean @ max likelihood
    eps = epsilon #abs(np.random.normal(epsilon))
    thet = theta #abs(np.random.normal(theta))
    sig = sigma #abs(np.random.normal(sigma))
    # if you're carrying max number of serotypes, then you don't feel any FOI for new serotype
    if susc <= 0:
        ind_lambdas[index] = 0
    # otherwise, your FOI changes based on your colonization history and competition
    else:
        ind_lambdas[index] = (1-(sig*(prev_col[index] > 0)))*ind_lambdas[index]*(1-thet*(susc == 1))
    # serotypes you are carrying change FOI based on colonization history and immunity
    ind_lambdas[non_index] = (1/250) + ((ind_lambdas[non_index]-(1/250))*math.exp(eps*sum(prev_col)))
    return(ind_lambdas)

carriers = carriers_children.copy(deep = True)
t = 0
stop_t = max_t #10000 # = 27 years  # t is in units of days
tot_lambda1 = sum(sero_events_df['initFOI'])
#(sum(carriers.iloc[:,0])*sero_events_df[sero_events_df['serotype']=='1'].initFOI) + 
#(indivs - sum(carriers.iloc[:,0]))*sero_rates_child[1,0]
    
# outputs for graph of time vs infected ppl
chunkSize = 100000
currentLength = chunkSize
plotsy = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = unique_sero)
plotsy.insert(loc = 0, column = 'time', value = np.nan)
    
next_ev = 1
start = time.time()
while t < stop_t and tot_lambda1 != 0:

    index_children = np.where(carriers.age <= 6205)[0]
    # add NAs to end of output dataframe if there is no space left to fill
    if next_ev > currentLength:
        plotsy = plotsy.append(pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = plotsy.columns), ignore_index = True)
        currentLength = currentLength + chunkSize
    
    # only record data when t > 9000 (i.e. last 3 years of time)
    #if t > 100: # record after 100th day
    plotsy.iloc[next_ev-1, 1:] = np.array(carriers.iloc[:,0:num_sero].sum(axis = 0))
    plotsy.time[next_ev-1] = t

    # calc the total FOI = contact rate * # of infected individuals for that serotype
    totFOIbysero = [0]*num_sero
    for i in range(len(seroparams)):
        totFOIbysero[i] = seroparams['beta'][i]*(carriers.iloc[index_children,0:num_sero].sum(axis = 0)[i])
    #midtime3_5 = time.time()
        
    # individual rates        
    ind_lambdas = carriers.iloc[:, 0:num_sero].copy(deep = True)
    
    for s in range(indivs):
        for j in range(num_sero):
            prevcol = prev_col_children.iloc[s,j]
            alpha = seroparams['alpha_emcee'][j]
            A = seroparams['A_emcee'][j]
            B = seroparams['B_emcee'][j]
            if ind_lambdas.iloc[s,j] == 0:
                ind_lambdas.iloc[s,j] = totFOIbysero[j]
            else:
                ind_lambdas.iloc[s,j] = 1/np.random.gamma(shape = alpha, scale = (A*np.exp(-prevcol*B))/alpha)
       
    
    #midtime4 = time.time()
        
    # tot.lambda = total rate of all/any events
    
    
    for i in range(10): # update 20 children at a time
        tot_lambda1 = ind_lambdas.to_numpy().sum()
        U1 = np.random.uniform(0,1)
        U1_tot_lambda1 = U1*tot_lambda1
        
        # choose person index ( = which person will this happen to?)
        cumsum_person = ind_lambdas.sum(axis = 1)
        which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        
        #midtime5 = time.time()

        # choose event that will take place (= which serotype will be toggled?)
        person_lambdas = new_lambdas(carriers.loc[which_person], ind_lambdas.loc[which_person], prev_col_children.loc[which_person], sigma)
        tot_lambda2 = sum(person_lambdas)
        U2 = np.random.uniform(0,1)
        U2_tot_lambda2 = U2*tot_lambda2

        which_indices = np.where(U2_tot_lambda2 < np.cumsum(person_lambdas))
        
        #midtime6 = time.time()
        
        # switching event (S --> I or I --> S)
        if which_indices[0].size > 0:
            which_sero = which_indices[0][0]
            carriers.iloc[which_person, which_sero] = ~carriers.iloc[which_person, which_sero] + 2

        # reupdate the susc score
        new_susc = max_carr_cap - sum(carriers.iloc[which_person, 0:num_sero])
    
        # update colonisation history if person acquired new serotype    
        if new_susc < carriers.susc.loc[which_person]:
            prev_col_children.iloc[which_person, which_sero] = prev_col_children.iloc[which_person, which_sero] + 1
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        
        #midtime7 = time.time()

        # choose next time step
    #U3 = np.random.uniform(0,1)
    #time_next = -math.log(U3)/tot_lambda1

    # update time and age
    t += 1#time_next
    next_ev += 1
    carriers.age += 1#time_next
        
    #midtime8 = time.time()
        
end = time.time()
total_time = end - start_time
plotsy = plotsy.dropna()
plotsy_melt= pd.melt(plotsy, id_vars = 'time')
pivotplots = plotsy_melt.pivot_table('value', 'time', 'variable', aggfunc='mean')
plot_fig = pivotplots.plot(kind='line')
#plot_fig.savefig("simulation2.pdf", bbox_inches = "tight")

#### checking that average derivatives at t = 0 is 0 i.e. that it's at equilibrium before the simulation runs since this ########
############# is what we assume from the parameters ########################################-----------------------------------------
arraysim = []
for x in range(1000):
    children_data = np.empty([indivs, num_sero])
    for i in range(num_sero):
        children_data[:,i] = np.array([np.random.choice([int(0),int(1)], size = indivs, replace = True, p = [0.999, 0.001])])

    carriers_children = pd.DataFrame(data = children_data, columns = unique_sero)
    carriers_children = carriers_children.astype(int) # change floating 0s and 1s to integers

    prevcolnull = carriers_children.iloc[:,0:num_sero].copy(deep = True)
    ind_lambdas1 = carriers_children.iloc[:, 0:num_sero].copy(deep = True)
    totFOIbysero1 = [0]*num_sero
    for i in range(num_sero):
         totFOIbysero1[i] = seroparams['beta'][i]*(carriers_children.iloc[:,0:num_sero].sum(axis = 0)[i])
    for s in range(indivs):
        for j in range(num_sero):
            prevcol = prevcolnull.iloc[s,j]
            alpha = seroparams['alpha_emcee'][j]
            A = seroparams['A_emcee'][j]
            B = seroparams['B_emcee'][j]
            if ind_lambdas1.iloc[s,j] == 0:
                ind_lambdas1.iloc[s,j] = totFOIbysero1[j]
            else:
                ind_lambdas1.iloc[s,j] = -(1/np.random.gamma(shape = alpha, scale = (A*np.exp(-prevcol*B))/alpha))
    arraysim.append(ind_lambdas1.iloc[:,0:num_sero].sum(axis = 0))
#max(ind_lambdas1.iloc[:,0:num_sero].sum(axis = 0)) # these should be around 0
    
arraysimdf = pd.DataFrame(arraysim, columns = arraysim[0].index)
meanrates = arraysimdf.mean() # 0.0007 average

# distribution of total previous colonisations in simulation vs in data
totprevcol = prev_col_children.sum(axis = 1)
plt.hist(totprevcol, bins = 15)

plt.hist(carr_dat['ID'].value_counts(), bins = 15)