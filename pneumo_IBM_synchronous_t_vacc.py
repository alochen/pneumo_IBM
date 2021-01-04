########## PNEUMO IBM with fitted sero params imported and with synchronous time where 20 events occur at each ############
########## time step, and a migration rate such that when a serotype goes extinct, it is reseeded into the     ############
########## population in the form of one infectious carrier                                                    ############

#import os
import pandas as pd
import numpy as np
import math
import random
import time
import glob
#import matplotlib
#import matplotlib.pyplot as plt

#os.getcwd()
#os.chdir('Q:\\Technical\\Python')
# function that computes new lambdas based on serotype switching event
def new_lambdas(carrier_person, ind_lambdas, prev_col, sigma):    
    susc = carrier_person.susc
    vacc = carrier_person.vacc
    
    index = np.where(carrier_person[0:num_sero] == 0) # index of serotypes not carried
    non_index = np.bool_(carrier_person[0:num_sero]) # index of serotypes carried
    vt_ind = index[0][np.isin(carrier_person.index[index], vt)] # index of VT serotypes not carried 
    # draw from distribution of parameters with mean @ max likelihood
    eps = epsilon#np.random.normal(epsilon)
    thet = theta#np.random.normal(theta)
    sig = sigma#np.random.normal(sigma)
    # if you're carrying max number of serotypes, then you don't feel any FOI for new serotype (regardless of vaccination)
    if susc <= 0:
        ind_lambdas[index[0]] = 0
        
    # otherwise, your FOI changes based on your colonization history and competition    
    else:
        ind_lambdas[index[0]] = (1-(sig*(prev_col[index[0]] > 0)))*ind_lambdas[index[0]]*(1-thet*(susc < max_carr_cap))   
        
    # if individual is vaccinated, VT not carried have rate = 1/vacc_dur, overwrites ind_lambdas = 0 from previous step
    if vacc == True:
        ind_lambdas[vt_ind] = 1/vacc_dur
    
    # serotypes you are carrying change recovery rate based on colonization history and immunity
    #ind_lambdas[non_index] = ind_lambdas[non_index]*math.exp(eps*sum(prev_col)) #(1/250) + ((ind_lambdas[non_index]-(1/250))*math.exp(eps*sum(prev_col)))
    ind_lambdas[non_index] = 1/(mincarrdur+(((1/ind_lambdas[non_index])-mincarrdur)*math.exp(-eps*sum(prev_col))))
    #if ind_lambdas[non_index] < (1/300): # max carriage duration is 300 days
    #    ind_lambdas[non_index] = 1/300
    #if np.any(ind_lambdas[non_index] > (1/5)): # min carriage duration, max recovery rate
    #    ind_lambdas[non_index] = 1/5
    return(ind_lambdas)

# import fitted carriage duration params from IBM_paramfitting.py
seroparams = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/seroparams.csv")
# DATA FROM 20 YR STEADY STATE SIMULATIIION
# carriers = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/IBM sims/tot population/carriers_fullpop_migr_20yrs_0111.csv")
# carriers = carriers.drop(['Unnamed: 0'], axis = 1)
# carriers.vacc = False
# prev_col = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/IBM sims/tot population/prevcol_fullpop_migr_20yrs_0111.csv")
# prev_col = prev_col.drop(['Unnamed: 0'], axis = 1)
# ind_lambdas = prev_col.copy(deep = True)
# output = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/IBM sims/tot population/output_fullpop_migr_20yrs_0111.csv")
# output.columns = ['key', 'value']
# sero_carr = output.drop([0, 72])
# t = 7299
#sero_carr.value.sort()
colnames = np.unique(seroparams['serotype']) #list(map(lambda x: 'serotype' + str(x+1), range(num_sero)))

### HIGHEST PREVALENT SEROTYPES (IN DISEASE) FOR VACCINATION
# COMMENTED OUT - top in carriage
#top13 = np.array(sero_carr.value)[np.argsort(sero_carr.value)][-13:]
#top_ind = np.isin(sero_carr.value, top13)
#sero_carr.loc[top_ind]

# TOP IN DISEASE
disinc = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/IBM sims/simulation analysis/20yrsimulation_diseaseincidence_0511.csv")
top13 = np.array(disinc.disease)[np.argsort(disinc.disease)][-13:]
top_ind = np.isin(disinc.disease, top13)
vt = disinc.loc[top_ind].Serotype
vt_ind = np.where(np.isin(colnames, vt))

# Initialise population
indivs = 1000
num_sero = len(np.unique(seroparams['serotype']))
max_carr_cap = 1
pop_breakdown = [0.2, 0.6, 0.2] # has to sum to 1
pop_child = int(pop_breakdown[0]*indivs)
pop_adult = int(pop_breakdown[1]*indivs)
pop_elder = int(pop_breakdown[2]*indivs)

# (physical AND non physical) contact matrix taken from Mossong et al Finland: https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050074#pmed-0050074-st005
contact_matrix = np.array([[21.21, 12.42, 0.85], # number contacts per day kids make with [0] kids, [1] adults, and [2] elderly
                  [16.44, 91.40, 6.42], # number contacts per day adults make with [0] kids, [1] adults, and [2] elderly
                  [1.04, 4.95, 4.15]]) # number contacts per day elderly make with [0] kids, [1] adults, and [2] elderly


full_pop = np.empty([indivs, num_sero])
for i in range(num_sero):
    full_pop[:,i] = np.array([np.random.choice([int(0),int(1)], size = indivs, replace = True, p = [0.99, 0.01])])
  
carriers_all = pd.DataFrame(data = full_pop, columns = colnames)
carriers_all = carriers_all.astype(int) # change floating 0s and 1s to integers
prev_col = carriers_all.copy(deep = True) # serotype-specific previous colonizations
ind_lambdas = carriers_all.copy(deep = True) # frame for individual lambda rates 
susc_full = max_carr_cap - carriers_all.sum(axis = 1)
carriers_all['susc'] = susc_full
agegroup = np.repeat([int(0), int(1), int(2)], [pop_child, pop_adult, pop_elder])
carriers_all['agegroup'] = agegroup
age = random.sample(range(0,6205), pop_child)+random.sample(range(6206,23360), pop_adult)+random.sample(range(23361,36500), pop_elder)
carriers_all['age'] = age
# # comment out for no vaccination
carriers_all['vacc'] = False

# beta (transmission or FOI) rates ## these will be different for children and adults but are the same for all serotypes
# values taken from Melegaro 2004 longitudinal carriage study in the UK
# Finland transmission param 0.63/month = .021/day https://bmcinfectdis.biomedcentral.com/articles/10.1186/1471-2334-13-180 
B_rates_child_val = 0.021#[0.012/indivs]*num_sero # per day
B_rates_adult_val = 0.0063 # children are 70% of FOI so adult FOI is reduced by this amount #[0.004/indivs]*num_sero # per day
B_rates_elder_val = 0.021#[0.012/indivs]*num_sero # per day
seroparams['beta_child'] = B_rates_child_val
seroparams['beta_adult'] = B_rates_adult_val
seroparams['beta_elder'] = B_rates_elder_val

#### IMMUNITY AND COMPETITION PARAMS 
sigma = 0.4 # specific immunity
epsilon = 0.1 # non-specific immunity
theta = 0.01 # competition
mincarrdur = 5 # minimum carriage duration

# real VT - NOT used
# vt7 = ['4','6B', '9V','14','18C','19F','23F']
# vt10 = ['1','5','7F']
# vt13 = ['3','6A', '19A']
# vt7_ind = np.where(np.isin(colnames, vt7))
# vt10_ind = np.where(np.isin(colnames, vt10))
# vt13_ind = np.where(np.isin(colnames, vt13))

# vaccine params
vacc_start = 7300 # start of vaccination period
vacc_dur = 10 # years of vaccine-induced protection
vacc_dur = vacc_dur*365 # protection converted to days 
vacc_eff = 1

#carriers = carriers_all.copy(deep = True)
#t = 0
stop_t = 9125 #9125 # 7300 = 20 years; 9125 = 5 years post vaccination at 20 yrs
# initial dummy tot_lambda
tot_lambda1 = (sum(carriers.iloc[:,0])*seroparams['beta'][0])
    
# outputs for graph of time vs infected ppl
chunkSize = 100000
currentLength = chunkSize
DS_inf_nams = np.unique(seroparams['serotype'])#list(map(lambda x: 'infxtd' + str(x+1), range(num_sero)))
plotsy = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = DS_inf_nams)
plotsy.insert(loc = 0, column = 'time', value = np.nan)
    
next_ev = 1
start = time.time()
while t <= stop_t and tot_lambda1 != 0:
    
    if len(np.where(carriers.age > 36500)[0]) > 0:
        ind_to_die = np.where(carriers.age > 36500)[0]
        carriers.loc[ind_to_die] = 0 # making this person a baby now
        carriers.loc[ind_to_die, 'susc'] = max_carr_cap # susc = max carrying capacity when nothing is carried
        carriers.loc[ind_to_die, 'vacc'] = False
        prev_col.loc[ind_to_die] = 0
        
    index_children = np.where(carriers.age <= 6205)[0]
    index_adults = np.where((carriers.age > 6205) & (carriers.age <= 23360))[0]
    index_elderly = np.where(carriers.age > 23360)[0]
    
    carriers.loc[index_children, 'agegroup'] = 0
    carriers.loc[index_adults, 'agegroup'] = 1
    carriers.loc[index_elderly, 'agegroup'] = 2
        
    # add NAs to end of output dataframe if there is no space left to fill
    if next_ev > currentLength:
        plotsy = plotsy.append(pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = plotsy.columns), ignore_index = True)
        currentLength = currentLength + chunkSize
    
    # only record data when t > 9000 (i.e. last 3 years of time)
    #if t > 9000: 100 # only record after 100th day
    plotsy.iloc[next_ev-1, 1:] = np.array(carriers.iloc[:,0:num_sero].sum(axis = 0))
    plotsy.time[next_ev-1] = t
    
    if t >= 7300:
        index_vacc = np.where((carriers.age >= 730) & (carriers.age <= 1095) & (carriers.vacc == False)) # index of ppl to be vaccinated
        num_newvacc = round(vacc_eff*len(index_vacc[0])) # number of new people to vaccinate (= vacc efficacy * tot number of people)
        adjusted_ind = index_vacc[0][0:num_newvacc] # indices of people for whom vaccination works/is protective
        carriers.loc[adjusted_ind, 'vacc'] = True
    
    #if t > 3650:
    #    carriers_all.vacc[(carriers_all.age >= 730) & (carriers_all.age <= 1095)] = True
        
    #migration: add 1 carrier to serotypes with zero carriers
    extinctsero_child = np.where(carriers.iloc[index_children,0:num_sero].sum(axis = 0) == 0)[0]
    totchildcarr = carriers.iloc[index_children,0:num_sero].sum(axis = 0)
    totchildcarr[extinctsero_child] = 1
    
    extinctsero_adult = np.where(carriers.iloc[index_adults,0:num_sero].sum(axis = 0) == 0)[0]
    totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
    totadultcarr[extinctsero_adult] = 1
    
    extinctsero_elder = np.where(carriers.iloc[index_elderly,0:num_sero].sum(axis = 0) == 0)[0]
    toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)
    toteldercarr[extinctsero_elder] = 1
    
    totFOIbysero_children = [0]*num_sero
    totFOIbysero_adults = [0]*num_sero
    totFOIbysero_elders = [0]*num_sero
    
    for i in range(len(seroparams)):
        totFOIbysero_children[i] = (contact_matrix[0][0]*seroparams['beta_child'][i]*totchildcarr[i]+ 
                                    contact_matrix[0][1]*seroparams['beta_adult'][i]*totadultcarr[i]+
                                    contact_matrix[0][2]*seroparams['beta_elder'][i]*toteldercarr[i])
        totFOIbysero_adults[i] = (contact_matrix[1][0]*seroparams['beta_child'][i]*totchildcarr[i]+  
                                  contact_matrix[1][1]*seroparams['beta_adult'][i]*totadultcarr[i]+ 
                                  contact_matrix[1][2]*seroparams['beta_elder'][i]*toteldercarr[i])
        totFOIbysero_elders[i] =  (contact_matrix[2][0]*seroparams['beta_child'][i]*totchildcarr[i]+ 
                                    contact_matrix[2][1]*seroparams['beta_adult'][i]*totadultcarr[i]+ 
                                    contact_matrix[2][2]*seroparams['beta_elder'][i]*toteldercarr[i])
   
    totFOIbysero = np.row_stack((totFOIbysero_children, totFOIbysero_adults, totFOIbysero_elders))
    
    #totalFOI_children = (contact_matrix[0][0]*B_rates_child_val*totchildcarr +
    #                     contact_matrix[0][1]*B_rates_adult_val*totadultcarr +
    #                     contact_matrix[0][2]*B_rates_elder_val*toteldercarr)
    #totalFOI_adults = (contact_matrix[1][0]*B_rates_child_val*totchildcarr + 
    #                   contact_matrix[1][1]*B_rates_adult_val*totadultcarr +
    #                   contact_matrix[1][2]*B_rates_elder_val*toteldercarr)
    #totalFOI_elderly = (contact_matrix[2][0]*B_rates_child_val*totchildcarr + 
    #                    contact_matrix[2][1]*B_rates_adult_val*totadultcarr +
    #                    contact_matrix[2][2]*B_rates_elder_val*toteldercarr)
    #totFOIbysero = np.row_stack((totalFOI_children, totalFOI_adults, totalFOI_elderly)) 
    
    # individual rates                
    for s in range(indivs):
        for j in range(num_sero):
            prevcol = prev_col.iloc[s,j]
            alpha = seroparams['alpha_emcee'][j]
            A = seroparams['A_emcee'][j]
            B = seroparams['B_emcee'][j]
            # if person not carrying
            if carriers.iloc[s,j] == 0:
                # if person is vaccinated & serotype is VT
                if carriers.vacc[s] == True and np.isin(j,vt_ind):
                    ind_lambdas.iloc[s,j] = 1/vacc_dur
                # if serotype is NVT or if person is not vaccinated
                else:
                    ind_lambdas.iloc[s,j] = totFOIbysero[carriers['agegroup'][s],j]
            # if person is carrying, then recovery rate as function of prevcol
            else:
                ind_lambdas.iloc[s,j] = (1/np.random.gamma(shape = alpha, scale = (A*np.exp(-prevcol*B))/alpha))                
    
    # tot.lambda = total rate of all/any events
    i = 0
    U0 = np.random.uniform(0,1, size = 20)
    U3 = np.random.uniform(0,1, size = 20)
    while i < 20: # update 50 individuals each time step
        tot_lambda1 = ind_lambdas.sum().sum()
        #print(tot_lambda1)
        U1 = U0[i]
        U1_tot_lambda1 = U1*tot_lambda1
        
        # choose person index ( = which person will this happen to?)
        cumsum_person = ind_lambdas.sum(axis = 1)
        which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        #try:
        #    which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        #except IndexError:
        #    which_person = np.where(np.random.uniform(0,1)*ind_lambdas.sum().sum() < np.cumsum(cumsum_person))[0][0]#which_person = indivs - 1 #
        #    break
        #print(which_person)
        
        # choose event that will take place (= which serotype will be toggled?)
        person_lambdas = new_lambdas(carriers.loc[which_person], ind_lambdas.loc[which_person], prev_col.loc[which_person], sigma)
            
        if np.any(person_lambdas < 0):
            print('NEGATIVE IND_LAMBDA')
            print("indiv:", which_person)
            print("carr status:", carriers.iloc[which_person, np.where(person_lambdas < 0)[0][0]])
            print("serotype:", colnames[np.where(person_lambdas < 0)])
            print("prevcol:",prev_col.iloc[which_person, np.where(person_lambdas < 0)[0][0]])
            break
        tot_lambda2 = sum(person_lambdas)
        U2 = U3[i]
        U2_tot_lambda2 = U2*tot_lambda2

        which_indices = np.where(U2_tot_lambda2 < np.cumsum(person_lambdas))
        
        #midtime6 = time.time()
        
        # switching event (S --> I or I --> S)
        if which_indices[0].size > 0:
            which_sero = which_indices[0][0]
            carriers.iloc[which_person, which_sero] = ~carriers.iloc[which_person, which_sero] + 2        

        #print(which_sero)
        # reupdate the susc score
        new_susc = max_carr_cap - sum(carriers.iloc[which_person, 0:num_sero])
    
        # update colonisation history if person acquired new serotype    
        if new_susc < carriers.susc.loc[which_person]:
            prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        i += 1

    # update time and age
    t += 1
    next_ev += 1
    carriers.age += 1
        
end = time.time()
total_time = end - start
plotsy = plotsy.dropna()
plotsy_melt= pd.melt(plotsy, id_vars = 'time')
plot = plotsy_melt.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

# output file
sim_output = plotsy.loc[len(plotsy) - 1]
sim_output['sim_time'] = total_time/3600

mylist = [f for f in glob.glob("//qdrive/homes/al817/Technical/Python/simulation*.csv")]
latestfile = max(mylist)
simnum = int(latestfile.rsplit('simulationrun', 1)[1].rsplit('.csv',1)[0])+1
newfiledest = "//qdrive/homes/al817/Technical/Python/simulationrun" + str(simnum) + ".csv"
sim_output.to_csv(path_or_buf= newfiledest, index = True)

# prevcol
prevcolnam = "//qdrive/homes/al817/Technical/Python/prevcol" + str(simnum) + ".csv"
prev_col.to_csv(path_or_buf = prevcolnam, index = True)

# carriers
carrnam =  "//qdrive/homes/al817/Technical/Python/carriers_tmax" + str(simnum) + ".csv"
carriers.to_csv(path_or_buf = carrnam, index = True)

# output figure
plotname = "simulation" + str(simnum) + ".pdf"
#plot.figure.savefig(plotname, bbox_inches = "tight")
#plt.plot(plotsy_melt).savefig(plotname, bbox_inches = "tight")

# comment this out for cluster
#plot.get_figure().savefig(plotname, bbox_inches = "tight")

# POST-VACCINATION

vacc_stop_t <- 4000

while t <= stop_t and tot_lambda1 != 0:
    
    if len(np.where(carriers.age > 36500)[0]) > 0:
        ind_to_die = np.where(carriers.age > 36500)[0]
        carriers.loc[ind_to_die] = 0 # making this person a baby now
        carriers.loc[ind_to_die, 'susc'] = max_carr_cap # susc = max carrying capacity when nothing is carried
        carriers.loc[ind_to_die, 'vacc'] = False
        prev_col.loc[ind_to_die] = 0
        
    index_children = np.where(carriers.age <= 6205)[0]
    index_adults = np.where((carriers.age > 6205) & (carriers.age <= 23360))[0]
    index_elderly = np.where(carriers.age > 23360)[0]
    
    carriers.loc[index_children, 'agegroup'] = 0
    carriers.loc[index_adults, 'agegroup'] = 1
    carriers.loc[index_elderly, 'agegroup'] = 2
        
    # add NAs to end of output dataframe if there is no space left to fill
    if next_ev > currentLength:
        plotsy = plotsy.append(pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = plotsy.columns), ignore_index = True)
        currentLength = currentLength + chunkSize
    
    # only record data when t > 9000 (i.e. last 3 years of time)
    #if t > 9000: 100 # only record after 100th day
    plotsy.iloc[next_ev-1, 1:] = np.array(carriers.iloc[:,0:num_sero].sum(axis = 0))
    plotsy.time[next_ev-1] = t
    
    if t >= 7300:
        index_vacc = np.where((carriers.age >= 730) & (carriers.age <= 1095) & (carriers.vacc == False)) # index of ppl to be vaccinated
        num_newvacc = round(vacc_eff*len(index_vacc[0])) # number of new people to vaccinate (= vacc efficacy * tot number of people)
        adjusted_ind = index_vacc[0][0:num_newvacc] # indices of people for whom vaccination works/is protective
        carriers.loc[adjusted_ind, 'vacc'] = True
    
    #if t > 3650:
    #    carriers_all.vacc[(carriers_all.age >= 730) & (carriers_all.age <= 1095)] = True
        
    #migration: add 1 carrier to serotypes with zero carriers
    extinctsero_child = np.where(carriers.iloc[index_children,0:num_sero].sum(axis = 0) == 0)[0]
    totchildcarr = carriers.iloc[index_children,0:num_sero].sum(axis = 0)
    totchildcarr[extinctsero_child] = 1
    
    extinctsero_adult = np.where(carriers.iloc[index_adults,0:num_sero].sum(axis = 0) == 0)[0]
    totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
    totadultcarr[extinctsero_adult] = 1
    
    extinctsero_elder = np.where(carriers.iloc[index_elderly,0:num_sero].sum(axis = 0) == 0)[0]
    toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)
    toteldercarr[extinctsero_elder] = 1
    
    totFOIbysero_children = [0]*num_sero
    totFOIbysero_adults = [0]*num_sero
    totFOIbysero_elders = [0]*num_sero
    
    for i in range(len(seroparams)):
        totFOIbysero_children[i] = (contact_matrix[0][0]*seroparams['beta_child'][i]*totchildcarr[i]+ 
                                    contact_matrix[0][1]*seroparams['beta_adult'][i]*totadultcarr[i]+
                                    contact_matrix[0][2]*seroparams['beta_elder'][i]*toteldercarr[i])
        totFOIbysero_adults[i] = (contact_matrix[1][0]*seroparams['beta_child'][i]*totchildcarr[i]+  
                                  contact_matrix[1][1]*seroparams['beta_adult'][i]*totadultcarr[i]+ 
                                  contact_matrix[1][2]*seroparams['beta_elder'][i]*toteldercarr[i])
        totFOIbysero_elders[i] =  (contact_matrix[2][0]*seroparams['beta_child'][i]*totchildcarr[i]+ 
                                    contact_matrix[2][1]*seroparams['beta_adult'][i]*totadultcarr[i]+ 
                                    contact_matrix[2][2]*seroparams['beta_elder'][i]*toteldercarr[i])
   
    totFOIbysero = np.row_stack((totFOIbysero_children, totFOIbysero_adults, totFOIbysero_elders))
    
    #totalFOI_children = (contact_matrix[0][0]*B_rates_child_val*totchildcarr +
    #                     contact_matrix[0][1]*B_rates_adult_val*totadultcarr +
    #                     contact_matrix[0][2]*B_rates_elder_val*toteldercarr)
    #totalFOI_adults = (contact_matrix[1][0]*B_rates_child_val*totchildcarr + 
    #                   contact_matrix[1][1]*B_rates_adult_val*totadultcarr +
    #                   contact_matrix[1][2]*B_rates_elder_val*toteldercarr)
    #totalFOI_elderly = (contact_matrix[2][0]*B_rates_child_val*totchildcarr + 
    #                    contact_matrix[2][1]*B_rates_adult_val*totadultcarr +
    #                    contact_matrix[2][2]*B_rates_elder_val*toteldercarr)
    #totFOIbysero = np.row_stack((totalFOI_children, totalFOI_adults, totalFOI_elderly)) 
    
    # individual rates                
    for s in range(indivs):
        for j in range(num_sero):
            prevcol = prev_col.iloc[s,j]
            alpha = seroparams['alpha_emcee'][j]
            A = seroparams['A_emcee'][j]
            B = seroparams['B_emcee'][j]
            # if person not carrying
            if carriers.iloc[s,j] == 0:
                # if person is vaccinated & serotype is VT
                if carriers.vacc[s] == True and np.isin(j,vt_ind):
                    ind_lambdas.iloc[s,j] = 1/vacc_dur
                # if serotype is NVT or if person is not vaccinated
                else:
                    ind_lambdas.iloc[s,j] = totFOIbysero[carriers['agegroup'][s],j]
            # if person is carrying, then recovery rate as function of prevcol
            else:
                ind_lambdas.iloc[s,j] = (1/np.random.gamma(shape = alpha, scale = (A*np.exp(-prevcol*B))/alpha))                
    
    # tot.lambda = total rate of all/any events
    i = 0
    U0 = np.random.uniform(0,1, size = 20)
    U3 = np.random.uniform(0,1, size = 20)
    while i < 20: # update 50 individuals each time step
        tot_lambda1 = ind_lambdas.sum().sum()
        #print(tot_lambda1)
        U1 = U0[i]
        U1_tot_lambda1 = U1*tot_lambda1
        
        # choose person index ( = which person will this happen to?)
        cumsum_person = ind_lambdas.sum(axis = 1)
        which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        #try:
        #    which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        #except IndexError:
        #    which_person = np.where(np.random.uniform(0,1)*ind_lambdas.sum().sum() < np.cumsum(cumsum_person))[0][0]#which_person = indivs - 1 #
        #    break
        #print(which_person)
        
        # choose event that will take place (= which serotype will be toggled?)
        person_lambdas = new_lambdas(carriers.loc[which_person], ind_lambdas.loc[which_person], prev_col.loc[which_person], sigma)
            
        if np.any(person_lambdas < 0):
            print('NEGATIVE IND_LAMBDA')
            print("indiv:", which_person)
            print("carr status:", carriers.iloc[which_person, np.where(person_lambdas < 0)[0][0]])
            print("serotype:", colnames[np.where(person_lambdas < 0)])
            print("prevcol:",prev_col.iloc[which_person, np.where(person_lambdas < 0)[0][0]])
            break
        tot_lambda2 = sum(person_lambdas)
        U2 = U3[i]
        U2_tot_lambda2 = U2*tot_lambda2

        which_indices = np.where(U2_tot_lambda2 < np.cumsum(person_lambdas))
        
        #midtime6 = time.time()
        
        # switching event (S --> I or I --> S)
        if which_indices[0].size > 0:
            which_sero = which_indices[0][0]
            carriers.iloc[which_person, which_sero] = ~carriers.iloc[which_person, which_sero] + 2        

        #print(which_sero)
        # reupdate the susc score
        new_susc = max_carr_cap - sum(carriers.iloc[which_person, 0:num_sero])
    
        # update colonisation history if person acquired new serotype    
        if new_susc < carriers.susc.loc[which_person]:
            prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        i += 1

    # update time and age
    t += 1
    next_ev += 1
    carriers.age += 1



####### Testing with dummy carriage data ########################################################################################################################

# taken from this study https://academic.oup.com/jid/article/184/4/451/809030 
#serogroup = np.array(['6B', '23F', '19F', '6A', '11', '14', '15', '35', 'R', 'NT'])
#num_isolates = np.array([213, 268, 258, 140, 116, 67, 56, 49, 55, 17])
#list_of_tuples = list(zip(serogroup, num_isolates))
#carrprevdf = pd.DataFrame(list_of_tuples, columns = ['Serotype', 'Carr_isolates'])
#carrprevdf['Carr_prev'] = carrprevdf.Carr_isolates/1530 #1530 is the total number of swabs from the study

###############################################################################################################################################################

