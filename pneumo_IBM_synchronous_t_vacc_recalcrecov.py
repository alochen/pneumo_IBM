########## PNEUMO IBM with fitted sero params imported and with synchronous time where 20 events occur at each ############
########## time step, and a migration rate such that when a serotype goes extinct, it is reseeded into the     ############
########## population in the form of one infectious carrier                                                    ############

import pandas as pd
import numpy as np
import math
import random
import time
import glob
import os
from builtins import dict

#import matplotlib
#import matplotlib.pyplot as plt

# import parameters
modparams = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/modparams.csv")
country = modparams.country[0]
globalinv = False
indivs = modparams.indivs[0]
sigma = modparams.sigma[0] # specific immunity
epsilon = modparams.epsilon[0] # non-specific immunity
theta = modparams.theta[0] # competition
mincarrdur = 5 # minimum carriage duration
vacc_dur = modparams.vacc_dur[0] # years of vaccine-induced protection
vacc_eff = modparams.vacc_eff[0] # vaccine efficacy
beta = modparams.beta[0] # transmission coefficient

# ## Scenario analysis params
# country = "fra"
# globalinv = False
# indivs = 1000
# sigma = 0.6 # specific immunity
# epsilon = 0.25 # non-specific immunity
# theta = 0.1 #0.01 # competition
# mincarrdur = 5 # minimum carriage duration
# vacc_dur = 10 # years of vaccine-induced protection
# vacc_eff = 0.95 # vaccine efficacy
# beta = 0.00012075 # transmission coefficient

#os.getcwd()
#os.chdir('Q:\\Technical\\Python')
phase = 'prevacc'

filedest = "//qdrive/homes/al817/Technical/Python/popruns/new May 2021 log files/" + country + '_dis/' 

seroparams = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/seroparams.csv")
seroparams = seroparams.drop(seroparams.index[seroparams.serotype.isin(['15B/C_old','15B','15C'])]).reset_index()
colnames = np.unique(seroparams['serotype']) #list(map(lambda x: 'serotype' + str(x+1), range(num_sero)))

mylist = [f for f in glob.glob(filedest + "*.csv")]
mylist = [fname.rsplit('.', 1)[0] for fname in os.listdir(filedest)]
filenums = [item.split('_',1)[0] for item in mylist]
filenums = [item.split('simulation',1)[-1] for item in filenums]
filenums = [int(i) for i in filenums] 
simnum = max(filenums) + 1
simnum = str(simnum) + "_pop" + str(indivs)
endstring = "simulation" + str(simnum) + "_dis_" + phase + "_" + country + "_" + 'initFOIbeta_'
seroparams.to_csv(path_or_buf = filedest + endstring + 'dummy.csv', index = False)

if country == "usa":
    country_inv = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/usalocalinv.csv")
if country == "fi":
    country_inv = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/finlandlocalinv.csv")
else:
    country_inv = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/fralocalinv.csv")

#################################################################################################################################################################
############# PART 1: PRE-VACCINATION PERIOD ####################################################################################################################
#################################################################################################################################################################

# function that computes new lambdas based on serotype switching event
def new_lambdas(carrier_person, ind_lambdas, prev_col, sigma):    
    susc = carrier_person.susc
    vacc = carrier_person.vacc
    
    index = np.where(carrier_person[0:num_sero] == 0) # index of serotypes not carried
    non_index = np.bool_(carrier_person[0:num_sero]) # index of serotypes carried
    
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

    # serotypes you are carrying change recovery rate based on colonization history and immunity
    ind_lambdas[non_index] = 1/(mincarrdur+(((1/ind_lambdas[non_index])-mincarrdur)*math.exp(-eps*sum(prev_col))))
    return(ind_lambdas)

# Initialise population
num_sero = len(np.unique(seroparams['serotype']))
max_carr_cap = 2 #1

# CONTACT MATRICES for children < 15yrs
# # FI (physical AND non physical) contact matrix taken from Mossong et al Finland: https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050074#pmed-0050074-st005
# contact_matrix_fi = np.array([[21.21, 12.42, 0.85], # avg number contacts per day kids make with [0] kids, [1] adults, and [2] elderly
#                   [16.44, 91.40, 6.42], # avg number contacts per day adults make with [0] kids, [1] adults, and [2] elderly
#                   [1.04, 4.95, 4.15]]) # avg number contacts per day elderly make with [0] kids, [1] adults, and [2] elderly
# 
# # USA (physical AND non physical) contact matrix taken from Prem et al: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697#sec002
# contact_matrix_usa = np.array([[22.2773775, 15.13427685, 0.67475806], # number contacts per day kids make with [0] kids, [1] adults, and [2] elderly
#                             [19.65830386, 121.0409943, 2.531342454], # number contacts per day adults make with [0] kids, [1] adults, and [2] elderly
#                             [2.730939972, 10.13063602, 3.481243491]])# number contacts per day elderly make with [0] kids, [1] adults, and [2] elderly
# 
# # France (physical AND non physical) contact matrix taken from Prem et al: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697#sec002
# contact_matrix_fr = np.array([[23.9255104, 28.52664938, 0.833780374], # number contacts per day kids make with [0] kids, [1] adults, and [2] elderly
#                            [14.17373913, 102.3279637, 2.785170733],# number contacts per day adults make with [0] kids, [1] adults, and [2] elderly
#                            [2.488764893, 9.938006315, 4.714598106]])# number contacts per day elderly make with [0] kids, [1] adults, and [2] elderly

# CONTACT MATRICES for children < 5yrs, all taken from Prem et al (physical and nonphysical) https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697#sec002
contact_matrix_fi = np.array([[4.155665732, 6.254676868, 0.220187955],
                              [5.689855384, 151.7515979, 3.152128926],
                              [0.851313659, 7.667759404, 6.515933663]])

contact_matrix_fr = np.array([[3.804920797, 6.691382506, 0.26971156],
                              [5.489241519, 152.9683178, 3.349239548],
                              [0.546826761, 12.75265203, 4.714598105]])

contact_matrix_usa = np.array([[2.598237474, 6.686712073, 0.248696214],
                               [5.675060688, 163.1509422, 2.957404299],
                               [0.59186623, 12.26970976, 3.481243491]])


if country == 'usa':
    contact_matrix = contact_matrix_usa
    maxage = 78.5 # life expectancy USA 78.5 years (2010) https://datacatalog.worldbank.org/dataset/world-development-indicators 
    pop_breakdown = [0.06, 0.81, 0.13]
if country == 'fi':
    contact_matrix = contact_matrix_fi
    maxage = 79.9 # life expectancy 79.9 years (2010) https://datacatalog.worldbank.org/dataset/world-development-indicators
    pop_breakdown = [0.05, 0.77, 0.18]
else:
    contact_matrix = contact_matrix_fr
    maxage = 81.7 # life expectancy 81.7 years (2010) https://datacatalog.worldbank.org/dataset/world-development-indicators 
    pop_breakdown = [0.06, 0.77, 0.17]

pop_child = int(pop_breakdown[0]*indivs)
pop_adult = int(pop_breakdown[1]*indivs)
pop_elder = int(pop_breakdown[2]*indivs)
agecutoffs = [5, 65, maxage]
age_breakdown = [round(i * 365) for i in agecutoffs] #[365*5, 23725, maxage]

full_pop = np.empty([indivs, num_sero])
for i in range(num_sero):
    full_pop[:,i] = np.array([np.random.choice([int(0),int(1)], size = indivs, replace = True, p = [0.999, 0.001])])

carriers_all = pd.DataFrame(data = full_pop, columns = colnames).astype(np.int8)
#carriers_all1 = carriers_all.astype(int) # change floating 0s and 1s to integers
prev_col = carriers_all.copy(deep = True) # serotype-specific previous colonizations
ind_lambdas = carriers_all.copy(deep = True) # frame for individual lambda rates 
#ind_lambdas = ind_lambdas.to_numpy(dtype = np.float16)
susc_full = max_carr_cap - carriers_all.sum(axis = 1)
carriers_all['susc'] = susc_full
agegroup = np.repeat([int(0), int(1), int(2)], [pop_child, pop_adult, pop_elder])
carriers_all['agegroup'] = agegroup
age = random.choices(range(0,age_breakdown[0]), k=pop_child)+random.choices(range(age_breakdown[0],age_breakdown[1]), k=pop_adult)+random.choices(range(age_breakdown[1],age_breakdown[2]), k=pop_elder)
#np.concatenate((np.random.choice(range(0, 6206), size = pop_child, replace = True), np.random.choice(range(6206,23361), size = pop_adult, replace = True),np.random.choice(range(23361,36500), size = pop_elder, replace = True)), axis=None)
carriers_all['age'] = age

# change previous colonisations for adults and elderly according to child to adult simulations' distributions
prevcolprob = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/Prevcolhist/prevcolprob_allsims.csv")
prevcolprob.loc[np.where(prevcolprob.serotype == '15B.C')[0], 'serotype'] = '15B/C'
prevcolprob.loc[np.where(prevcolprob.serotype == '6A.C')[0], 'serotype'] = '6A/C'
for i in colnames:
    df = prevcolprob[prevcolprob.serotype == i]
    #prev_col.loc[np.where((carriers_all['age'] >= 6205) & (carriers_all['age'] < 9125))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 6205) & (carriers_all['age'] < 9125))[0]), replace = True, p = df.m20) 
    prev_col.loc[np.where((carriers_all['age'] >= 1825) & (carriers_all['age'] < 9125))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 1825) & (carriers_all['age'] < 9125))[0]), replace = True, p = df.m20) 
    prev_col.loc[np.where((carriers_all['age'] >= 9125) & (carriers_all['age'] < 10950))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 9125) & (carriers_all['age'] < 10950))[0]), replace = True, p = df.m25)
    prev_col.loc[np.where((carriers_all['age'] >= 10950) & (carriers_all['age'] < 12775))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 10950) & (carriers_all['age'] < 12775))[0]), replace = True, p = df.m30)
    prev_col.loc[np.where((carriers_all['age'] >= 12775) & (carriers_all['age'] < 14600))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 12775) & (carriers_all['age'] < 14600))[0]), replace = True, p = df.m35)
    prev_col.loc[np.where((carriers_all['age'] >= 14600) & (carriers_all['age'] < 16425))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 14600) & (carriers_all['age'] < 16425))[0]), replace = True, p = df.m40)
    prev_col.loc[np.where((carriers_all['age'] >= 16425) & (carriers_all['age'] < 18250))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 16425) & (carriers_all['age'] < 18250))[0]), replace = True, p = df.m45)
    prev_col.loc[np.where((carriers_all['age'] >= 18250) & (carriers_all['age'] < 20075))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 18250) & (carriers_all['age'] < 20075))[0]), replace = True, p = df.m50)
    prev_col.loc[np.where((carriers_all['age'] >= 20075) & (carriers_all['age'] < 21900))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 20075) & (carriers_all['age'] < 21900))[0]), replace = True, p = df.m55)
    prev_col.loc[np.where((carriers_all['age'] >= 21900) & (carriers_all['age'] < 23725))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 21900) & (carriers_all['age'] < 23725))[0]), replace = True, p = df.m60)
    prev_col.loc[np.where((carriers_all['age'] >= 23725) & (carriers_all['age'] < 25550))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 23725) & (carriers_all['age'] < 25550))[0]), replace = True, p = df.m65)
    prev_col.loc[np.where((carriers_all['age'] >= 25550) & (carriers_all['age'] < 27375))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 25550) & (carriers_all['age'] < 27375))[0]), replace = True, p = df.m70)
    prev_col.loc[np.where(carriers_all['age'] >= 27375)[0], i] = np.random.choice(df.prevcol, size = len(np.where(carriers_all['age'] >= 27375)[0]), replace = True, p = df.m75) #) & (carriers_all['age'] < 29200) ... ) & (carriers_all['age'] < 29200)
    #prev_col.loc[np.where((carriers_all['age'] >= 29200) & (carriers_all['age'] < 31025))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 29200) & (carriers_all['age'] < 31025))[0]), replace = True, p = df.m80)
    #prev_col.loc[np.where((carriers_all['age'] >= 31025) & (carriers_all['age'] < 32850))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 31025) & (carriers_all['age'] < 32850))[0]), replace = True, p = df.m85)
    #prev_col.loc[np.where(carriers_all['age'] >= 32850)[0], i] = np.random.choice(df.prevcol, size = len(np.where(carriers_all['age'] >= 32850)[0]), replace = True, p = df.m90)


carriers_all['vacc'] = False
carriers = carriers_all.copy(deep = True)

# beta (transmission or FOI) rates ## these will be different for children and adults but are the same for all serotypes
# values taken from Melegaro 2004 longitudinal carriage study in the UK https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2870123/pdf/15188713.pdf
# Finland transmission param 0.63/month = .021/day https://bmcinfectdis.biomedcentral.com/articles/10.1186/1471-2334-13-180 
B_rates_child_val = 0.012/indivs #0.021#[0.012/indivs]*num_sero # per day
B_rates_adult_val = 0.004/indivs #0.0063 # children are 70% of FOI so adult FOI is reduced by this amount #[0.004/indivs]*num_sero # per day
B_rates_elder_val = 0.012/indivs #0.021#[0.012/indivs]*num_sero # per day
seroparams['beta_child'] = beta#0.00010169 #0.012 #0.021 #B_rates_child_val
seroparams['beta_adult'] = beta#0.00010169 #0.012 #0.021 #B_rates_adult_val
seroparams['beta_elder'] = beta#0.00010169 #0.012 #0.021 #B_rates_elder_val

# index of children and sum of all children, adult and elderly carriers
index_children = np.where(carriers.age < age_breakdown[0])[0]
index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
index_elderly = np.where(carriers.age >= age_breakdown[1])[0]

totchildcarr = carriers.iloc[index_children,0:num_sero].sum(axis = 0)
totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)

# calculation of initial rates of individuals
totFOIbysero_children = [0]*num_sero
totFOIbysero_adults = [0]*num_sero
totFOIbysero_elders = [0]*num_sero
    
for i in range(len(seroparams)):
    totFOIbysero_children[i] = ((contact_matrix[0][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                (contact_matrix[0][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+
                                (contact_matrix[0][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
    totFOIbysero_adults[i] = ((contact_matrix[1][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                              (contact_matrix[1][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                              (contact_matrix[1][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
    totFOIbysero_elders[i] =  ((contact_matrix[2][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                               (contact_matrix[2][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                               (contact_matrix[2][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
   
totFOIbysero = np.row_stack((totFOIbysero_children, totFOIbysero_adults, totFOIbysero_elders))

for s in range(indivs):
    for j in range(num_sero):
        prevcol = prev_col.iloc[s,j]
        alpha = seroparams['alpha_emcee'][j]
        A = seroparams['A_emcee'][j]
        B = seroparams['B_emcee'][j]
        # if person not carrying
        if carriers.iloc[s,j] == 0:
            # if serotype is NVT or if person is not vaccinated
            ind_lambdas.iloc[s,j] = totFOIbysero[carriers['agegroup'][s],j]
            #ind_lambdas[s,j] = totFOIbysero[carriers['agegroup'][s],j] # if np array
        # if person is carrying, then recovery rate as function of prevcol
        else:
            #ind_lambdas[s,j] = (1/np.random.gamma(shape = alpha, scale = (A*np.exp(-prevcol*B))/alpha))  # np array
            ind_lambdas.iloc[s,j] = (1/np.random.gamma(shape = alpha, scale = (A*np.exp(-prevcol*B))/alpha))

#### IMMUNITY AND COMPETITION PARAMS taken from Cobey & Lipsitch
#sigma = 0.4 # specific immunity
#epsilon = 0.1 # non-specific immunity
#theta = 0.01 # competition
#mincarrdur = 5 # minimum carriage duration

# vaccine params
vacc_start = 7300 # start of vaccination period
#vacc_dur = 10 # years of vaccine-induced protection
vacc_dur = vacc_dur*365 # protection converted to days 
#vacc_eff = 1

# log file of main params
paramlogfile = dict(simnum = max(filenums) + 1, country = country, beta = seroparams['beta_child'][0], indivs = indivs, vacc_dur = vacc_dur, vacc_eff = vacc_eff,
                    sigma = sigma, epsilon = epsilon, theta = theta)
paramdf = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in paramlogfile.items() ])) 
filedestparams = filedest + endstring + 'paramlogfile.csv'
paramdf.to_csv(path_or_buf = filedestparams, index = False)

t = 0
stop_t = vacc_start # 20 years to reach SS without vacc
# initial tot_lambda
#tot_lambda1 = ind_lambdas.sum().sum()#(sum(carriers.iloc[:,0])*seroparams['beta'][0])
    
# outputs for graph of time vs infected ppl
chunkSize = 100000
currentLength = chunkSize
DS_inf_nams = np.unique(seroparams['serotype'])#list(map(lambda x: 'infxtd' + str(x+1), range(num_sero)))
plotsy = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = DS_inf_nams)
plotsy.insert(loc = 0, column = 'time', value = np.nan)
    
next_ev = 1
start = time.time()

for t in range(stop_t+1):
#while t <= stop_t and tot_lambda1 != 0:
    
    if len(np.where(carriers.age > age_breakdown[2])[0]) > 0:
        ind_to_die = np.where(carriers.age > age_breakdown[2])[0]
        carriers.loc[ind_to_die] = 0 # making this person a baby now
        carriers.loc[ind_to_die, 'susc'] = max_carr_cap # susc = max carrying capacity when nothing is carried
        carriers.loc[ind_to_die, 'vacc'] = False
        prev_col.loc[ind_to_die] = 0
        
    index_children = np.where(carriers.age < age_breakdown[0])[0]
    index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
    index_elderly = np.where(carriers.age >= age_breakdown[1])[0]
    
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
    
    # tot.lambda = total rate of all/any events
    numindivspertime = round(indivs/50)
    i = 0
    U0 = np.random.uniform(0,1, size = numindivspertime)
    U3 = np.random.uniform(0,1, size = numindivspertime)
    
    #while i < numindivspertime: # update 20 individuals each time step for indivs = 1000
    for i in range(numindivspertime):
        tot_lambda1 = ind_lambdas.sum().sum()
        #print(tot_lambda1)
        U1 = U0[i]
        U1_tot_lambda1 = U1*tot_lambda1
        
        # choose person index ( = which person will this happen to?)
        cumsum_person = ind_lambdas.sum(axis = 1)
        which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        #try: # DEBUGGING
        #    which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        #except IndexError:
        #    which_person = np.where(np.random.uniform(0,1)*ind_lambdas.sum().sum() < np.cumsum(cumsum_person))[0][0]#which_person = indivs - 1 #
        #    break
        #print(which_person)
        
        # choose event that will take place (= which serotype will be toggled?) #ind_lambdas[which_person,:], 
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
        
        # switching event (S --> I or I --> S)
        if which_indices[0].size > 0:
            which_sero = which_indices[0][0]
            carriers.iloc[which_person, which_sero] = ~carriers.iloc[which_person, which_sero] + 2    

        #print(which_sero)
        # reupdate the susc score
        new_susc = max_carr_cap - sum(carriers.iloc[which_person, 0:num_sero])
        
        # if person acquired new serotype, update prevcol unless person is elderly and updatet recov rate
        if (new_susc < carriers.susc.loc[which_person]):
            if (carriers.age.loc[which_person] < age_breakdown[1]): 
                prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1
            alphatemp = seroparams['alpha_emcee'][which_sero]
            Atemp = seroparams['A_emcee'][which_sero]
            Btemp = seroparams['B_emcee'][which_sero]
            ind_lambdas.iloc[which_person, which_sero] = (1/np.random.gamma(shape = alphatemp,
                                                                            scale = (Atemp*np.exp(-prev_col.iloc[which_person, which_sero]*Btemp))/alphatemp))
            #ind_lambdas[which_person, which_sero] = (1/np.random.gamma(shape = alphatemp,
        
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        
    # migration: add 1 carrier to serotypes with zero carriers
    extinctsero = np.where(carriers.iloc[:,0:num_sero].sum(axis = 0) == 0)[0]
    
    #extinctsero_child = np.where(carriers.iloc[index_children,0:num_sero].sum(axis = 0) == 0)[0]
    totchildcarr = carriers.iloc[index_children,0:num_sero].sum(axis = 0)
    #totchildcarr[extinctsero_child] = 1
    totchildcarr[extinctsero] = 1
    
    #extinctsero_adult = np.where(carriers.iloc[index_adults,0:num_sero].sum(axis = 0) == 0)[0]
    totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
    #totadultcarr[extinctsero_adult] = 1
    
    #extinctsero_elder = np.where(carriers.iloc[index_elderly,0:num_sero].sum(axis = 0) == 0)[0]
    toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)
    #toteldercarr[extinctsero_elder] = 1
    
    totFOIbysero_children = [0]*num_sero
    totFOIbysero_adults = [0]*num_sero
    totFOIbysero_elders = [0]*num_sero
    
    for i in range(len(seroparams)):
        totFOIbysero_children[i] = ((contact_matrix[0][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                    (contact_matrix[0][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+
                                    (contact_matrix[0][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
        totFOIbysero_adults[i] = ((contact_matrix[1][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                  (contact_matrix[1][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                                  (contact_matrix[1][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
        totFOIbysero_elders[i] =  ((contact_matrix[2][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                   (contact_matrix[2][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                                   (contact_matrix[2][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
   
    totFOIbysero = np.row_stack((totFOIbysero_children, totFOIbysero_adults, totFOIbysero_elders))
    
    if np.any(totFOIbysero < 0):
            print('NEGATIVE FOI')
            print("sero:", seroparams.serotype[np.where(totFOIbysero < 0)[1]])
            print("indiv:", np.where(carriers_all[seroparams.serotype[np.where(totFOIbysero < 0)[1]]] > 0))
    
    # update FOI            
    for s in range(indivs):
        for j in range(num_sero):
            # if person not carrying
            if carriers.iloc[s,j] == 0:
                #ind_lambdas[s,j] = totFOIbysero[carriers['agegroup'][s],j]
                ind_lambdas.iloc[s,j] = totFOIbysero[carriers['agegroup'][s],j]


    # update time and age
    #t += 1
    next_ev += 1
    carriers.age += 1
        
end = time.time()
total_time = end - start
plotsy = plotsy.dropna()
plotsy_melt= pd.melt(plotsy, id_vars = 'time')
plot = plotsy_melt.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

plotsy_meltdest = filedest + endstring + 'carrdt.csv'
plotsy_melt.to_csv(path_or_buf= plotsy_meltdest, index = True)

# output file
sim_output = plotsy.loc[len(plotsy) - 1]
sim_output['sim_time'] = total_time/3600

newfiledest = filedest + endstring +  "output" +".csv"
sim_output.to_csv(path_or_buf= newfiledest, index = True)

# prevcol
prevcolnam = filedest +  endstring + "prevcol" +".csv"
prev_col.to_csv(path_or_buf = prevcolnam, index = True)

# carriers
carrnam =  filedest + endstring + "carrierstmax" +".csv"
carriers.to_csv(path_or_buf = carrnam, index = True)

# output figure
plotname = filedest + endstring + "simulation" + ".pdf"
plot.get_figure().savefig(plotname, bbox_inches = "tight")

sim_output = pd.DataFrame({'key':sim_output.index, 'value':sim_output.values})
seroprev = sim_output.drop([0, num_sero+1])
seroprev = seroprev.rename(columns = {'key': 'Serotype', 'value': 'carrprev'})
seroprev['carrprev'] = seroprev['carrprev']/indivs

### DISEASE INCIDENCE USING INVASIVENESS
disinc = seroprev.merge(country_inv, on = 'Serotype')
disease_cas = disinc['carrprev']*disinc['invasiveness']*indivs*(t/365)
disinc['disease_cas'] = disease_cas
disinc['disease_inc'] = (disease_cas/indivs)*100000

# Disease Incidence
index_children = np.where(carriers.age < age_breakdown[0])[0]
index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
index_elderly = np.where(carriers.age >= age_breakdown[1])[0]

pop_child = len(index_children)
pop_adult = len(index_adults)
pop_elder = len(index_elderly)

# Carriage by serotype
carrprevchildsero = np.array(carriers[carriers.agegroup == 0].iloc[:,0:num_sero].sum(axis = 0))/pop_child
carrprevadultsero = np.array(carriers[carriers.agegroup == 1].iloc[:,0:num_sero].sum(axis = 0))/pop_adult
carrpreveldersero = np.array(carriers[carriers.agegroup == 2].iloc[:,0:num_sero].sum(axis = 0))/pop_elder

disinc['carrprevchild'] = carrprevchildsero
disinc['carrprevadult'] = carrprevadultsero
disinc['carrprevelder'] = carrpreveldersero

# disease cases by serotype
disease_child = disinc['carrprevchild']*disinc['invasiveness']*(t/365)*pop_child
disease_adult = disinc['carrprevadult']*disinc['invasiveness']*(t/365)*pop_adult
disease_elder = disinc['carrprevelder']*disinc['invasiveness']*(t/365)*pop_elder

# save to file destination:
newfiledest = filedest + endstring + "disinc" + ".csv"
disinc.to_csv(path_or_buf= newfiledest, index = True)

# overall IPD incidence by age group
overalldis = (disinc['disease_cas'].sum())/(indivs*(t/365))*100000
overalldischild = (sum(disease_child)/(pop_child*(t/365)))*100000
overalldisadult = (sum(disease_adult)/(pop_adult*(t/365)))*100000
overalldiselder = (sum(disease_elder)/(pop_elder*(t/365)))*100000

# overall carr prev by age group
# overallcarrprevchild = carriers[agegroup == 0].iloc[:,0:num_sero].sum().sum()/pop_child #.605
# overallcarrprevadult = carriers[agegroup == 1].iloc[:,0:num_sero].sum().sum()/pop_adult #.607
# overallcarrprevelder = carriers[agegroup == 2].iloc[:,0:num_sero].sum().sum()/pop_elder #.5

overallcarrprev = len(np.where(carriers.susc < max_carr_cap)[0])/indivs
overallcarrprevchild = len(np.where(carriers[carriers.agegroup == 0].susc < max_carr_cap)[0])/pop_child
overallcarrprevadult = len(np.where(carriers[carriers.agegroup == 1].susc < max_carr_cap)[0])/pop_adult
overallcarrprevelder = len(np.where(carriers[carriers.agegroup == 2].susc < max_carr_cap)[0])/pop_elder

agegrpstats = np.array([[overallcarrprev, overallcarrprevchild, overallcarrprevadult, overallcarrprevelder],
                        [overalldis, overalldischild, overalldisadult, overalldiselder]])

agegrpstats = pd.DataFrame(agegrpstats, columns = ['overall','child', 'adult', 'elder'], index = ['carrprev', 'disease inc'])

newfiledest = filedest + endstring + "agegrpstats" + ".csv"
agegrpstats.to_csv(path_or_buf= newfiledest, index = True)

#################################################################################################################################################################
############# PART 2: POST-VACCINATION PERIOD 1 #################################################################################################################
#################################################################################################################################################################


def new_lambdas2(carrier_person, ind_lambdas, prev_col, sigma):    
    susc = carrier_person.susc
    vacc = carrier_person.vacc
    
    index = np.where(carrier_person[0:num_sero] == 0) # index of serotypes not carried
    non_index = np.bool_(carrier_person[0:num_sero]) # index of serotypes carried
    vt_index = index[0][np.isin(carrier_person.index[index], vt)] # index of VT serotypes not carried 
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
        ind_lambdas[vt_index] = 1/vacc_dur
    
    # serotypes you are carrying change recovery rate based on colonization history and immunity
    ind_lambdas[non_index] = 1/(mincarrdur+(((1/ind_lambdas[non_index])-mincarrdur)*math.exp(-eps*sum(prev_col))))
    return(ind_lambdas)



### Top 13 disease-causing serotypes get vaccinated against
vt = np.array(disinc.nlargest(13, 'disease_inc').Serotype)
vt_ind = np.where(np.isin(colnames, vt))[0]

vt0 = vt # for dataframe concatenation later

# top13 = np.array(disinc.disease_inc)[np.argsort(disinc.disease_inc)][-13:]
# top_ind = np.isin(disinc.disease_inc, top13)
# vt = disinc.loc[top_ind].Serotype
# vt_ind = np.where(np.isin(colnames, vt))[0]

vacc_stop_t = t + 3650 #1825 # current time + 5 years

# outputs for graph of time vs infected ppl
currentLength = chunkSize
plotsy2 = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = DS_inf_nams)
plotsy2.insert(loc = 0, column = 'time', value = np.nan)

next_ev = 1
start2 = time.time()

#while t <= vacc_stop_t and tot_lambda1 != 0:
for t in range(t+1, vacc_stop_t + 1):
    
    if len(np.where(carriers.age > age_breakdown[2])[0]) > 0:
        ind_to_die = np.where(carriers.age > age_breakdown[2])[0]
        carriers.loc[ind_to_die] = 0 # making this person a baby now
        carriers.loc[ind_to_die, 'susc'] = max_carr_cap # susc = max carrying capacity when nothing is carried
        carriers.loc[ind_to_die, 'vacc'] = False
        prev_col.loc[ind_to_die] = 0
        
    index_children = np.where(carriers.age < age_breakdown[0])[0]
    index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
    index_elderly = np.where(carriers.age >= age_breakdown[1])[0]
    
    carriers.loc[index_children, 'agegroup'] = 0
    carriers.loc[index_adults, 'agegroup'] = 1
    carriers.loc[index_elderly, 'agegroup'] = 2
        
    # add NAs to end of output dataframe if there is no space left to fill
    if next_ev > currentLength:
        plotsy2 = plotsy2.append(pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = plotsy2.columns), ignore_index = True)
        currentLength = currentLength + chunkSize
    
    # only record data when t > 9000 (i.e. last 3 years of time)
    plotsy2.iloc[next_ev-1, 1:] = np.array(carriers.iloc[:,0:num_sero].sum(axis = 0))
    plotsy2.time[next_ev-1] = t
    
    #if t >= 7300:
    index_vacc = np.where((carriers.age >= 730) & (carriers.age <= 1095) & (carriers.vacc == False)) # index of ppl to be vaccinated
    num_newvacc = round(vacc_eff*len(index_vacc[0])) # number of new people to vaccinate (= vacc efficacy * tot number of people)
    adjusted_ind = index_vacc[0][0:num_newvacc] # indices of people for whom vaccination works/is protective
    carriers.loc[adjusted_ind, 'vacc'] = True
        
        
    # tot.lambda = total rate of all/any events
    i = 0
    U0 = np.random.uniform(0,1, size = numindivspertime)#20)
    U3 = np.random.uniform(0,1, size = numindivspertime)
    #while i < numindivspertime:#20: # update 50 individuals each time step
    for i in range(numindivspertime):
        tot_lambda1 = ind_lambdas.sum().sum()
        #print(tot_lambda1)
        U1 = U0[i]
        U1_tot_lambda1 = U1*tot_lambda1
        
        # choose person index ( = which person will this happen to?)
        cumsum_person = ind_lambdas.sum(axis = 1)
        which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        
        # choose event that will take place (= which serotype will be toggled?) #ind_lambdas[which_person,],
        person_lambdas = new_lambdas2(carriers.loc[which_person], ind_lambdas.loc[which_person], 
                                      prev_col.loc[which_person], sigma)
            
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
        
        # switching event (S --> I or I --> S)
        if which_indices[0].size > 0:
            which_sero = which_indices[0][0]
            carriers.iloc[which_person, which_sero] = ~carriers.iloc[which_person, which_sero] + 2        

        #print(which_sero)
        # reupdate the susc score
        new_susc = max_carr_cap - sum(carriers.iloc[which_person, 0:num_sero])
    
        # update colonisation history and recovery rate if person acquired new serotype    
        if new_susc < carriers.susc.loc[which_person]: 
            if (carriers.age.loc[which_person] < age_breakdown[1]):
                prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1
            alphatemp = seroparams['alpha_emcee'][which_sero]
            Atemp = seroparams['A_emcee'][which_sero]
            Btemp = seroparams['B_emcee'][which_sero]
            #ind_lambdas[which_person, which_sero] = (1/np.random.gamma(shape = alphatemp, # np array                                                 
            ind_lambdas.iloc[which_person, which_sero] = (1/np.random.gamma(shape = alphatemp, 
                                                                            scale = (Atemp*np.exp(-prev_col.iloc[which_person, which_sero]*Btemp))/alphatemp))  
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        #i += 1
    
    #migration: add 1 carrier to serotypes with zero carriers
    extinctsero = np.where(carriers.iloc[:,0:num_sero].sum(axis = 0) == 0)[0]
    totchildcarr = carriers.iloc[index_children,0:num_sero].sum(axis = 0)
    totchildcarr[extinctsero] = 1
    
    #extinctsero_adult = np.where(carriers.iloc[index_adults,0:num_sero].sum(axis = 0) == 0)[0]
    totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
    #totadultcarr[extinctsero_adult] = 1
    
    #extinctsero_elder = np.where(carriers.iloc[index_elderly,0:num_sero].sum(axis = 0) == 0)[0]
    toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)
    #toteldercarr[extinctsero_elder] = 1
    
    totFOIbysero_children = [0]*num_sero
    totFOIbysero_adults = [0]*num_sero
    totFOIbysero_elders = [0]*num_sero
    
    for i in range(len(seroparams)):
        totFOIbysero_children[i] = ((contact_matrix[0][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                    (contact_matrix[0][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+
                                    (contact_matrix[0][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
        totFOIbysero_adults[i] = ((contact_matrix[1][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                  (contact_matrix[1][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                                  (contact_matrix[1][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
        totFOIbysero_elders[i] =  ((contact_matrix[2][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                   (contact_matrix[2][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                                   (contact_matrix[2][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
   
    totFOIbysero = np.row_stack((totFOIbysero_children, totFOIbysero_adults, totFOIbysero_elders))
    
    if np.any(totFOIbysero < 0):
        print('NEGATIVE FOI')
        print("sero:", seroparams.serotype[np.where(totFOIbysero < 0)[1]])
        print("indiv:", np.where(carriers_all[seroparams.serotype[np.where(totFOIbysero < 0)[1]]] > 0))
    
    # individual rates                
    for s in range(indivs):
        for j in range(num_sero):
            # if person not carrying
            if carriers.iloc[s,j] == 0:
                # if person is vaccinated & serotype is VT
                if carriers.vacc[s] == True and np.isin(j,vt_ind):
                    #ind_lambdas[s,j] = 1/vacc_dur
                    ind_lambdas.iloc[s,j] = 1/vacc_dur
                # if serotype is NVT or if person is not vaccinated
                else:
                    #ind_lambdas[s,j] = totFOIbysero[carriers['agegroup'][s],j] 
                    ind_lambdas.iloc[s,j] = totFOIbysero[carriers['agegroup'][s],j]             

    # update time and age
    #t += 1
    next_ev += 1
    carriers.age += 1

end2 = time.time()
total_time2 = end2 - start2
plotsy2 = plotsy2.dropna()
plotsy_melt2= pd.melt(plotsy2, id_vars = 'time')
plot2 = plotsy_melt2.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

# output file
sim_output2 = plotsy2.loc[len(plotsy2) - 1]
sim_output2['sim_time'] = total_time2/3600

phase = 'postPCV13'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

plotsy_melt2dest = filedest + endstring + 'carrdt.csv'
plotsy_melt2.to_csv(path_or_buf= plotsy_melt2dest, index = True)

newfiledest2 = filedest + endstring + "output" + ".csv"
sim_output2.to_csv(path_or_buf= newfiledest2, index = True)

# prevcol
prevcolnam2 = filedest + endstring + "prevcol" + ".csv"
prev_col.to_csv(path_or_buf = prevcolnam2, index = True)

# carriers
carrnam2 =  filedest + endstring + "carrierstmax" + ".csv"
carriers.to_csv(path_or_buf = carrnam2, index = True)

# output figure
plotname2 = filedest + endstring + "simulation" + ".pdf"
plot2.get_figure().savefig(plotname2, bbox_inches = "tight")

# Age group stats
index_children = np.where(carriers.age < age_breakdown[0])[0]
index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
index_elderly = np.where(carriers.age >= age_breakdown[1])[0]

pop_child = len(index_children)
pop_adult = len(index_adults)
pop_elder = len(index_elderly)

sim_output2 = pd.DataFrame({'key':sim_output2.index, 'value':sim_output2.values})
seroprev2 = sim_output2.drop([0, num_sero+1])
seroprev2 = seroprev2.rename(columns = {'key': 'Serotype', 'value': 'carrprev'})
seroprev2['carrprev'] = seroprev2['carrprev']/indivs # overall carr prev by serotype (regardless of age group)

### Get disease incidence using invasiveness
disinc2 = seroprev2.merge(country_inv, on = 'Serotype')
tprime = vacc_stop_t - stop_t
disease_cas2 = disinc2['carrprev']*disinc2['invasiveness']*indivs*(tprime/365)
disinc2['disease_cas'] = disease_cas2
disinc2['disease_inc'] = (disease_cas2/indivs)*100000 # overall disease incidence by serotype (regardless of age group)

# carriage prevalence by serotype and age grp
carrprevchildsero2 = np.array(carriers[carriers.agegroup == 0].iloc[:,0:num_sero].sum(axis = 0))/pop_child
carrprevadultsero2 = np.array(carriers[carriers.agegroup == 1].iloc[:,0:num_sero].sum(axis = 0))/pop_adult
carrpreveldersero2 = np.array(carriers[carriers.agegroup == 2].iloc[:,0:num_sero].sum(axis = 0))/pop_elder

disinc2['carrprevchild'] = carrprevchildsero2
disinc2['carrprevadult'] = carrprevadultsero2
disinc2['carrprevelder'] = carrpreveldersero2

# disease cases per serotype and age grp
disease_child2 = disinc2['carrprevchild']*disinc2['invasiveness']*(tprime/365)*pop_child
disease_adult2 = disinc2['carrprevadult']*disinc2['invasiveness']*(tprime/365)*pop_adult
disease_elder2 = disinc2['carrprevelder']*disinc2['invasiveness']*(tprime/365)*pop_elder

# save to file destination:
newfiledest2 = filedest + endstring + "disinc" + ".csv"
disinc2.to_csv(path_or_buf= newfiledest2, index = True)

# overall age group stats (all serotypes): carriage prevalence and IPD incidence in each age group
overalldis2 = (disinc2['disease_cas'].sum())/(indivs*(tprime/365))*100000
overalldischild2 = (sum(disease_child2)/(pop_child*(tprime/365)))*100000
overalldisadult2 = (sum(disease_adult2)/(pop_adult*(tprime/365)))*100000
overalldiselder2 = (sum(disease_elder2)/(pop_elder*(tprime/365)))*100000

overallcarrprev2 = len(np.where(carriers.susc < max_carr_cap)[0])/indivs
overallcarrprevchild2 = len(np.where(carriers[carriers.agegroup == 0].susc < max_carr_cap)[0])/pop_child
overallcarrprevadult2 = len(np.where(carriers[carriers.agegroup == 1].susc < max_carr_cap)[0])/pop_adult
overallcarrprevelder2 = len(np.where(carriers[carriers.agegroup == 2].susc < max_carr_cap)[0])/pop_elder

agegrpstats2 = np.array([[overallcarrprev2, overallcarrprevchild2, overallcarrprevadult2, overallcarrprevelder2],
                        [overalldis2, overalldischild2, overalldisadult2, overalldiselder2]])

agegrpstats2 = pd.DataFrame(agegrpstats2, columns = ['overall','child', 'adult', 'elder'], index = ['carrprev', 'disease inc'])

newfiledest = filedest + endstring + "agegrpstats" + ".csv"
agegrpstats2.to_csv(path_or_buf= newfiledest, index = True)

#################################################################################################################################################################
############# PART 3: POST-VACCINATION PERIOD 2 - extension of PCV ##############################################################################################
#################################################################################################################################################################

### Top 20 disease-causing serotypes get vaccinated against
new_vt = np.array(disinc2.nlargest(10, 'disease_inc').Serotype)
new_vt_ind = np.where(np.isin(colnames, new_vt))[0]

# top7add = np.array(disinc2.disease_inc)[np.argsort(disinc2.disease_inc)][-20:-13]
# new_vt = np.array(disinc2.Serotype)[np.argsort(disinc2.disease_inc)][-20:-13]
# new_vt_ind = np.where(np.isin(colnames, new_vt))[0]
    
vt = np.concatenate((vt, new_vt))
vt_ind = np.concatenate((vt_ind, new_vt_ind))

vacc2_stop_t = t + 3650 # current time + 10 years

# outputs for graph of time vs infected ppl
currentLength = chunkSize
plotsy3 = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = DS_inf_nams)
plotsy3.insert(loc = 0, column = 'time', value = np.nan)

next_ev = 1
start3 = time.time()
#while t <= vacc2_stop_t and tot_lambda1 != 0:
for t in range(t+1, vacc2_stop_t+1):
    
    if len(np.where(carriers.age > age_breakdown[2])[0]) > 0:
        ind_to_die = np.where(carriers.age > age_breakdown[2])[0]
        carriers.loc[ind_to_die] = 0 # making this person a baby now
        carriers.loc[ind_to_die, 'susc'] = max_carr_cap # susc = max carrying capacity when nothing is carried
        carriers.loc[ind_to_die, 'vacc'] = False
        prev_col.loc[ind_to_die] = 0
        
    index_children = np.where(carriers.age < age_breakdown[0])[0]
    index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
    index_elderly = np.where(carriers.age >= age_breakdown[1])[0]
    
    carriers.loc[index_children, 'agegroup'] = 0
    carriers.loc[index_adults, 'agegroup'] = 1
    carriers.loc[index_elderly, 'agegroup'] = 2
        
    # add NAs to end of output dataframe if there is no space left to fill
    if next_ev > currentLength:
        plotsy3 = plotsy3.append(pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = plotsy3.columns), ignore_index = True)
        currentLength = currentLength + chunkSize
    
    # only record data when t > 9000 (i.e. last 3 years of time)
    plotsy3.iloc[next_ev-1, 1:] = np.array(carriers.iloc[:,0:num_sero].sum(axis = 0))
    plotsy3.time[next_ev-1] = t
    
    #if t >= 7300:
    index_vacc = np.where((carriers.age >= 730) & (carriers.age <= 1095) & (carriers.vacc == False)) # index of ppl to be vaccinated
    num_newvacc = round(vacc_eff*len(index_vacc[0])) # number of new people to vaccinate (= vacc efficacy * tot number of people)
    adjusted_ind = index_vacc[0][0:num_newvacc] # indices of people for whom vaccination works/is protective
    carriers.loc[adjusted_ind, 'vacc'] = True
        
               
    
    # tot.lambda = total rate of all/any events
    #i = 0
    U0 = np.random.uniform(0,1, size = numindivspertime) #20)
    U3 = np.random.uniform(0,1, size = numindivspertime) #20)
    #while i < numindivspertime: #20: # update 20 individuals each time step
    for i in range(numindivspertime):
        tot_lambda1 = ind_lambdas.sum().sum()
        #print(tot_lambda1)
        U1 = U0[i]
        U1_tot_lambda1 = U1*tot_lambda1
        
        # choose person index ( = which person will this happen to?)
        cumsum_person = ind_lambdas.sum(axis = 1)
        which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        
        # choose event that will take place (= which serotype will be toggled?)
        person_lambdas = new_lambdas2(carriers.loc[which_person], ind_lambdas.loc[which_person], prev_col.loc[which_person], sigma)
            
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
        
        # switching event (S --> I or I --> S)
        if which_indices[0].size > 0:
            which_sero = which_indices[0][0]
            carriers.iloc[which_person, which_sero] = ~carriers.iloc[which_person, which_sero] + 2        

        #print(which_sero)
        # reupdate the susc score
        new_susc = max_carr_cap - sum(carriers.iloc[which_person, 0:num_sero])
    
        # update colonisation history if person acquired new serotype    
        if new_susc < carriers.susc.loc[which_person]:
            if (carriers.age.loc[which_person] <= 23360): # if children or adult, update prev col
                prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1
            # update recovery rate
            alphatemp = seroparams['alpha_emcee'][which_sero]
            Atemp = seroparams['A_emcee'][which_sero]
            Btemp = seroparams['B_emcee'][which_sero]
            ind_lambdas.iloc[which_person, which_sero] = (1/np.random.gamma(shape = alphatemp, 
                                                                            scale = (Atemp*np.exp(-prev_col.iloc[which_person, which_sero]*Btemp))/alphatemp))  
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        #i += 1
    
    #migration: add 1 carrier to serotypes with zero carriers
    extinctsero = np.where(carriers.iloc[:,0:num_sero].sum(axis = 0) == 0)[0]
    totchildcarr = carriers.iloc[index_children,0:num_sero].sum(axis = 0)
    totchildcarr[extinctsero] = 1
    
    #extinctsero_adult = np.where(carriers.iloc[index_adults,0:num_sero].sum(axis = 0) == 0)[0]
    totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
    #totadultcarr[extinctsero_adult] = 1
    
    #extinctsero_elder = np.where(carriers.iloc[index_elderly,0:num_sero].sum(axis = 0) == 0)[0]
    toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)
    #toteldercarr[extinctsero_elder] = 1
    
    totFOIbysero_children = [0]*num_sero
    totFOIbysero_adults = [0]*num_sero
    totFOIbysero_elders = [0]*num_sero
    
    for i in range(len(seroparams)):
        totFOIbysero_children[i] = ((contact_matrix[0][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                    (contact_matrix[0][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+
                                    (contact_matrix[0][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
        totFOIbysero_adults[i] = ((contact_matrix[1][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                  (contact_matrix[1][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                                  (contact_matrix[1][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
        totFOIbysero_elders[i] =  ((contact_matrix[2][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                   (contact_matrix[2][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                                   (contact_matrix[2][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
   
    totFOIbysero = np.row_stack((totFOIbysero_children, totFOIbysero_adults, totFOIbysero_elders))
    
    # individual rates                
    for s in range(indivs):
        for j in range(num_sero):
            # if person not carrying
            if carriers.iloc[s,j] == 0:
                # if person is vaccinated & serotype is VT
                if carriers.vacc[s] == True and np.isin(j,vt_ind):
                    ind_lambdas.iloc[s,j] = 1/vacc_dur
                # if serotype is NVT or if person is not vaccinated
                else:
                    ind_lambdas.iloc[s,j] = totFOIbysero[carriers['agegroup'][s],j]
                    
    # update time and age
    #t += 1
    next_ev += 1
    carriers.age += 1

end3 = time.time()
total_time3 = end3 - start3
plotsy3 = plotsy3.dropna()
plotsy_melt3= pd.melt(plotsy3, id_vars = 'time')
plot3 = plotsy_melt3.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

# output file
sim_output3 = plotsy3.loc[len(plotsy3) - 1]
sim_output3['sim_time'] = total_time3/3600

phase = 'postPCV20'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

plotsy_melt3dest = filedest + endstring + 'carrdt.csv'
plotsy_melt3.to_csv(path_or_buf= plotsy_melt3dest, index = True)

newfiledest3 = filedest + endstring + "output" + ".csv"
sim_output3.to_csv(path_or_buf= newfiledest3, index = True)

# prevcol
prevcolnam3 = filedest + endstring + "prevcol" + ".csv"
prev_col.to_csv(path_or_buf = prevcolnam3, index = True)

# carriers
carrnam3 =  filedest + endstring + "carrierstmax" + ".csv"
carriers.to_csv(path_or_buf = carrnam3, index = True)

# output figure
plotname3 = filedest + endstring + "simulation" + ".pdf"
plot3.get_figure().savefig(plotname3, bbox_inches = "tight")

# SERO CAUSING IN DISEASE
index_children = np.where(carriers.age < age_breakdown[0])[0]
index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
index_elderly = np.where(carriers.age >= age_breakdown[1])[0]

pop_child = len(index_children)
pop_adult = len(index_adults)
pop_elder = len(index_elderly)

sim_output3 = pd.DataFrame({'key':sim_output3.index, 'value':sim_output3.values})
seroprev3 = sim_output3.drop([0, num_sero+1])
seroprev3 = seroprev3.rename(columns = {'key': 'Serotype', 'value': 'carrprev'})
seroprev3['carrprev'] = seroprev3['carrprev']/indivs

### Get disease incidence using invasiveness
disinc3 = seroprev3.merge(country_inv, on = 'Serotype')
tprime = vacc2_stop_t - vacc_stop_t
disease_cas3 = disinc3['carrprev']*disinc3['invasiveness']*indivs*(tprime/365)
disinc3['disease_cas'] = disease_cas3
disinc3['disease_inc'] = (disease_cas3/indivs)*100000

# overall carriage prevalence and IPD by age group
carrprevchildsero3 = np.array(carriers[carriers.agegroup == 0].iloc[:,0:num_sero].sum(axis = 0))/pop_child
carrprevadultsero3 = np.array(carriers[carriers.agegroup == 1].iloc[:,0:num_sero].sum(axis = 0))/pop_adult
carrpreveldersero3 = np.array(carriers[carriers.agegroup == 2].iloc[:,0:num_sero].sum(axis = 0))/pop_elder

disinc3['carrprevchild'] = carrprevchildsero3
disinc3['carrprevadult'] = carrprevadultsero3
disinc3['carrprevelder'] = carrpreveldersero3

# number of disease cases of each serotype
disease_child3 = disinc3['carrprevchild']*disinc3['invasiveness']*(tprime/365)*pop_child
disease_adult3 = disinc3['carrprevadult']*disinc3['invasiveness']*(tprime/365)*pop_adult
disease_elder3 = disinc3['carrprevelder']*disinc3['invasiveness']*(tprime/365)*pop_elder

# save to file destination:
newfiledest = filedest + endstring + "disinc" + ".csv"
disinc3.to_csv(path_or_buf= newfiledest, index = True)

# overall age group stats (allserotypes): carr prev and IPD incidence for each age group
overalldis3 = (disinc3['disease_cas'].sum())/(indivs*(tprime/365))*100000
overalldischild3 = (sum(disease_child3)/(pop_child*(tprime/365)))*100000
overalldisadult3 = (sum(disease_adult3)/(pop_adult*(tprime/365)))*100000
overalldiselder3 = (sum(disease_elder3)/(pop_elder*(tprime/365)))*100000

overallcarrprev3 = len(np.where(carriers.susc < max_carr_cap)[0])/indivs
overallcarrprevchild3 = len(np.where(carriers[carriers.agegroup == 0].susc < max_carr_cap)[0])/pop_child
overallcarrprevadult3 = len(np.where(carriers[carriers.agegroup == 1].susc < max_carr_cap)[0])/pop_adult
overallcarrprevelder3 = len(np.where(carriers[carriers.agegroup == 2].susc < max_carr_cap)[0])/pop_elder

agegrpstats3 = np.array([[overallcarrprev3, overallcarrprevchild3, overallcarrprevadult3, overallcarrprevelder3],
                        [overalldis3, overalldischild3, overalldisadult3, overalldiselder3]])

agegrpstats3 = pd.DataFrame(agegrpstats3, columns = ['overall','child', 'adult', 'elder'], index = ['carrprev', 'disease inc'])

newfiledest3 = filedest + endstring + "agegrpstats" + ".csv"
agegrpstats3.to_csv(path_or_buf= newfiledest3, index = True)

#################################################################################################################################################################
############# PART 4: POST-VACCINATION PERIOD 3 - extension of PCV ##############################################################################################
#################################################################################################################################################################

### Top 30 disease-causing serotypes get vaccinated against
new_vt2 = np.array(disinc3.nlargest(10, 'disease_inc').Serotype)
new_vt_ind2 = np.where(np.isin(colnames, new_vt2))[0]

# top10add = np.array(disinc3.disease_inc)[np.argsort(disinc3.disease_inc)][30:-20]
# new_vt2 = np.array(disinc3.Serotype)[np.argsort(disinc3.disease_inc)][-30:-20]
# new_vt_ind2 = np.where(np.isin(colnames, new_vt))[0]
    
# print VT into dataframe
vt_csv =  filedest + endstring + "VTsero" + ".csv"
d = dict(vt_1 = vt0, vt_2 = new_vt, vt_3 = new_vt2)  
vtdf = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in d.items() ])) 
vtdf.to_csv(path_or_buf = vt_csv, index = False)

vt = np.concatenate((vt, new_vt2))
vt_ind = np.concatenate((vt_ind, new_vt_ind2))

vacc3_stop_t = t + 3650 # current time + 10 years

# outputs for graph of time vs infected ppl
currentLength = chunkSize
plotsy4 = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = DS_inf_nams)
plotsy4.insert(loc = 0, column = 'time', value = np.nan)

next_ev = 1
start4 = time.time()
#while t <= vacc2_stop_t and tot_lambda1 != 0:
for t in range(t+1, vacc3_stop_t+1):
    
    if len(np.where(carriers.age > age_breakdown[2])[0]) > 0:
        ind_to_die = np.where(carriers.age > age_breakdown[2])[0]
        carriers.loc[ind_to_die] = 0 # making this person a baby now
        carriers.loc[ind_to_die, 'susc'] = max_carr_cap # susc = max carrying capacity when nothing is carried
        carriers.loc[ind_to_die, 'vacc'] = False
        prev_col.loc[ind_to_die] = 0
        
    index_children = np.where(carriers.age < age_breakdown[0])[0]
    index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
    index_elderly = np.where(carriers.age >= age_breakdown[1])[0]
    
    carriers.loc[index_children, 'agegroup'] = 0
    carriers.loc[index_adults, 'agegroup'] = 1
    carriers.loc[index_elderly, 'agegroup'] = 2
        
    # add NAs to end of output dataframe if there is no space left to fill
    if next_ev > currentLength:
        plotsy4 = plotsy4.append(pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = plotsy4.columns), ignore_index = True)
        currentLength = currentLength + chunkSize
    
    # only record data when t > 9000 (i.e. last 3 years of time)
    plotsy4.iloc[next_ev-1, 1:] = np.array(carriers.iloc[:,0:num_sero].sum(axis = 0))
    plotsy4.time[next_ev-1] = t
    
    #if t >= 7300:
    index_vacc = np.where((carriers.age >= 730) & (carriers.age <= 1095) & (carriers.vacc == False)) # index of ppl to be vaccinated
    num_newvacc = round(vacc_eff*len(index_vacc[0])) # number of new people to vaccinate (= vacc efficacy * tot number of people)
    adjusted_ind = index_vacc[0][0:num_newvacc] # indices of people for whom vaccination works/is protective
    carriers.loc[adjusted_ind, 'vacc'] = True
        
    
    # tot.lambda = total rate of all/any events
    #i = 0
    U0 = np.random.uniform(0,1, size = numindivspertime) #20)
    U3 = np.random.uniform(0,1, size = numindivspertime) #20)
    #while i < numindivspertime: #20: # update 20 individuals each time step
    for i in range(numindivspertime):
        tot_lambda1 = ind_lambdas.sum().sum()
        #print(tot_lambda1)
        U1 = U0[i]
        U1_tot_lambda1 = U1*tot_lambda1
        
        # choose person index ( = which person will this happen to?)
        cumsum_person = ind_lambdas.sum(axis = 1)
        which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        
        # choose event that will take place (= which serotype will be toggled?)
        person_lambdas = new_lambdas2(carriers.loc[which_person], ind_lambdas.loc[which_person], prev_col.loc[which_person], sigma)
            
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
        
        # switching event (S --> I or I --> S)
        if which_indices[0].size > 0:
            which_sero = which_indices[0][0]
            carriers.iloc[which_person, which_sero] = ~carriers.iloc[which_person, which_sero] + 2        

        #print(which_sero)
        # reupdate the susc score
        new_susc = max_carr_cap - sum(carriers.iloc[which_person, 0:num_sero])
    
        # update colonisation history if person acquired new serotype    
        if new_susc < carriers.susc.loc[which_person]:
            if (carriers.age.loc[which_person] <= 23360): # if children or adult, update prev col
                prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1
            # update recovery rate
            alphatemp = seroparams['alpha_emcee'][which_sero]
            Atemp = seroparams['A_emcee'][which_sero]
            Btemp = seroparams['B_emcee'][which_sero]
            ind_lambdas.iloc[which_person, which_sero] = (1/np.random.gamma(shape = alphatemp, 
                                                                            scale = (Atemp*np.exp(-prev_col.iloc[which_person, which_sero]*Btemp))/alphatemp))  
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        #i += 1
    
    #migration: add 1 carrier to serotypes with zero carriers
    extinctsero = np.where(carriers.iloc[:,0:num_sero].sum(axis = 0) == 0)[0]
    totchildcarr = carriers.iloc[index_children,0:num_sero].sum(axis = 0)
    totchildcarr[extinctsero] = 1
    
    #extinctsero_adult = np.where(carriers.iloc[index_adults,0:num_sero].sum(axis = 0) == 0)[0]
    totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
    #totadultcarr[extinctsero_adult] = 1
    
    #extinctsero_elder = np.where(carriers.iloc[index_elderly,0:num_sero].sum(axis = 0) == 0)[0]
    toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)
    #toteldercarr[extinctsero_elder] = 1
    
    totFOIbysero_children = [0]*num_sero
    totFOIbysero_adults = [0]*num_sero
    totFOIbysero_elders = [0]*num_sero
    
    for i in range(len(seroparams)):
        totFOIbysero_children[i] = ((contact_matrix[0][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                    (contact_matrix[0][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+
                                    (contact_matrix[0][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
        totFOIbysero_adults[i] = ((contact_matrix[1][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                  (contact_matrix[1][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                                  (contact_matrix[1][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
        totFOIbysero_elders[i] =  ((contact_matrix[2][0]*seroparams['beta_child'][i]*totchildcarr[i]/pop_child)+ 
                                   (contact_matrix[2][1]*seroparams['beta_adult'][i]*totadultcarr[i]/pop_adult)+ 
                                   (contact_matrix[2][2]*seroparams['beta_elder'][i]*toteldercarr[i]/pop_elder))
   
    totFOIbysero = np.row_stack((totFOIbysero_children, totFOIbysero_adults, totFOIbysero_elders))
    
    # individual rates                
    for s in range(indivs):
        for j in range(num_sero):
            # if person not carrying
            if carriers.iloc[s,j] == 0:
                # if person is vaccinated & serotype is VT
                if carriers.vacc[s] == True and np.isin(j,vt_ind):
                    ind_lambdas.iloc[s,j] = 1/vacc_dur
                # if serotype is NVT or if person is not vaccinated
                else:
                    ind_lambdas.iloc[s,j] = totFOIbysero[carriers['agegroup'][s],j]
                    
    # update time and age
    #t += 1
    next_ev += 1
    carriers.age += 1

end4 = time.time()
total_time4 = end4 - start4
plotsy4 = plotsy4.dropna()
plotsy_melt4= pd.melt(plotsy4, id_vars = 'time')
plot4 = plotsy_melt4.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

# output file
sim_output4 = plotsy4.loc[len(plotsy4) - 1]
sim_output4['sim_time'] = total_time4/3600

phase = 'postPCV30'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

plotsy_melt4dest = filedest + endstring + 'carrdt.csv'
plotsy_melt4.to_csv(path_or_buf= plotsy_melt4dest, index = True)

newfiledest4 = filedest + endstring + "output" + ".csv"
sim_output4.to_csv(path_or_buf= newfiledest4, index = True)

# prevcol
prevcolnam4 = filedest + endstring + "prevcol" + ".csv"
prev_col.to_csv(path_or_buf = prevcolnam4, index = True)

# carriers
carrnam4 =  filedest + endstring + "carrierstmax" + ".csv"
carriers.to_csv(path_or_buf = carrnam4, index = True)

# output figure
plotname4 = filedest + endstring + "simulation" + ".pdf"
plot4.get_figure().savefig(plotname4, bbox_inches = "tight")

# SERO CAUSING IN DISEASE
index_children = np.where(carriers.age < age_breakdown[0])[0]
index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
index_elderly = np.where(carriers.age >= age_breakdown[1])[0]

pop_child = len(index_children)
pop_adult = len(index_adults)
pop_elder = len(index_elderly)

sim_output4 = pd.DataFrame({'key':sim_output4.index, 'value':sim_output4.values})
seroprev4 = sim_output4.drop([0, num_sero+1])
seroprev4 = seroprev4.rename(columns = {'key': 'Serotype', 'value': 'carrprev'})
seroprev4['carrprev'] = seroprev4['carrprev']/indivs

### Get disease incidence using invasiveness
disinc4 = seroprev4.merge(country_inv, on = 'Serotype')
tprime = vacc3_stop_t - vacc2_stop_t
disease_cas4 = disinc4['carrprev']*disinc4['invasiveness']*indivs*(tprime/365)
disinc4['disease_cas'] = disease_cas4
disinc4['disease_inc'] = (disease_cas4/indivs)*100000

# overall carriage prevalence and IPD by age group
carrprevchildsero4 = np.array(carriers[carriers.agegroup == 0].iloc[:,0:num_sero].sum(axis = 0))/pop_child
carrprevadultsero4 = np.array(carriers[carriers.agegroup == 1].iloc[:,0:num_sero].sum(axis = 0))/pop_adult
carrpreveldersero4 = np.array(carriers[carriers.agegroup == 2].iloc[:,0:num_sero].sum(axis = 0))/pop_elder

disinc4['carrprevchild'] = carrprevchildsero4
disinc4['carrprevadult'] = carrprevadultsero4
disinc4['carrprevelder'] = carrpreveldersero4

# number of disease cases of each serotype
disease_child4 = round(disinc4['carrprevchild']*disinc4['invasiveness']*(tprime/365)*pop_child)
disease_adult4 = round(disinc4['carrprevadult']*disinc4['invasiveness']*(tprime/365)*pop_adult)
disease_elder4 = round(disinc4['carrprevelder']*disinc4['invasiveness']*(tprime/365)*pop_elder)

# save to file destination:
newfiledest = filedest + endstring + "disinc" + ".csv"
disinc4.to_csv(path_or_buf= newfiledest, index = True)

# overall age group stats (allserotypes): carr prev and IPD incidence for each age group
overalldis4 = (disinc4['disease_cas'].sum())/(indivs*(tprime/365))*100000
overalldischild4 = (sum(disease_child4)/(pop_child*(tprime/365)))*100000
overalldisadult4 = (sum(disease_adult4)/(pop_adult*(tprime/365)))*100000
overalldiselder4 = (sum(disease_elder4)/(pop_elder*(tprime/365)))*100000

overallcarrprev4 = len(np.where(carriers.susc < max_carr_cap)[0])/indivs
overallcarrprevchild4 = len(np.where(carriers[carriers.agegroup == 0].susc < max_carr_cap)[0])/pop_child
overallcarrprevadult4 = len(np.where(carriers[carriers.agegroup == 1].susc < max_carr_cap)[0])/pop_adult
overallcarrprevelder4 = len(np.where(carriers[carriers.agegroup == 2].susc < max_carr_cap)[0])/pop_elder

agegrpstats4 = np.array([[overallcarrprev4, overallcarrprevchild4, overallcarrprevadult4, overallcarrprevelder4],
                        [overalldis4, overalldischild4, overalldisadult4, overalldiselder4]])

agegrpstats4 = pd.DataFrame(agegrpstats4, columns = ['overall','child', 'adult', 'elder'], index = ['carrprev', 'disease inc'])

newfiledest4 = filedest + endstring + "agegrpstats" + ".csv"
agegrpstats4.to_csv(path_or_buf= newfiledest4, index = True)

###############################################################################################################################
######## Post-sim calc ########################################################################################################
###############################################################################################################################

# INCIDENCE RATE RATIO

# children
post13_child = overalldischild2/overalldischild
post20_child = overalldischild3/overalldischild2
post30_child = overalldischild4/overalldischild3
all_child = overalldischild4/overalldischild

IRR_child = np.array((post13_child, post20_child, post30_child, all_child))

# adults
post13_adult = overalldisadult2/overalldisadult
post20_adult = overalldisadult3/overalldisadult2
post30_adult = overalldisadult4/overalldisadult3
all_adult = overalldisadult4/overalldisadult

IRR_adult = np.array((post13_adult, post20_adult, post30_adult, all_adult))

# elderly
post13_elder = overalldiselder2/overalldiselder
post20_elder = overalldiselder3/overalldiselder2
post30_elder = overalldiselder4/overalldiselder3
all_elder = overalldiselder4/overalldiselder

IRR_elder = np.array((post13_elder, post20_elder, post30_elder, all_elder))

# overall
post13_all = overalldis2/overalldis
post20_all = overalldis3/overalldis2
post30_all = overalldis4/overalldis3
all_all = overalldis4/overalldis

IRR_all = np.array((post13_all, post20_all, post30_all, all_all))

IRR_df = pd.DataFrame({'children': IRR_child, 'adults': IRR_adult, 'elder': IRR_elder, 'overall': IRR_all}, 
                      index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])

newfiledest5 = filedest + endstring + "IRR_allperiods" + ".csv" 
IRR_df.to_csv(path_or_buf = newfiledest5, index = True) # incidence rate ratio dataframe

### INCIDENCE RATES

pre = agegrpstats.loc['disease inc']
post13 = agegrpstats2.loc['disease inc']
post20 = agegrpstats3.loc['disease inc']
post30 = agegrpstats4.loc['disease inc']

agespec = pd.concat([pre, post13, post20, post30], axis=1)
agespec.columns = ['pre', 'post13', 'post20', 'post30']

newfiledest5 = filedest + endstring + "disinc_eachperiod" + ".csv" 
agespec.to_csv(path_or_buf = newfiledest5, index = True) # incidence rate dataframe

### CARRIAGE PREVALENCE

pre_cp = agegrpstats.loc['carrprev']
post13_cp = agegrpstats2.loc['carrprev']
post20_cp = agegrpstats3.loc['carrprev']
post30_cp = agegrpstats4.loc['carrprev']

agespec_cp = pd.concat([pre_cp, post13_cp, post20_cp, post30_cp], axis=1)
agespec_cp.columns = ['pre', 'post13', 'post20', 'post30']

newfiledest5 = filedest + endstring + "carrprev_eachperiod" + ".csv" 
agespec_cp.to_csv(path_or_buf = newfiledest5, index = True) # carr prev dataframe

####### Testing with dummy carriage data ########################################################################################################################

# taken from this study https://academic.oup.com/jid/article/184/4/451/809030 
#serogroup = np.array(['6B', '23F', '19F', '6A', '11', '14', '15', '35', 'R', 'NT'])
#num_isolates = np.array([213, 268, 258, 140, 116, 67, 56, 49, 55, 17])
#list_of_tuples = list(zip(serogroup, num_isolates))
#carrprevdf = pd.DataFrame(list_of_tuples, columns = ['Serotype', 'Carr_isolates'])
#carrprevdf['Carr_prev'] = carrprevdf.Carr_isolates/1530 #1530 is the total number of swabs from the study

#################################################################################################################################################################

# import win32com.client as win32
# from win32 import win32api
# outlook = win32.Dispatch('outlook.application')
# mail = outlook.CreateItem(0)
# mail.To = 'a.lochen17@imperial.ac.uk'
# mail.Subject = country + 'simulation ' + str(simnum) + ' has finished running'
# mail.Body = 'Congrats!'
# mail.Send()