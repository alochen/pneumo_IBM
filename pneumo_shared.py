########## PNEUMO IBM with fitted sero params imported and with synchronous time where 20 events occur at each ############
########## time step, and a migration rate such that when a serotype goes extinct, it is reseeded into the     ############
########## population in the form of one infectious carrier                                                    ############

import pandas as pd
import numpy as np
# import math
import random
import time
# import glob
# import os
# from builtins import dict

# import parameters
modparams = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/modparams.csv")
country = modparams.country[0]
indivs = modparams.indivs[0]
sigma = modparams.sigma[0] # specific immunity
epsilon = modparams.epsilon[0] # non-specific immunity
theta = modparams.theta[0] # competition
mincarrdur = 5 # minimum carriage duration
vacc_dur = modparams.vacc_dur[0] # years of vaccine-induced protection
vacc_eff = modparams.vacc_eff[0] # vaccine efficacy
beta = modparams.beta[0] # transmission coefficient
max_carr_cap = 2 #1

seroparams = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/seroparams.csv")
seroparams = seroparams.drop(seroparams.index[seroparams.serotype.isin(['15B/C_old','15B','15C'])]).reset_index()
colnames = np.unique(seroparams['serotype']) #list(map(lambda x: 'serotype' + str(x+1), range(num_sero)))
num_sero = len(np.unique(seroparams['serotype']))


if country == "usa":
    country_inv = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/usalocalinv.csv")
if country == "fi":
    country_inv = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/finlandlocalinv.csv")
else:
    country_inv = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/fralocalinv.csv")

# socialmixr contact matrices that are symmetric; finland taken from POLYMOD, france from french survey, US from UK polymod
# age groups: 0-4 yrs, 5-64 yrs, 65+
contact_matrix_fi = np.array([[1.9670330,  7.003670, 0.3501001],
                              [0.4826382, 10.922424, 0.6010904],
                              [0.1190129,  2.965144, 1.8888889]])

contact_matrix_fr = np.array([[4.6754967, 12.81338, 0.720343],
                              [1.0412869, 23.92071, 2.223987],
                              [0.2638627, 10.02454, 5.290244]])

contact_matrix_usa = np.array([[1.9157895,  6.348054, 0.3648945],
                               [0.4649120, 11.308140, 0.9242947],
                               [0.1307877,  4.523555, 1.7142857]])


if country == 'usa':
    contact_matrix = contact_matrix_usa
    maxage = 78.5 # life expectancy USA 78.5 years (2010) https://datacatalog.worldbank.org/dataset/world-development-indicators 
    pop_breakdown = [0.06, 0.81, 0.13] # pop breakdown from UN
if country == 'fi':
    contact_matrix = contact_matrix_fi
    maxage = 79.9 # life expectancy 79.9 years (2010) https://datacatalog.worldbank.org/dataset/world-development-indicators 
    pop_breakdown = [0.05, 0.77, 0.18] # pop breakdown from UN
else:
    contact_matrix = contact_matrix_fr
    maxage = 81.7 # life expectancy 81.7 years (2010) https://datacatalog.worldbank.org/dataset/world-development-indicators 
    pop_breakdown = [0.06, 0.77, 0.17] # pop breakdown from UN


pop_child = int(pop_breakdown[0]*indivs)
pop_adult = int(pop_breakdown[1]*indivs)
pop_elder = int(pop_breakdown[2]*indivs)
agecutoffs = [5, 65, maxage]
age_breakdown = [round(i * 365) for i in agecutoffs] #[365*5, 23725, maxage*365]

# # weighted contacts like Cobey and Lipsitch??
# alpha_matrix = np.array([[contact_matrix[0][0]/pop_child, contact_matrix[0][1]/pop_child, contact_matrix[0][2]/pop_child],
#                          [contact_matrix[1][0]/pop_adult, contact_matrix[1][1]/pop_adult, contact_matrix[1][2]/pop_adult],
#                          [contact_matrix[2][0]/pop_elder, contact_matrix[2][1]/pop_elder, contact_matrix[2][2]/pop_elder]])

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

# beta (transmission or FOI) rates ## these will be the same for all serotypes
seroparams['beta_child'] = beta # 0.00010169#0.012 #0.021 #B_rates_child_val
seroparams['beta_adult'] = beta # 0.00010169#0.012 #0.021 #B_rates_adult_val
seroparams['beta_elder'] = beta # 0.00010169#0.012 #0.021 #B_rates_elder_val

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

# vaccine params
vacc_start = 7300 # start of vaccination period
vacc_dur = vacc_dur*365 # protection converted to days 

