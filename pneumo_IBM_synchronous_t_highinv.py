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
import os
#import matplotlib
#import matplotlib.pyplot as plt

#os.getcwd()
#os.chdir('Q:\\Technical\\Python')
country = "fi"
globalinv = False
phase = 'prevacc'
indivs = 1000

filedest = "//qdrive/homes/al817/Technical/Python/1000popruns/" + country + '_' + globalinv*("glob") + "inv/" 

seroparams = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/seroparams.csv")
seroparams = seroparams.drop(seroparams.index[seroparams.serotype.isin(['15B/C_old','15B','15C'])]).reset_index()
colnames = np.unique(seroparams['serotype']) #list(map(lambda x: 'serotype' + str(x+1), range(num_sero)))

mylist = [f for f in glob.glob(filedest + "*.csv")]
mylist = [fname.rsplit('.', 1)[0] for fname in os.listdir(filedest)]
filenums = [item.split('_',1)[0] for item in mylist]
filenums = [item.split('simulation',1)[-1] for item in filenums]
filenums = [int(i) for i in filenums] 
simnum = max(filenums) + 1
simnum = str(simnum) + '_' + globalinv*("glob") + "inv" 
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"
seroparams.to_csv(path_or_buf = filedest + endstring + 'dummy.csv', index = False)

#################################################################################################################################################################
############# Invasive serotypes ################################################################################################################################
#################################################################################################################################################################

inv = pd.read_csv(r"//qdrive/homes/al817/Technical/R/Case-to-carrier/abs.child-new.csv")

finland_inv = inv.loc[np.where(inv.DS == 'Finland.pre.PCV')[0]]
usa_inv = inv.loc[np.where(inv.DS == 'Massachusetts.post.PCV7')[0]]  
fra_inv = inv.loc[np.where(inv.DS == 'France.post.PCV7')[0]] 

if country == "usa":
    country_inv = usa_inv
if country == "fi":
    country_inv = finland_inv
else:
    country_inv = fra_inv

country_inv = country_inv.drop(['Unnamed: 0', 'DS', 'Serogroup', 'carriage', 'disease', 
                                'n.swab', 'N', 'time.int', 'lambda', 'lambda.low', 'lambda.high', 
                                'carr.prev', 'carr.prev.low', 'carr.prev.high', 'agegrp'], axis = 1)

global_inv = pd.read_csv(r"//qdrive/homes/al817/Technical/R/Case-to-carrier/consol.child-new.csv")
global_inv = global_inv.drop(['Unnamed: 0', 'Serogroup'], axis = 1)
global_inv.columns = country_inv.columns

sero_nocountryinv = np.setdiff1d(colnames,country_inv.Serotype)

sero_globalinvadd = sero_nocountryinv[pd.Series(sero_nocountryinv).isin(global_inv.Serotype)]
country_inv = country_inv.append(global_inv.loc[np.where(global_inv.Serotype.isin(sero_globalinvadd))]).reset_index().drop(['index'],axis = 1)

sero_nocountryinv2 = np.setdiff1d(colnames,country_inv.Serotype)

setinv = np.array([0.0001, 0, 0.0001]) # inv for serotypes without inv estimate
setinvdf = pd.DataFrame(np.tile(setinv, (len(sero_nocountryinv2),1)), columns = country_inv.columns[1:])
setinvdf['Serotype'] = sero_nocountryinv2

country_inv = country_inv.append(setinvdf).reset_index().drop(['index'],axis = 1)
country_inv.iloc[np.where(country_inv.Serotype == 'NT')[0], 1:4] = 0

serotypesnotincarrdat = np.setdiff1d(country_inv.Serotype, colnames)
country_inv = country_inv.drop(np.where(country_inv.Serotype.isin(serotypesnotincarrdat))[0])

# top 13 serotypes get vaccinated against
top13 = np.array(country_inv.invasiveness)[np.argsort(country_inv.invasiveness)][-13:]
vt = np.array(country_inv.Serotype)[np.argsort(country_inv.invasiveness)][-13:]
vt_ind = np.where(np.isin(colnames, vt))[0]

# top 13 invasive serotypes get vaccinated against
top7add = np.array(country_inv.invasiveness)[np.argsort(country_inv.invasiveness)][-20:-13]
new_vt = np.array(country_inv.Serotype)[np.argsort(country_inv.invasiveness)][-20:-13]
new_vt_ind = np.where(np.isin(colnames, new_vt))[0]

# top 13 globally invasive serotypes get vaccinated against
if globalinv == True:
    top13 = np.array(global_inv.invasiveness)[np.argsort(global_inv.invasiveness)][-13:]
    vt = np.array(global_inv.Serotype)[np.argsort(global_inv.invasiveness)][-13:]
    vt_ind = np.where(np.isin(colnames, vt))[0]
    top7add = np.array(global_inv.invasiveness)[np.argsort(global_inv.invasiveness)][-20:-13]
    new_vt = np.array(global_inv.Serotype)[np.argsort(global_inv.invasiveness)][-20:-13]
    new_vt_ind = np.where(np.isin(colnames, new_vt))[0]
    
# print VT into dataframe
vt_csv =  filedest + endstring + "VTsero" + ".csv"
d = dict(vt_1 = vt, vt_2 = new_vt)  
vtdf = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in d.items() ])) 
vtdf.to_csv(path_or_buf = vt_csv, index = False)

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

# import fitted carriage duration params from IBM_paramfitting.py
# Initialise population
num_sero = len(np.unique(seroparams['serotype']))
max_carr_cap = 1
pop_breakdown = [0.2, 0.6, 0.2] # has to sum to 1
pop_child = int(pop_breakdown[0]*indivs)
pop_adult = int(pop_breakdown[1]*indivs)
pop_elder = int(pop_breakdown[2]*indivs)

# FI (physical AND non physical) contact matrix taken from Mossong et al Finland: https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050074#pmed-0050074-st005
contact_matrix_fi = np.array([[21.21, 12.42, 0.85], # avg number contacts per day kids make with [0] kids, [1] adults, and [2] elderly
                  [16.44, 91.40, 6.42], # avg number contacts per day adults make with [0] kids, [1] adults, and [2] elderly
                  [1.04, 4.95, 4.15]]) # avg number contacts per day elderly make with [0] kids, [1] adults, and [2] elderly

# USA (physical AND non physical) contact matrix taken from Prem et al: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697#sec002
contact_matrix_usa = np.array([[22.2773775, 15.13427685, 0.67475806], # number contacts per day kids make with [0] kids, [1] adults, and [2] elderly
                            [19.65830386, 121.0409943, 2.531342454], # number contacts per day adults make with [0] kids, [1] adults, and [2] elderly
                            [2.730939972, 10.13063602, 3.481243491]])# number contacts per day elderly make with [0] kids, [1] adults, and [2] elderly

# France (physical AND non physical) contact matrix taken from Prem et al: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005697#sec002
contact_matrix_fr = np.array([[23.9255104, 28.52664938, 0.833780374], # number contacts per day kids make with [0] kids, [1] adults, and [2] elderly
                           [14.17373913, 102.3279637, 2.785170733],# number contacts per day adults make with [0] kids, [1] adults, and [2] elderly
                           [2.488764893, 9.938006315, 4.714598106]])# number contacts per day elderly make with [0] kids, [1] adults, and [2] elderly


if country == 'usa':
    contact_matrix = contact_matrix_usa
if country == 'fi':
    contact_matrix = contact_matrix_fi
else:
    contact_matrix = contact_matrix_fr

full_pop = np.empty([indivs, num_sero])
for i in range(num_sero):
    full_pop[:,i] = np.array([np.random.choice([int(0),int(1)], size = indivs, replace = True, p = [0.99, 0.01])])
  
carriers_all = pd.DataFrame(data = full_pop, columns = colnames).astype(np.int8)
#carriers_all1 = carriers_all.astype(int) # change floating 0s and 1s to integers
prev_col = carriers_all.copy(deep = True) # serotype-specific previous colonizations
ind_lambdas = carriers_all.copy(deep = True) # frame for individual lambda rates 
susc_full = max_carr_cap - carriers_all.sum(axis = 1)
carriers_all['susc'] = susc_full
agegroup = np.repeat([int(0), int(1), int(2)], [pop_child, pop_adult, pop_elder])
carriers_all['agegroup'] = agegroup
age = random.choices(range(0,6206), k=pop_child)+random.choices(range(6206,23361), k=pop_adult)+random.choices(range(23361,36500), k=pop_elder)
#np.concatenate((np.random.choice(range(0, 6206), size = pop_child, replace = True), np.random.choice(range(6206,23361), size = pop_adult, replace = True),np.random.choice(range(23361,36500), size = pop_elder, replace = True)), axis=None)
carriers_all['age'] = age

# change previous colonisations for adults and elderly according to child to adult simulations' distributions
prevcolprob = pd.read_csv(r"//qdrive.dide.ic.ac.uk/homes/al817/Technical/Python/Prevcolhist/prevcolprob_allsims.csv")
prevcolprob.loc[np.where(prevcolprob.serotype == '15B.C')[0], 'serotype'] = '15B/C'
prevcolprob.loc[np.where(prevcolprob.serotype == '6A.C')[0], 'serotype'] = '6A/C'
for i in colnames:
    df = prevcolprob[prevcolprob.serotype == i]
    prev_col.loc[np.where((carriers_all['age'] >= 6205) & (carriers_all['age'] < 9125))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 6205) & (carriers_all['age'] < 9125))[0]), replace = True, p = df.m20) 
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
    prev_col.loc[np.where((carriers_all['age'] >= 27375) & (carriers_all['age'] < 29200))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 27375) & (carriers_all['age'] < 29200))[0]), replace = True, p = df.m75)
    prev_col.loc[np.where((carriers_all['age'] >= 29200) & (carriers_all['age'] < 31025))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 29200) & (carriers_all['age'] < 31025))[0]), replace = True, p = df.m80)
    prev_col.loc[np.where((carriers_all['age'] >= 31025) & (carriers_all['age'] < 32850))[0], i] = np.random.choice(df.prevcol, size = len(np.where((carriers_all['age'] >= 31025) & (carriers_all['age'] < 32850))[0]), replace = True, p = df.m85)
    prev_col.loc[np.where(carriers_all['age'] >= 32850)[0], i] = np.random.choice(df.prevcol, size = len(np.where(carriers_all['age'] >= 32850)[0]), replace = True, p = df.m90)

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

#### IMMUNITY AND COMPETITION PARAMS taken from Cobey & Lipsitch
sigma = 0.4 # specific immunity
epsilon = 0.1 # non-specific immunity
theta = 0.01 # competition
mincarrdur = 5 # minimum carriage duration

# vaccine params
vacc_start = 7300 # start of vaccination period
vacc_dur = 10 # years of vaccine-induced protection
vacc_dur = vacc_dur*365 # protection converted to days 
vacc_eff = 1

carriers = carriers_all.copy(deep = True)
t = 0
stop_t = vacc_start # 20 years to reach SS without vacc
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

for t in range(stop_t+1):
#while t <= stop_t and tot_lambda1 != 0:
    
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
    
    # turn on vaccination after 20 years SS
#     if t >= 7300:
#         index_vacc = np.where((carriers.age >= 730) & (carriers.age <= 1095) & (carriers.vacc == False)) # index of ppl to be vaccinated
#         num_newvacc = round(vacc_eff*len(index_vacc[0])) # number of new people to vaccinate (= vacc efficacy * tot number of people)
#         adjusted_ind = index_vacc[0][0:num_newvacc] # indices of people for whom vaccination works/is protective
#         carriers.loc[adjusted_ind, 'vacc'] = True
        
    # migration: add 1 carrier to serotypes with zero carriers
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
    
    if np.any(totFOIbysero < 0):
            print('NEGATIVE FOI')
            print("sero:", seroparams.serotype[np.where(totFOIbysero < 0)[1]])
            print("indiv:", np.where(carriers_all[seroparams.serotype[np.where(totFOIbysero < 0)[1]]] > 0))
    
    # individual rates                
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
            # if person is carrying, then recovery rate as function of prevcol
            else:
                ind_lambdas.iloc[s,j] = (1/np.random.gamma(shape = alpha, scale = (A*np.exp(-prevcol*B))/alpha))                
    
    # tot.lambda = total rate of all/any events
    numindivspertime = round(indivs/50)
    i = 0
    U0 = np.random.uniform(0,1, size = numindivspertime)#20)
    U3 = np.random.uniform(0,1, size = numindivspertime)#20)
    
    #while i < numindivspertime: # update 20 individuals each time step for indivs = 1000
    for i in range(numindivspertime):
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
        
        # switching event (S --> I or I --> S)
        if which_indices[0].size > 0:
            which_sero = which_indices[0][0]
            carriers.iloc[which_person, which_sero] = ~carriers.iloc[which_person, which_sero] + 2    

        #print(which_sero)
        # reupdate the susc score
        new_susc = max_carr_cap - sum(carriers.iloc[which_person, 0:num_sero])
    
        # update colonisation history if person acquired new serotype unless person is elderly  
        if (new_susc < carriers.susc.loc[which_person]) & (carriers.age.loc[which_person] <= 23360): 
            prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        #i += 1

    # update time and age
    #t += 1
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
disinc['disease'] = disease_cas/100000

# Disease Incidence
index_children = np.where(carriers.age <= 6205)[0]
index_adults = np.where((carriers.age > 6205) & (carriers.age <= 23360))[0]
index_elderly = np.where(carriers.age > 23360)[0]

pop_child = len(index_children)
pop_adult = len(index_adults)
pop_elder = len(index_elderly)

# Carriage by serotype
carrprevchildsero = np.array(carriers[agegroup == 0].iloc[:,0:num_sero].sum(axis = 0))/pop_child
carrprevadultsero = np.array(carriers[agegroup == 1].iloc[:,0:num_sero].sum(axis = 0))/pop_adult
carrpreveldersero = np.array(carriers[agegroup == 2].iloc[:,0:num_sero].sum(axis = 0))/pop_elder

disinc['carrprevchild'] = carrprevchildsero
disinc['carrprevadult'] = carrprevadultsero
disinc['carrprevelder'] = carrpreveldersero

# disease cases by serotype
disease_child = round(disinc['carrprevchild']*disinc['invasiveness']*(t/365)*pop_child)
disease_adult = round(disinc['carrprevadult']*disinc['invasiveness']*(t/365)*pop_adult)
disease_elder = round(disinc['carrprevelder']*disinc['invasiveness']*(t/365)*pop_elder)

# save to file destination:
newfiledest = filedest + endstring + "disinc" + ".csv"
disinc.to_csv(path_or_buf= newfiledest, index = True)

# overall IPD incidence by age group
overalldischild = sum(disease_child)/(pop_child*(t/365))*100000
overalldisadult = sum(disease_adult)/(pop_adult*(t/365))*100000
overalldiselder = sum(disease_elder)/(pop_elder*(t/365))*100000

# overall carr prev by age group
overallcarrprevchild = carriers[agegroup == 0].iloc[:,0:num_sero].sum().sum()/pop_child #.605
overallcarrprevadult = carriers[agegroup == 1].iloc[:,0:num_sero].sum().sum()/pop_adult #.607
overallcarrprevelder = carriers[agegroup == 2].iloc[:,0:num_sero].sum().sum()/pop_elder #.5

agegrpstats = np.array([[overallcarrprevchild, overallcarrprevadult, overallcarrprevelder],
                        [overalldischild, overalldisadult, overalldiselder]])

agegrpstats = pd.DataFrame(agegrpstats, columns = ['child', 'adult', 'elder'], index = ['carrprev', 'disease inc'])

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
    ind_lambdas[non_index] = 1/(mincarrdur+(((1/ind_lambdas[non_index])-mincarrdur)*math.exp(-eps*sum(prev_col))))
    return(ind_lambdas)


vacc_stop_t = t + 1825 # current time + 20 years

# outputs for graph of time vs infected ppl
currentLength = chunkSize
plotsy2 = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = DS_inf_nams)
plotsy2.insert(loc = 0, column = 'time', value = np.nan)

next_ev = 1
start2 = time.time()

#while t <= vacc_stop_t and tot_lambda1 != 0:
for t in range(t+1, vacc_stop_t + 1):
    
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
            prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        #i += 1

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

phase = 'postvacc'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

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
index_children = np.where(carriers.age <= 6205)[0]
index_adults = np.where((carriers.age > 6205) & (carriers.age <= 23360))[0]
index_elderly = np.where(carriers.age > 23360)[0]

pop_child = len(index_children)
pop_adult = len(index_adults)
pop_elder = len(index_elderly)

sim_output2 = pd.DataFrame({'key':sim_output2.index, 'value':sim_output2.values})
seroprev2 = sim_output2.drop([0, num_sero+1])
seroprev2 = seroprev2.rename(columns = {'key': 'Serotype', 'value': 'carrprev'})
seroprev2['carrprev'] = seroprev2['carrprev']/indivs # overall carr prev by serotype (regardless of age group)

### Get disease incidence using invasiveness
disinc2 = seroprev2.merge(country_inv, on = 'Serotype')
disease_cas2 = disinc2['carrprev']*disinc2['invasiveness']
disinc2['disease'] = disease_cas2*indivs*(t/365)/100000 # overall disease incidence by serotype (regardless of age group)

# carriage prevalence by serotype and age grp
carrprevchildsero2 = np.array(carriers[agegroup == 0].iloc[:,0:num_sero].sum(axis = 0))/pop_child
carrprevadultsero2 = np.array(carriers[agegroup == 1].iloc[:,0:num_sero].sum(axis = 0))/pop_adult
carrpreveldersero2 = np.array(carriers[agegroup == 2].iloc[:,0:num_sero].sum(axis = 0))/pop_elder

disinc2['carrprevchild'] = carrprevchildsero2
disinc2['carrprevadult'] = carrprevadultsero2
disinc2['carrprevelder'] = carrpreveldersero2

# disease cases per serotype and age grp
disease_child2 = round(disinc2['carrprevchild']*disinc2['invasiveness']*(t/365)*pop_child)
disease_adult2 = round(disinc2['carrprevadult']*disinc2['invasiveness']*(t/365)*pop_adult)
disease_elder2 = round(disinc2['carrprevelder']*disinc2['invasiveness']*(t/365)*pop_elder)

# save to file destination:
newfiledest2 = filedest + endstring + "disinc" + ".csv"
disinc2.to_csv(path_or_buf= newfiledest2, index = True)

# overall age group stats (all serotypes): carriage prevalence and IPD incidence in each age group
overalldischild2 = sum(disease_child2)/(pop_child*(t/365))*100000
overalldisadult2 = sum(disease_adult2)/(pop_adult*(t/365))*100000
overalldiselder2 = sum(disease_elder2)/(pop_elder*(t/365))*100000

overallcarrprevchild2 = carriers[agegroup == 0].iloc[:,0:num_sero].sum().sum()/pop_child #.605
overallcarrprevadult2 = carriers[agegroup == 1].iloc[:,0:num_sero].sum().sum()/pop_adult #.607
overallcarrprevelder2 = carriers[agegroup == 2].iloc[:,0:num_sero].sum().sum()/pop_elder #.5

agegrpstats2 = np.array([[overallcarrprevchild2, overallcarrprevadult2, overallcarrprevelder2],
                        [overalldischild2, overalldisadult2, overalldiselder2]])

agegrpstats2 = pd.DataFrame(agegrpstats2, columns = ['child', 'adult', 'elder'], index = ['carrprev', 'disease inc'])

newfiledest = filedest + endstring + "agegrpstats" + ".csv"
agegrpstats2.to_csv(path_or_buf= newfiledest, index = True)

#################################################################################################################################################################
############# PART 3: POST-VACCINATION PERIOD 2 - extension of PCV ##############################################################################################
#################################################################################################################################################################

vt = np.concatenate((vt, new_vt))
vt_ind = np.concatenate((vt_ind, new_vt_ind))

vacc2_stop_t = t + 1825 # current time + 20 years

# outputs for graph of time vs infected ppl
currentLength = chunkSize
plotsy3 = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = DS_inf_nams)
plotsy3.insert(loc = 0, column = 'time', value = np.nan)

next_ev = 1
start3 = time.time()
#while t <= vacc2_stop_t and tot_lambda1 != 0:
for t in range(t+1, vacc2_stop_t+1):
    
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
            prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
        #i += 1

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

phase = 'postvacc_extendedPCV'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

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
index_children = np.where(carriers.age <= 6205)[0]
index_adults = np.where((carriers.age > 6205) & (carriers.age <= 23360))[0]
index_elderly = np.where(carriers.age > 23360)[0]

pop_child = len(index_children)
pop_adult = len(index_adults)
pop_elder = len(index_elderly)

sim_output3 = pd.DataFrame({'key':sim_output3.index, 'value':sim_output3.values})
seroprev3 = sim_output3.drop([0, num_sero+1])
seroprev3 = seroprev3.rename(columns = {'key': 'Serotype', 'value': 'carrprev'})
seroprev3['carrprev'] = seroprev3['carrprev']/indivs

### Get disease incidence using invasiveness
disinc3 = seroprev3.merge(country_inv, on = 'Serotype')
disease_cas3 = disinc3['carrprev']*disinc3['invasiveness']
disinc3['disease'] = disease_cas3*indivs*(t/365)/100000

# overall carriage prevalence and IPD by age group
carrprevchildsero3 = np.array(carriers[agegroup == 0].iloc[:,0:num_sero].sum(axis = 0))/pop_child
carrprevadultsero3 = np.array(carriers[agegroup == 1].iloc[:,0:num_sero].sum(axis = 0))/pop_adult
carrpreveldersero3 = np.array(carriers[agegroup == 2].iloc[:,0:num_sero].sum(axis = 0))/pop_elder

disinc3['carrprevchild'] = carrprevchildsero3
disinc3['carrprevadult'] = carrprevadultsero3
disinc3['carrprevelder'] = carrpreveldersero3

# number of disease cases of each serotype
disease_child3 = round(disinc3['carrprevchild']*disinc3['invasiveness']*(t/365)*pop_child)
disease_adult3 = round(disinc3['carrprevadult']*disinc3['invasiveness']*(t/365)*pop_adult)
disease_elder3 = round(disinc3['carrprevelder']*disinc3['invasiveness']*(t/365)*pop_elder)

# save to file destination:
newfiledest = filedest + endstring + "disinc" + ".csv"
disinc3.to_csv(path_or_buf= newfiledest, index = True)

# overall age group stats (allserotypes): carr prev and IPD incidence for each age group
overalldischild3 = sum(disease_child3)/(pop_child*(t/365))*100000
overalldisadult3 = sum(disease_adult3)/(pop_adult*(t/365))*100000
overalldiselder3 = sum(disease_elder3)/(pop_elder*(t/365))*100000

overallcarrprevchild3 = carriers[agegroup == 0].iloc[:,0:num_sero].sum().sum()/pop_child #.605
overallcarrprevadult3 = carriers[agegroup == 1].iloc[:,0:num_sero].sum().sum()/pop_adult #.607
overallcarrprevelder3 = carriers[agegroup == 2].iloc[:,0:num_sero].sum().sum()/pop_elder #.5

agegrpstats3 = np.array([[overallcarrprevchild3, overallcarrprevadult3, overallcarrprevelder3],
                        [overalldischild3, overalldisadult3, overalldiselder3]])

agegrpstats3 = pd.DataFrame(agegrpstats3, columns = ['child', 'adult', 'elder'], index = ['carrprev', 'disease inc'])

newfiledest3 = filedest + endstring + "agegrpstats" + ".csv"
agegrpstats3.to_csv(path_or_buf= newfiledest3, index = True)

####### Testing with dummy carriage data ########################################################################################################################

# taken from this study https://academic.oup.com/jid/article/184/4/451/809030 
#serogroup = np.array(['6B', '23F', '19F', '6A', '11', '14', '15', '35', 'R', 'NT'])
#num_isolates = np.array([213, 268, 258, 140, 116, 67, 56, 49, 55, 17])
#list_of_tuples = list(zip(serogroup, num_isolates))
#carrprevdf = pd.DataFrame(list_of_tuples, columns = ['Serotype', 'Carr_isolates'])
#carrprevdf['Carr_prev'] = carrprevdf.Carr_isolates/1530 #1530 is the total number of swabs from the study

#################################################################################################################################################################

