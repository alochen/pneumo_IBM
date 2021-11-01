"""
Running pneumo IBM children to become adults

@author: al817
"""

import pandas as pd
import numpy as np
import math
import random
import time
import glob

mylist = [f for f in glob.glob("//qdrive/homes/al817/Technical/Python/Prevcolhist/sim*.csv")]
latestfile = max(mylist)
simnum = int(latestfile.rsplit('sim', 1)[1].rsplit('.csv',1)[0])+1
path = "//qdrive/homes/al817/Technical/Python/Prevcolhist/"
pd.Series(simnum).to_csv(path + 'sim' + str(simnum) + '.csv')

def new_lambdas(carrier_person, ind_lambdas, prev_col, sigma):    
    susc = carrier_person.susc

    index = np.where(carrier_person[0:num_sero] == 0) # index of serotypes not carried
    non_index = np.bool_(carrier_person[0:num_sero]) # index of serotypes carried
    # draw from distribution of parameters with mean @ max likelihood
    eps = epsilon#np.random.normal(epsilon)
    thet = theta#np.random.normal(theta)
    sig = sigma#np.random.normal(sigma)
    # if you're carrying max number of serotypes, then you don't feel any FOI for new serotype
    if susc <= 0:
        ind_lambdas[index[0]] = 0
    # otherwise, your FOI changes based on your colonization history and competition
    else:
        ind_lambdas[index[0]] = (1-(sig*(prev_col[index[0]] > 0)))*ind_lambdas[index[0]]*(1-thet*(susc < max_carr_cap))
    # serotypes you are carrying change FOI based on colonization history and immunity
    #ind_lambdas[non_index] = (1/250) + ((ind_lambdas[non_index]-(1/250))*math.exp(eps*sum(prev_col)))
    ind_lambdas[non_index] = 1/(mincarrdur+(((1/ind_lambdas[non_index])-mincarrdur)*math.exp(-eps*sum(prev_col))))
    return(ind_lambdas)


carr_dat = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/carriage_episodes_AL.csv")
#paramsdf = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/carrdurMCMCparams.csv")
seroparams = pd.read_csv(r"//qdrive/homes/al817/Technical/Python/seroparams.csv")
seroparams = seroparams.drop(seroparams.index[seroparams.serotype.isin(['15B/C_old','15B','15C'])]).reset_index().drop(['index'], axis = 1)

# # remove any rows with non positive length of carr
carr_dat = carr_dat.loc[np.where(carr_dat.length > 0)].reset_index(drop = True) 
carr_dat.t_f = carr_dat.t_i + carr_dat.length
# remove dates columns`
carr_dat = carr_dat.iloc[:,2:len(carr_dat.columns)+1]
# create ID instead of bday
carr_dat['ID'] = (carr_dat.sort_index().groupby(['bday'], sort=False).ngroup()+1)
carr_dat = carr_dat.drop(['bday','lanes'], axis = 1)
carr_dat = carr_dat.sort_values('t_i')
carr_dat['coln_num'] = carr_dat.groupby('ID').cumcount()
t = np.arange(0,max(carr_dat.t_f)+1)
max_t = max(t)
unique_sero = seroparams.serotype#np.unique(carr_dat.serotype)

indivs = len(np.unique(carr_dat['ID']))
pop_child = indivs
pop_adult = 0
pop_elder = 0
if pop_adult == 0:
    pop_adult = 0.0000000000000000000000000000000001
if pop_elder == 0:
    pop_elder = 0.0000000000000000000000000000000001
num_sero = len(seroparams.serotype)#len(unique_sero)
max_carr_cap = 2
children_data = np.empty([indivs, num_sero])
for i in range(num_sero):
    children_data[:,i] = np.array([np.random.choice([int(0),int(1)], size = indivs, replace = True, p = [0.99, 0.01])])

carriers_children = pd.DataFrame(data = children_data, columns = seroparams.serotype)
carriers_children = carriers_children.astype(int) # change floating 0s and 1s to integers
prev_col_children = carriers_children.copy(deep = True) # serotype-specific previous colonizations
susc_children = max_carr_cap - carriers_children.sum(axis = 1)
carriers_children['susc'] = susc_children
carriers_children['agegroup'] = int(0)
#carriers_children['age'] = carr_dat['age_d']
age = random.sample(range(0,365), indivs)
carriers_children['age'] = age
ind_lambdas = prev_col_children.copy(deep = True)

# contact_matrix = np.array([[21.21, 12.42, 0.85], # number contacts per day kids make with [0] kids, [1] adults, and [2] elderly
#                   [16.44, 91.40, 6.42], # number contacts per day adults make with [0] kids, [1] adults, and [2] elderly
#                   [1.04, 4.95, 4.15]]) # number contacts per day elderly make with [0] kids, [1] adults, and [2] elderly

# based on france
contact_matrix = np.array([[3.804920797, 6.691382506, 0.26971156],
                           [5.489241519, 152.9683178, 3.349239548],
                           [0.546826761, 12.75265203, 4.714598105]])


#B_rates_child_val = 0.021#[0.012/indivs]*num_sero # per day
#B_rates_adult_val = 0.0063 # children are 70% of FOI so adult FOI is reduced by this amount #[0.004/indivs]*num_sero # per day
#B_rates_elder_val = 0.021#[0.012/indivs]*num_sero # per day
#seroparams['beta_child'] = B_rates_child_val
#seroparams['beta_adult'] = B_rates_adult_val
#seroparams['beta_elder'] = B_rates_elder_val

sigma = 0.4 # specific immunity
epsilon = 0.1 # non-specific immunity
theta = 0.1 # competition
mincarrdur = 5
beta = 0.00012075

agecutoffs = [5, 65, 81.7]
age_breakdown = [round(i * 365) for i in agecutoffs] #[365*5, 23725, maxage*365]

carriers = carriers_children.copy(deep = True)
prev_col = prev_col_children.copy(deep = True) # serotype-specific previous colonizations
t = 0
stop_t = 33000#32850
tot_lambda1 = (sum(carriers.iloc[:,0])*seroparams['beta'][0])

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
    totFOIbysero_children[i] = ((contact_matrix[0][0]*beta*totchildcarr[i]/pop_child)+ 
                                (contact_matrix[0][1]*beta*totadultcarr[i]/pop_adult)+
                                (contact_matrix[0][2]*beta*toteldercarr[i]/pop_elder))
    totFOIbysero_adults[i] = ((contact_matrix[1][0]*beta*totchildcarr[i]/pop_child)+ 
                              (contact_matrix[1][1]*beta*totadultcarr[i]/pop_adult)+ 
                              (contact_matrix[1][2]*beta*toteldercarr[i]/pop_elder))
    totFOIbysero_elders[i] =  ((contact_matrix[2][0]*beta*totchildcarr[i]/pop_child)+ 
                               (contact_matrix[2][1]*beta*totadultcarr[i]/pop_adult)+ 
                               (contact_matrix[2][2]*beta*toteldercarr[i]/pop_elder))
   
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
    
# outputs for graph of time vs infected ppl
chunkSize = 100000
currentLength = chunkSize
plotsy = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = unique_sero)
plotsy.insert(loc = 0, column = 'time', value = np.nan)
    
next_ev = 1
start = time.time()
for t in range(stop_t+1):
    
    if len(np.where(carriers.age > age_breakdown[2])[0]) > 0:
        ind_to_die = np.where(carriers.age > age_breakdown[2])[0]
        carriers.loc[ind_to_die] = 0 # making this person a baby now
        carriers.loc[ind_to_die, 'susc'] = max_carr_cap # susc = max carrying capacity when nothing is carried
        prev_col.loc[ind_to_die] = 0
        
    index_children = np.where(carriers.age <= age_breakdown[0])[0]
    index_adults = np.where((carriers.age > age_breakdown[0]) & (carriers.age <= age_breakdown[1]))[0]
    index_elderly = np.where(carriers.age > age_breakdown[1])[0]
    
    pop_child = len(index_children)
    pop_adult = len(index_adults)
    pop_elder = len(index_elderly)
    if pop_adult == 0:
        pop_adult = 0.0000000000000000000000000000000001
    if pop_elder == 0:
        pop_elder = 0.0000000000000000000000000000000001
    
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
    i = 0
    U0 = np.random.uniform(0,1, size = 3)
    U3 = np.random.uniform(0,1, size = 3)
    while i < 3:
     #for i in range(3): # update 50 individuals each time step
        tot_lambda1 = ind_lambdas.sum().sum()
        U1 = U0[i]
        U1_tot_lambda1 = U1*tot_lambda1
        
        # choose person index ( = which person will this happen to?)
        cumsum_person = ind_lambdas.sum(axis = 1)
        try:
            which_person = np.where(U1_tot_lambda1 < np.cumsum(cumsum_person))[0][0]
        except IndexError:
            which_person = indivs - 1
        #max(np.where(U1_tot_lambda1 > np.cumsum(cumsum_person))[0]) 
        #print(which_person)

        # choose event that will take place (= which serotype will be toggled?)
        person_lambdas = new_lambdas(carriers.loc[which_person], ind_lambdas.loc[which_person], prev_col.loc[which_person], sigma)
        tot_lambda2 = sum(person_lambdas)
        U2 = U3[i]
        U2_tot_lambda2 = U2*tot_lambda2

        which_indices = np.where(U2_tot_lambda2 < np.cumsum(person_lambdas))
        
        # switching event (S --> I or I --> S)
        if which_indices[0].size > 0:
            which_sero = which_indices[0][0]
            carriers.iloc[which_person, which_sero] = ~carriers.iloc[which_person, which_sero] + 2

        # reupdate the susc score
        new_susc = max_carr_cap - sum(carriers.iloc[which_person, 0:num_sero])
    
        # update colonisation history if person acquired new serotype              
        if (new_susc < carriers.susc.loc[which_person]):
            #### the line below shouldn't be here which is why the prev col was so weird
                #if (carriers.age.loc[which_person] <= age_breakdown[1]):
                prev_col.iloc[which_person, which_sero] = prev_col.iloc[which_person, which_sero] + 1 # end of if statement
                alphatemp = seroparams['alpha_emcee'][which_sero]
                Atemp = seroparams['A_emcee'][which_sero]
                Btemp = seroparams['B_emcee'][which_sero]
                ind_lambdas.iloc[which_person, which_sero] = (1/np.random.gamma(shape = alphatemp,
                                                                                scale = (Atemp*np.exp(-prev_col.iloc[which_person, which_sero]*Btemp))/alphatemp))
                #ind_lambdas[which_person, which_sero] = (1/np.random.gamma(shape = alphatemp,
            
        carriers.loc[carriers.index == which_person,'susc'] = new_susc
 
        i += 1
    
    #migration: add 1 carrier to serotypes with zero carriers
    extinctsero = np.where(carriers.iloc[:,0:num_sero].sum(axis = 0) == 0)[0]
    totchildcarr = carriers.iloc[index_children,0:num_sero].sum(axis = 0)
    totchildcarr[extinctsero] = 1
    totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
    totadultcarr[extinctsero] = 1
    toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)
    toteldercarr[extinctsero] = 1
    
    totFOIbysero_children = [0]*num_sero
    totFOIbysero_adults = [0]*num_sero
    totFOIbysero_elders = [0]*num_sero
    
    for i in range(len(seroparams)):
        totFOIbysero_children[i] = ((contact_matrix[0][0]*beta*totchildcarr[i]/pop_child)+ 
                                    (contact_matrix[0][1]*beta*totadultcarr[i]/pop_adult)+
                                    (contact_matrix[0][2]*beta*toteldercarr[i]/pop_elder))
        totFOIbysero_adults[i] = ((contact_matrix[1][0]*beta*totchildcarr[i]/pop_child)+ 
                                  (contact_matrix[1][1]*beta*totadultcarr[i]/pop_adult)+ 
                                  (contact_matrix[1][2]*beta*toteldercarr[i]/pop_elder))
        totFOIbysero_elders[i] =  ((contact_matrix[2][0]*beta*totchildcarr[i]/pop_child)+ 
                                   (contact_matrix[2][1]*beta*totadultcarr[i]/pop_adult)+ 
                                   (contact_matrix[2][2]*beta*toteldercarr[i]/pop_elder))
       
    totFOIbysero = np.row_stack((totFOIbysero_children, totFOIbysero_adults, totFOIbysero_elders)) 
    
    # individual rates                
    for s in range(indivs):
        for j in range(num_sero):
            # if person not carrying
            if carriers.iloc[s,j] == 0:
                    ind_lambdas.iloc[s,j] = totFOIbysero[carriers['agegroup'][s],j]
        
    if t > 6205 and (t/365) % 5 == 0:
        newprevcol = prev_col.assign(age = carriers_children['age'], time = t)
        newprevcol.to_csv(path_or_buf = "//qdrive/homes/al817/Technical/Python/Prevcolhist/sim" + str(simnum) + '_prevcolhist_randage_' + str(t) + ".csv", index = True)

    # update age
    next_ev += 1
    carriers.age += 1#time_next
        
end = time.time()
total_time = end - start
plotsy = plotsy.dropna()
plotsy_melt= pd.melt(plotsy, id_vars = 'time')
#plot = plotsy_melt.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

# output file
sim_output = plotsy.loc[len(plotsy) - 1]
sim_output['sim_time'] = total_time/3600

newfiledest = "//qdrive/homes/al817/Technical/Python/Prevcolhist/child2adultsim" + str(simnum) + ".csv"
sim_output.to_csv(path_or_buf= newfiledest, index = True)

# prevcol
prevcolnam = "//qdrive/homes/al817/Technical/Python/Prevcolhist/child2adultprevcol" + str(simnum) + ".csv"
prev_col.to_csv(path_or_buf = prevcolnam, index = True)

# carriers
carrnam =  "//qdrive/homes/al817/Technical/Python/Prevcolhist/child2adult_carriers" + str(simnum) + ".csv"
carriers.to_csv(path_or_buf = carrnam, index = True)

# output figure
#plotname = "child2adultsim" + str(simnum) + ".pdf"

# comment this out for cluster
#plot.get_figure().savefig(plotname, bbox_inches = "tight")


