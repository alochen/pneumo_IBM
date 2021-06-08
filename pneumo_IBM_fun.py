########## PNEUMO IBM with fitted sero params imported and with synchronous time where 20 events occur at each ############
########## time step, and a migration rate such that when a serotype goes extinct, it is reseeded into the     ############
########## population in the form of one infectious carrier                                                    ############

import pandas as pd
import numpy as np
import math
import random
import time
import pneumo_shared as sh
# import glob
# import os
from builtins import dict

# function that computes new lambdas based on serotype switching event = prevax period
def new_lambdas_nv(carrier_person, ind_lambdas, prev_col):       
    susc = carrier_person.susc
    vacc = carrier_person.vacc
    num_sero = sh.num_sero
    max_carr_cap = sh.max_carr_cap
    
    index = np.where(carrier_person[0:num_sero] == 0) # index of serotypes not carried
    non_index = np.bool_(carrier_person[0:num_sero]) # index of serotypes carried
    
    eps = sh.epsilon#np.random.normal(epsilon)
    thet = sh.theta#np.random.normal(theta)
    sig = sh.sigma#np.random.normal(sigma)
    # if you're carrying max number of serotypes, then you don't feel any FOI for new serotype (regardless of vaccination)
    if susc <= 0:
        ind_lambdas[index[0]] = 0
        
    # otherwise, your FOI changes based on your colonization history and competition    
    else:
        ind_lambdas[index[0]] = (1-(sig*(prev_col[index[0]] > 0)))*ind_lambdas[index[0]]*(1-thet*(susc < max_carr_cap))   

    # serotypes you are carrying change recovery rate based on colonization history and immunity
    ind_lambdas[non_index] = 1/(mincarrdur+(((1/ind_lambdas[non_index])-mincarrdur)*math.exp(-eps*sum(prev_col))))
    return(ind_lambdas)


# function that computes new lambdas based on serotype switching event - vaccine period

def new_lambdas_v(carrier_person, ind_lambdas, prev_col):      
    susc = carrier_person.susc
    vacc = carrier_person.vacc
    num_sero = sh.num_sero
    max_carr_cap = sh.max_carr_cap
    mincarrdur = sh.mincarrdur
    
    index = np.where(carrier_person[0:num_sero] == 0) # index of serotypes not carried
    non_index = np.bool_(carrier_person[0:num_sero]) # index of serotypes carried
    vt_index = index[0][np.isin(carrier_person.index[index], vt)] # index of VT serotypes not carried 
    # draw from distribution of parameters with mean @ max likelihood
    eps = sh.epsilon#np.random.normal(epsilon)
    thet = sh.theta#np.random.normal(theta)
    sig = sh.sigma#np.random.normal(sigma)
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

def pneumoIBMrun(t0, stop_t, carriers, ind_lambdas, prev_col):
    global seroparams, age_breakdown, sigma, epsilon, theta, mincarrdur, vacc_dur, beta
    seroparams = sh.seroparams
    age_breakdown = sh.age_breakdown
    #carriers = sh.carriers
    #ind_lambdas = sh.ind_lambdas
    #prev_col = sh.prev_col
    num_sero = sh.num_sero
    vacc_start = sh.vacc_start
    vacc_dur = sh.vacc_dur
    indivs = sh.indivs
    sigma = sh.sigma
    epsilon = sh.epsilon # non-specific immunity
    theta = sh.theta # competition
    mincarrdur = sh.mincarrdur # minimum carriage duration
    vacc_dur = sh.vacc_dur # years of vaccine-induced protection
    vacc_eff = sh.vacc_eff # vaccine efficacy
    beta = sh.beta # transmission coefficient
    max_carr_cap = sh.max_carr_cap
    contact_matrix = sh.contact_matrix
    
    
   
    # outputs for graph of time vs infected ppl
    chunkSize = 100000
    currentLength = chunkSize
    DS_inf_nams = np.unique(seroparams['serotype'])#list(map(lambda x: 'infxtd' + str(x+1), range(num_sero)))
    plotsy = pd.DataFrame(data = np.nan, index = range(0,chunkSize), columns = DS_inf_nams)
    plotsy.insert(loc = 0, column = 'time', value = np.nan)
        
    next_ev = 1
    start = time.time()

    for t in range(t0, stop_t+1):
      
        # death rate = birth rate
        if len(np.where(carriers.age > age_breakdown[2])[0]) > 0:
            ind_to_die = np.where(carriers.age > age_breakdown[2])[0]
            carriers.loc[ind_to_die] = 0 # making this person a baby now
            carriers.loc[ind_to_die, 'susc'] = max_carr_cap # susc = max carrying capacity when nothing is carried
            carriers.loc[ind_to_die, 'vacc'] = False
            prev_col.loc[ind_to_die] = 0
            
        index_children = np.where(carriers.age < age_breakdown[0])[0]
        index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
        index_elderly = np.where(carriers.age >= age_breakdown[1])[0]
        
        pop_child = len(index_children)
        pop_adult = len(index_adults)
        pop_elder = len(index_elderly)
        
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
        
            #if t >= 7300:
        if t >= vacc_start:
            index_vacc = np.where((carriers.age >= 730) & (carriers.age <= 1095) & (carriers.vacc == False)) # index of ppl to be vaccinated
            num_newvacc = round(vacc_eff*len(index_vacc[0])) # number of new people to vaccinate (= vacc efficacy * tot number of people)
            adjusted_ind = index_vacc[0][0:num_newvacc] # indices of people for whom vaccination works/is protective
            carriers.loc[adjusted_ind, 'vacc'] = True
        
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
            if t < vacc_start:
                person_lambdas = new_lambdas_nv(carriers.loc[which_person], ind_lambdas.loc[which_person], prev_col.loc[which_person])
            else:
                person_lambdas = new_lambdas_v(carriers.loc[which_person], ind_lambdas.loc[which_person], prev_col.loc[which_person])          
            
            
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
            if (new_susc < carriers.susc.loc[which_person]):
                if (carriers.age.loc[which_person] <= age_breakdown[1]): 
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
        totchildcarr[extinctsero] = 1
        totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
        toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)
        
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
                    if t >= vacc_start and carriers.vacc[s] == True and np.isin(j,vt_ind):
                        ind_lambdas.iloc[s,j] = 1/vacc_dur
                    # if serotype is NVT or if person is not vaccinated
                    else:
                        ind_lambdas.iloc[s,j] = totFOIbysero[carriers['agegroup'][s],j]
    
    
        # update time and age
        next_ev += 1
        carriers.age += 1
        
    end = time.time()
    total_time = end - start
    plotsy = plotsy.dropna()
    sim_output = plotsy.loc[len(plotsy) - 1]
    sim_output['sim_time'] = total_time/3600    
    d = {'carriers': carriers, 'prev_col': prev_col, 'ind_lambdas':ind_lambdas, 'plotsy': plotsy, 'sim_output': sim_output}
    return (d)

