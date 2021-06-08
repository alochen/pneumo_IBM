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
import pneumo_IBM_fun
import pneumo_shared as sh

country = sh.country
indivs = sh.indivs
globalinv = True
country_inv = sh.country_inv
colnames = sh.colnames
num_sero = sh.num_sero
age_breakdown = sh.age_breakdown
max_carr_cap = sh.max_carr_cap

filedest = "//qdrive/homes/al817/Technical/Python/popruns/new May 2021 log files/" + country + '_' + globalinv*("glob") + "inv/" 

phase = 'prevacc'
mylist = [f for f in glob.glob(filedest + "*.csv")]
mylist = [fname.rsplit('.', 1)[0] for fname in os.listdir(filedest)]
filenums = [item.split('_',1)[0] for item in mylist]
filenums = [item.split('simulation',1)[-1] for item in filenums]
filenums = [int(i) for i in filenums] 
simnum = max(filenums) + 1
simnum = str(simnum) + '_' + globalinv*("glob") + "inv" + "pop" + str(indivs)
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_" + 'initFOIbeta_' #+ 'vaccdur15_'
#seroparams.to_csv(path_or_buf = filedest + endstring + 'dummy.csv', index = False)

# log file of main params
paramlogfile = dict(simnum = max(filenums) + 1, country = sh.country, beta = sh.beta, indivs = sh.indivs, vacc_dur = sh.vacc_dur, vacc_eff = sh.vacc_eff,
                    sigma = sh.sigma, epsilon = sh.epsilon, theta = sh.theta)
paramdf = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in paramlogfile.items() ])) 
filedestparams = filedest + endstring + 'paramlogfile.csv'
paramdf.to_csv(path_or_buf = filedestparams, index = False)

# PCV13
vt = np.array(country_inv.nlargest(13, 'invasiveness').Serotype)
vt_ind = np.where(np.isin(colnames, vt))[0]
# PCV20
PCV20df = country_inv[country_inv['Serotype'].isin(vt) == False]
new_vt = np.array(PCV20df.nlargest(7, 'invasiveness').Serotype)
new_vt_ind = np.where(np.isin(colnames, new_vt))[0]
# PCV30
PCV30df = PCV20df[PCV20df['Serotype'].isin(new_vt) == False]
new_vt2 = np.array(PCV30df.nlargest(10, 'invasiveness').Serotype)
new_vt_ind2 = np.where(np.isin(colnames, new_vt2))[0]

# top 13 globally invasive serotypes get vaccinated against
if globalinv == True:
    global_inv = pd.read_csv(r"//qdrive/homes/al817/Technical/R/Case-to-carrier/consol.child-new.csv")
    global_inv = global_inv.drop(['Unnamed: 0', 'Serogroup'], axis = 1)
    global_inv.columns = country_inv.columns
    # PCV13
    vt = np.array(global_inv.nlargest(13, 'invasiveness').Serotype)
    vt_ind = np.where(np.isin(colnames, vt))[0]
    # PCV20
    PCV20df = global_inv[global_inv['Serotype'].isin(vt) == False]
    new_vt = np.array(PCV20df.nlargest(7, 'invasiveness').Serotype)
    new_vt_ind = np.where(np.isin(colnames, new_vt))[0]
    # PCV30
    PCV30df = PCV20df[PCV20df['Serotype'].isin(new_vt) == False]
    new_vt2 = np.array(PCV30df.nlargest(10, 'invasiveness').Serotype)
    new_vt_ind2 = np.where(np.isin(colnames, new_vt2))[0]
    
# print VT into dataframe
vt_csv =  filedest + endstring + "VTsero" + ".csv"
d = dict(vt_1 = vt, vt_2 = new_vt, vt_3 = new_vt2)  
vtdf = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in d.items() ])) 
vtdf.to_csv(path_or_buf = vt_csv, index = False)

#################################################################################################################################################################
############# PART 1: PRE-VACCINATION PERIOD ####################################################################################################################
#################################################################################################################################################################

t = 0
stop_t = sh.vacc_start

# run IBM simulation
prevacc = pneumo_IBM_fun.pneumoIBMrun(t0 = t, stop_t = stop_t, carriers = sh.carriers, 
                                      ind_lambdas = sh.ind_lambdas, prev_col = sh.prev_col)

carriers = prevacc['carriers']
prev_col = prevacc['prev_col']
ind_lambdas = prevacc['ind_lambdas']
plotsy = prevacc['plotsy']
sim_output = prevacc['sim_output']

plotsy_melt= pd.melt(plotsy, id_vars = 'time')
plot = plotsy_melt.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

# carriage of each sero over time
plotsy_meltdest = filedest + endstring + 'carrdt.csv'
plotsy_melt.to_csv(path_or_buf= plotsy_meltdest, index = True)

# tot serotypes being carried at final time
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
disinc['disease_inc'] = (disease_cas/indivs)*100000 #dis per 100,000 ppl

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

# save to file destination:
newfiledest = filedest + endstring + "disinc" + ".csv"
disinc.to_csv(path_or_buf= newfiledest, index = True)

# disease cases by serotype
disease_child = disinc['carrprevchild']*disinc['invasiveness']*(stop_t/365)*pop_child
disease_adult = disinc['carrprevadult']*disinc['invasiveness']*(stop_t/365)*pop_adult
disease_elder = disinc['carrprevelder']*disinc['invasiveness']*(stop_t/365)*pop_elder

# overall IPD incidence by age group
overalldis = ((disease_child + disease_adult + disease_elder).sum()/(indivs*(stop_t/365)))*100000
overalldischild = (sum(disease_child)/(pop_child*(stop_t/365)))*100000
overalldisadult = (sum(disease_adult)/(pop_adult*(stop_t/365)))*100000
overalldiselder = (sum(disease_elder)/(pop_elder*(stop_t/365)))*100000

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
############# PART 2: POST-VACCINATION PERIOD PCV13 #################################################################################################################
#################################################################################################################################################################

vacc_stop_t = stop_t + 3650#1825 # current time + 5 years

# run IBM simulation
postPCV13 = pneumo_IBM_fun.pneumoIBMrun(t0 = stop_t + 1, stop_t = vacc_stop_t, carriers = carriers, 
                                      ind_lambdas = ind_lambdas, prev_col = prev_col)

carriers = postPCV13['carriers']
prev_col = postPCV13['prev_col']
ind_lambdas = postPCV13['ind_lambdas']
plotsy2 = postPCV13['plotsy']
sim_output2 = postPCV13['sim_output']

phase = 'postPCV13'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

plotsy_melt2= pd.melt(plotsy2, id_vars = 'time')
plot2 = plotsy_melt2.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

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
tprime = vacc_stop_t - stop_t # new period
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

# save to file destination:
newfiledest2 = filedest + endstring + "disinc" + ".csv"
disinc2.to_csv(path_or_buf= newfiledest2, index = True)

# disease cases by serotype and age group
disease_child2 = disinc2['carrprevchild']*disinc2['invasiveness']*(tprime/365)*pop_child
disease_adult2 = disinc2['carrprevadult']*disinc2['invasiveness']*(tprime/365)*pop_adult
disease_elder2 = disinc2['carrprevelder']*disinc2['invasiveness']*(tprime/365)*pop_elder

# overall IPD incidence by age group
overalldis2 = ((disease_child2 + disease_adult2 + disease_elder2).sum()/(indivs*(tprime/365)))*100000
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
############# PART 3: POST-VACCINATION PERIOD 2 - PCV20 ##############################################################################################
#################################################################################################################################################################

vt = np.concatenate((vt, new_vt))
vt_ind = np.concatenate((vt_ind, new_vt_ind))

vacc2_stop_t = vacc_stop_t + 3650 # current time + 10 years

# run IBM simulation
postPCV20 = pneumo_IBM_fun.pneumoIBMrun(t0 = vacc_stop_t + 1, stop_t = vacc2_stop_t, carriers = carriers, 
                                      ind_lambdas = ind_lambdas, prev_col = prev_col)

carriers = postPCV20['carriers']
prev_col = postPCV20['prev_col']
ind_lambdas = postPCV20['ind_lambdas']
plotsy3 = postPCV20['plotsy']
sim_output3 = postPCV20['sim_output']

plotsy_melt3= pd.melt(plotsy3, id_vars = 'time')
plot3 = plotsy_melt3.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

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

# save to file destination:
newfiledest = filedest + endstring + "disinc" + ".csv"
disinc3.to_csv(path_or_buf= newfiledest, index = True)

# disease cases by serotype and age group
disease_child3 = disinc3['carrprevchild']*disinc3['invasiveness']*(tprime/365)*pop_child
disease_adult3 = disinc3['carrprevadult']*disinc3['invasiveness']*(tprime/365)*pop_adult
disease_elder3 = disinc3['carrprevelder']*disinc3['invasiveness']*(tprime/365)*pop_elder

# overall IPD incidence by age group
overalldis3 = ((disease_child3 + disease_adult3 + disease_elder3).sum()/(indivs*(tprime/365)))*100000
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
############# PART 4: POST-VACCINATION PERIOD PCV30 ##############################################################################################
#################################################################################################################################################################

vt = np.concatenate((vt, new_vt2))
vt_ind = np.concatenate((vt_ind, new_vt_ind2))

vacc3_stop_t = t + 3650 # current time + 10 years

# run IBM simulation
postPCV30 = pneumo_IBM_fun.pneumoIBMrun(t0 = vacc2_stop_t + 1, stop_t = vacc3_stop_t, carriers = carriers, 
                                      ind_lambdas = ind_lambdas, prev_col = prev_col)

carriers = postPCV30['carriers']
prev_col = postPCV30['prev_col']
ind_lambdas = postPCV30['ind_lambdas']
plotsy4 = postPCV30['plotsy']
sim_output4 = postPCV30['sim_output']

plotsy_melt4= pd.melt(plotsy4, id_vars = 'time')
plot4 = plotsy_melt4.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

phase = 'postPCV30'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

plotsy_melt4dest = filedest + endstring + 'carrdt.csv'
plotsy_melt4.to_csv(path_or_buf= plotsy_meltdest, index = True)

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
disinc4['disease'] = (disease_cas4/indivs)*100000

# overall carriage prevalence and IPD by age group
carrprevchildsero4 = np.array(carriers[carriers.agegroup == 0].iloc[:,0:num_sero].sum(axis = 0))/pop_child
carrprevadultsero4 = np.array(carriers[carriers.agegroup == 1].iloc[:,0:num_sero].sum(axis = 0))/pop_adult
carrpreveldersero4 = np.array(carriers[carriers.agegroup == 2].iloc[:,0:num_sero].sum(axis = 0))/pop_elder

disinc4['carrprevchild'] = carrprevchildsero4
disinc4['carrprevadult'] = carrprevadultsero4
disinc4['carrprevelder'] = carrpreveldersero4

# save to file destination:
newfiledest = filedest + endstring + "disinc" + ".csv"
disinc4.to_csv(path_or_buf= newfiledest, index = True)

# disease cases by serotype
disease_child4 = disinc4['carrprevchild']*disinc4['invasiveness']*(tprime/365)*pop_child
disease_adult4 = disinc4['carrprevadult']*disinc4['invasiveness']*(tprime/365)*pop_adult
disease_elder4 = disinc4['carrprevelder']*disinc4['invasiveness']*(tprime/365)*pop_elder

# overall IPD incidence by age group
overalldis4 =  ((disease_child4 + disease_adult4 + disease_elder4).sum()/(indivs*(tprime/365)))*100000
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

#################################################################################################################################################################

# import win32com.client as win32
# from win32 import win32api
# outlook = win32.Dispatch('outlook.application')
# mail = outlook.CreateItem(0)
# mail.To = 'a.lochen17@imperial.ac.uk'
# mail.Subject = country + 'simulation ' + str(simnum) + ' has finished running'
# mail.Body = 'Congrats!'
# mail.Send()