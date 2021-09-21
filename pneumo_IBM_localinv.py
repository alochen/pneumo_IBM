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
import traceback

# import params
country = sh.country
indivs = sh.indivs
country_inv = sh.country_inv
colnames = sh.colnames
num_sero = sh.num_sero
age_breakdown = sh.age_breakdown
max_carr_cap = sh.max_carr_cap
seroparams = sh.seroparams
contact_matrix = sh.contact_matrix
simnum = sh.simnum
beta = sh.beta
sigma = sh.sigma
epsilon = sh.epsilon
theta = sh.theta
vacc_dur = sh.vacc_dur
vacc_eff = sh.vacc_eff

filedest = "//qdrive/homes/al817/Technical/Python/popsims/" + country + "_inv/" 

phase = 'prevacc'
# mylist = [f for f in glob.glob(filedest + "*.csv")]
# mylist = [fname.rsplit('.', 1)[0] for fname in os.listdir(filedest)]
# filenums = [item.split('_',1)[0] for item in mylist]
# filenums = [item.split('simulation',1)[-1] for item in filenums]
# filenums = [int(i) for i in filenums] 
# simnum = max(filenums) + 1
string1 = str(simnum) + '_' + "inv" + "pop" + str(indivs)
endstring = "simulation" + string1 + "_" + phase + "_" + country + "_" #+ 'initFOIbeta_' #+ 'vaccdur15_'
#seroparams.to_csv(path_or_buf = filedest + endstring + 'dummy.csv', index = False)

# log file of main params
paramlogfile = dict(simnum = simnum, country = country, beta = beta, indivs = indivs, vacc_dur = vacc_dur, vacc_eff = vacc_eff,
                    sigma = sigma, epsilon = epsilon, theta = theta)
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
try: 
    prevacc = pneumo_IBM_fun.pneumoIBMrun(t0 = t, stop_t = stop_t, carriers = sh.carriers, 
                                          ind_lambdas = sh.ind_lambdas, prev_col = sh.prev_col, vt = np.nan)
except Exception: 
    error_message = traceback.print_exc()
    print('simulation '+ str(simnum))

carriers = prevacc['carriers']
prev_col = prevacc['prev_col']
ind_lambdas = prevacc['ind_lambdas']
plotsy = prevacc['plotsy']
sim_output = prevacc['sim_output']

# commented out bc of agg fun
plotsy_melt= pd.melt(plotsy, id_vars = 'time')
plot = plotsy.plot(kind = 'line')
#plot = plotsy_melt.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

# carriage of each sero over time
#plotsy_meltdest = filedest + endstring + 'carrdt.csv'
#plotsy_melt.to_csv(path_or_buf= plotsy_meltdest, index = True)

# tot serotypes being carried at final time
newfiledest = filedest + endstring +  "output" +".csv"
sim_output.to_csv(path_or_buf= newfiledest, index = True)

# prevcol
prevcolnam = filedest +  endstring + "prevcol" +".csv"
prev_col.to_csv(path_or_buf = prevcolnam, index = True)

# carriers
carrnam =  filedest + endstring + "carrierstmax" +".csv"
carriers.to_csv(path_or_buf = carrnam, index = True)

# output figure - COMMENTED OUT bc of agg fun error
plotname = filedest + endstring + "simulation" + ".pdf"
plot.get_figure().savefig(plotname, bbox_inches = "tight")

sim_output = pd.DataFrame({'key':sim_output.index, 'value':sim_output.values})
seroprev = sim_output.drop([0, num_sero+1])
seroprev = seroprev.rename(columns = {'key': 'Serotype', 'value': 'carrprev'})
seroprev['carrprev'] = seroprev['carrprev']/indivs

### DISEASE INCIDENCE USING INVASIVENESS
disinc = seroprev.merge(country_inv, on = 'Serotype')
disease_cas = disinc['carrprev']*disinc['invasiveness']*indivs*(stop_t/365)
disinc['disease_cas'] = disease_cas
disinc['disease_inc'] = (disease_cas/indivs)*100000 #dis per 100,000 ppl

disease_cas_lo = disinc['carrprev']*disinc['invasiveness.low']*indivs*(stop_t/365)
disease_cas_hi = disinc['carrprev']*disinc['invasiveness.high']*indivs*(stop_t/365)
disinc['disease_cas_lo'] = disease_cas_lo
disinc['disease_cas_hi'] = disease_cas_hi
disinc['disease_inc_lo'] = (disease_cas_lo/indivs)*100000 #dis per 100,000 ppl
disinc['disease_inc_hi'] = (disease_cas_hi/indivs)*100000 #dis per 100,000 ppl

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

disease_child_lo = disinc['carrprevchild']*disinc['invasiveness.low']*(stop_t/365)*pop_child
disease_adult_lo = disinc['carrprevadult']*disinc['invasiveness.low']*(stop_t/365)*pop_adult
disease_elder_lo = disinc['carrprevelder']*disinc['invasiveness.low']*(stop_t/365)*pop_elder
disease_child_hi = disinc['carrprevchild']*disinc['invasiveness.high']*(stop_t/365)*pop_child
disease_adult_hi = disinc['carrprevadult']*disinc['invasiveness.high']*(stop_t/365)*pop_adult
disease_elder_hi = disinc['carrprevelder']*disinc['invasiveness.high']*(stop_t/365)*pop_elder

# overall IPD incidence by age group
overalldis = ((disease_child + disease_adult + disease_elder).sum()/(indivs*(stop_t/365)))*100000
overalldischild = (sum(disease_child)/(pop_child*(stop_t/365)))*100000
overalldisadult = (sum(disease_adult)/(pop_adult*(stop_t/365)))*100000
overalldiselder = (sum(disease_elder)/(pop_elder*(stop_t/365)))*100000

overalldis_lo = ((disease_child_lo + disease_adult_lo + disease_elder_lo).sum()/(indivs*(stop_t/365)))*100000
overalldischild_lo = (sum(disease_child_lo)/(pop_child*(stop_t/365)))*100000
overalldisadult_lo = (sum(disease_adult_lo)/(pop_adult*(stop_t/365)))*100000
overalldiselder_lo = (sum(disease_elder_lo)/(pop_elder*(stop_t/365)))*100000
overalldis_hi = ((disease_child_hi + disease_adult_hi + disease_elder_hi).sum()/(indivs*(stop_t/365)))*100000
overalldischild_hi = (sum(disease_child_hi)/(pop_child*(stop_t/365)))*100000
overalldisadult_hi = (sum(disease_adult_hi)/(pop_adult*(stop_t/365)))*100000
overalldiselder_hi = (sum(disease_elder_hi)/(pop_elder*(stop_t/365)))*100000

overallcarrprev = len(np.where(carriers.susc < max_carr_cap)[0])/indivs
overallcarrprevchild = len(np.where(carriers[carriers.agegroup == 0].susc < max_carr_cap)[0])/pop_child
overallcarrprevadult = len(np.where(carriers[carriers.agegroup == 1].susc < max_carr_cap)[0])/pop_adult
overallcarrprevelder = len(np.where(carriers[carriers.agegroup == 2].susc < max_carr_cap)[0])/pop_elder

agegrpstats = np.array([[overallcarrprev, overallcarrprevchild, overallcarrprevadult, overallcarrprevelder],
                        [overalldis, overalldischild, overalldisadult, overalldiselder],
                        [overalldis_lo, overalldischild_lo, overalldisadult_lo, overalldiselder_lo],
                        [overalldis_hi, overalldischild_hi, overalldisadult_hi, overalldiselder_hi]])

agegrpstats = pd.DataFrame(agegrpstats, columns = ['overall','child', 'adult', 'elder'], 
                           index = ['carrprev', 'disease inc', 'disinc_lo', 'disinc_hi'])

newfiledest = filedest + endstring + "agegrpstats" + ".csv"
agegrpstats.to_csv(path_or_buf= newfiledest, index = True)

#################################################################################################################################################################
############# PART 2: POST-VACCINATION PERIOD PCV13 #################################################################################################################
#################################################################################################################################################################

vacc_stop_t = stop_t + 3650#1825 # current time + 5 years

# run IBM simulation
try:
    postPCV13 = pneumo_IBM_fun.pneumoIBMrun(t0 = stop_t + 1, stop_t = vacc_stop_t, carriers = carriers, 
                                            ind_lambdas = ind_lambdas, prev_col = prev_col, vt = vt)
except Exception: 
    error_message = traceback.print_exc()
    print('simulation '+ str(simnum))

carriers = postPCV13['carriers']
prev_col = postPCV13['prev_col']
ind_lambdas = postPCV13['ind_lambdas']
plotsy2 = postPCV13['plotsy']
sim_output2 = postPCV13['sim_output']

phase = 'postPCV13'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

plotsy_melt2= pd.melt(plotsy2, id_vars = 'time')
plot2 = plotsy2.plot(kind = 'line')
#plot2 = plotsy_melt2.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

#plotsy_melt2dest = filedest + endstring + 'carrdt.csv'
#plotsy_melt2.to_csv(path_or_buf= plotsy_melt2dest, index = True)

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

disease_cas_lo2 = disinc2['carrprev']*disinc2['invasiveness.low']*indivs*(tprime/365)
disease_cas_hi2 = disinc2['carrprev']*disinc2['invasiveness.high']*indivs*(tprime/365)
disinc2['disease_cas_lo'] = disease_cas_lo2
disinc2['disease_cas_hi'] = disease_cas_hi2
disinc2['disease_inc_lo'] = (disease_cas_lo2/indivs)*100000 #dis per 100,000 ppl
disinc2['disease_inc_hi'] = (disease_cas_hi2/indivs)*100000 #dis per 100,000 ppl

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

disease_child2_lo = disinc2['carrprevchild']*disinc2['invasiveness.low']*(tprime/365)*pop_child
disease_adult2_lo = disinc2['carrprevadult']*disinc2['invasiveness.low']*(tprime/365)*pop_adult
disease_elder2_lo = disinc2['carrprevelder']*disinc2['invasiveness.low']*(tprime/365)*pop_elder
disease_child2_hi = disinc2['carrprevchild']*disinc2['invasiveness.high']*(tprime/365)*pop_child
disease_adult2_hi = disinc2['carrprevadult']*disinc2['invasiveness.high']*(tprime/365)*pop_adult
disease_elder2_hi = disinc2['carrprevelder']*disinc2['invasiveness.high']*(tprime/365)*pop_elder

# overall IPD incidence by age group
overalldis2 = ((disease_child2 + disease_adult2 + disease_elder2).sum()/(indivs*(tprime/365)))*100000
overalldischild2 = (sum(disease_child2)/(pop_child*(tprime/365)))*100000
overalldisadult2 = (sum(disease_adult2)/(pop_adult*(tprime/365)))*100000
overalldiselder2 = (sum(disease_elder2)/(pop_elder*(tprime/365)))*100000

overalldis2_lo = ((disease_child2_lo + disease_adult2_lo + disease_elder2_lo).sum()/(indivs*(tprime/365)))*100000
overalldischild2_lo = (sum(disease_child2_lo)/(pop_child*(tprime/365)))*100000
overalldisadult2_lo = (sum(disease_adult2_lo)/(pop_adult*(tprime/365)))*100000
overalldiselder2_lo = (sum(disease_elder2_lo)/(pop_elder*(tprime/365)))*100000
overalldis2_hi = ((disease_child2_hi + disease_adult2_hi + disease_elder2_hi).sum()/(indivs*(tprime/365)))*100000
overalldischild2_hi = (sum(disease_child2_hi)/(pop_child*(tprime/365)))*100000
overalldisadult2_hi = (sum(disease_adult2_hi)/(pop_adult*(tprime/365)))*100000
overalldiselder2_hi = (sum(disease_elder2_hi)/(pop_elder*(tprime/365)))*100000

overallcarrprev2 = len(np.where(carriers.susc < max_carr_cap)[0])/indivs
overallcarrprevchild2 = len(np.where(carriers[carriers.agegroup == 0].susc < max_carr_cap)[0])/pop_child
overallcarrprevadult2 = len(np.where(carriers[carriers.agegroup == 1].susc < max_carr_cap)[0])/pop_adult
overallcarrprevelder2 = len(np.where(carriers[carriers.agegroup == 2].susc < max_carr_cap)[0])/pop_elder

agegrpstats2 = np.array([[overallcarrprev2, overallcarrprevchild2, overallcarrprevadult2, overallcarrprevelder2],
                        [overalldis2, overalldischild2, overalldisadult2, overalldiselder2],
                        [overalldis2_lo, overalldischild2_lo, overalldisadult2_lo, overalldiselder2_lo],
                        [overalldis2_hi, overalldischild2_hi, overalldisadult2_hi, overalldiselder2_hi]])

agegrpstats2 = pd.DataFrame(agegrpstats2, columns = ['overall','child', 'adult', 'elder'], 
                           index = ['carrprev', 'disease inc', 'disinc_lo', 'disinc_hi'])

newfiledest = filedest + endstring + "agegrpstats" + ".csv"
agegrpstats2.to_csv(path_or_buf= newfiledest, index = True)

#################################################################################################################################################################
############# PART 3: POST-VACCINATION PERIOD 2 - PCV20 ##############################################################################################
#################################################################################################################################################################

vt = np.concatenate((vt, new_vt))
vt_ind = np.concatenate((vt_ind, new_vt_ind))

vacc2_stop_t = vacc_stop_t + 3650 # current time + 10 years

# run IBM simulation
try:
    postPCV20 = pneumo_IBM_fun.pneumoIBMrun(t0 = vacc_stop_t + 1, stop_t = vacc2_stop_t, carriers = carriers, 
                                            ind_lambdas = ind_lambdas, prev_col = prev_col, vt = vt)
except Exception: 
    error_message = traceback.print_exc()
    print('simulation '+ str(simnum))

carriers = postPCV20['carriers']
prev_col = postPCV20['prev_col']
ind_lambdas = postPCV20['ind_lambdas']
plotsy3 = postPCV20['plotsy']
sim_output3 = postPCV20['sim_output']

plotsy_melt3= pd.melt(plotsy3, id_vars = 'time')
plot3 = plotsy3.plot(kind = 'line')
#plot3 = plotsy_melt3.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

phase = 'postPCV20'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

#plotsy_melt3dest = filedest + endstring + 'carrdt.csv'
#plotsy_melt3.to_csv(path_or_buf= plotsy_melt3dest, index = True)

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

disease_cas_lo3 = disinc3['carrprev']*disinc3['invasiveness.low']*indivs*(tprime/365)
disease_cas_hi3 = disinc3['carrprev']*disinc3['invasiveness.high']*indivs*(tprime/365)
disinc3['disease_cas_lo'] = disease_cas_lo3
disinc3['disease_cas_hi'] = disease_cas_hi3
disinc3['disease_inc_lo'] = (disease_cas_lo3/indivs)*100000 #dis per 100,000 ppl
disinc3['disease_inc_hi'] = (disease_cas_hi3/indivs)*100000 #dis per 100,000 ppl

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

disease_child3_lo = disinc3['carrprevchild']*disinc3['invasiveness.low']*(tprime/365)*pop_child
disease_adult3_lo = disinc3['carrprevadult']*disinc3['invasiveness.low']*(tprime/365)*pop_adult
disease_elder3_lo = disinc3['carrprevelder']*disinc3['invasiveness.low']*(tprime/365)*pop_elder
disease_child3_hi = disinc3['carrprevchild']*disinc3['invasiveness.high']*(tprime/365)*pop_child
disease_adult3_hi = disinc3['carrprevadult']*disinc3['invasiveness.high']*(tprime/365)*pop_adult
disease_elder3_hi = disinc3['carrprevelder']*disinc3['invasiveness.high']*(tprime/365)*pop_elder

# overall IPD incidence by age group
overalldis3 = ((disease_child3 + disease_adult3 + disease_elder3).sum()/(indivs*(tprime/365)))*100000
overalldischild3 = (sum(disease_child3)/(pop_child*(tprime/365)))*100000
overalldisadult3 = (sum(disease_adult3)/(pop_adult*(tprime/365)))*100000
overalldiselder3 = (sum(disease_elder3)/(pop_elder*(tprime/365)))*100000

overalldis3_lo = ((disease_child3_lo + disease_adult3_lo + disease_elder3_lo).sum()/(indivs*(tprime/365)))*100000
overalldischild3_lo = (sum(disease_child3_lo)/(pop_child*(tprime/365)))*100000
overalldisadult3_lo = (sum(disease_adult3_lo)/(pop_adult*(tprime/365)))*100000
overalldiselder3_lo = (sum(disease_elder3_lo)/(pop_elder*(tprime/365)))*100000
overalldis3_hi = ((disease_child3_hi + disease_adult3_hi + disease_elder3_hi).sum()/(indivs*(tprime/365)))*100000
overalldischild3_hi = (sum(disease_child3_hi)/(pop_child*(tprime/365)))*100000
overalldisadult3_hi = (sum(disease_adult3_hi)/(pop_adult*(tprime/365)))*100000
overalldiselder3_hi = (sum(disease_elder3_hi)/(pop_elder*(tprime/365)))*100000

overallcarrprev3 = len(np.where(carriers.susc < max_carr_cap)[0])/indivs
overallcarrprevchild3 = len(np.where(carriers[carriers.agegroup == 0].susc < max_carr_cap)[0])/pop_child
overallcarrprevadult3 = len(np.where(carriers[carriers.agegroup == 1].susc < max_carr_cap)[0])/pop_adult
overallcarrprevelder3 = len(np.where(carriers[carriers.agegroup == 2].susc < max_carr_cap)[0])/pop_elder

agegrpstats3 = np.array([[overallcarrprev3, overallcarrprevchild3, overallcarrprevadult3, overallcarrprevelder3],
                        [overalldis3, overalldischild3, overalldisadult3, overalldiselder3],
                        [overalldis3_lo, overalldischild3_lo, overalldisadult3_lo, overalldiselder3_lo],
                        [overalldis3_hi, overalldischild3_hi, overalldisadult3_hi, overalldiselder3_hi]])

agegrpstats3 = pd.DataFrame(agegrpstats3, columns = ['overall','child', 'adult', 'elder'], 
                           index = ['carrprev', 'disease inc', 'disinc_lo', 'disinc_hi'])

newfiledest3 = filedest + endstring + "agegrpstats" + ".csv"
agegrpstats3.to_csv(path_or_buf= newfiledest3, index = True)

#################################################################################################################################################################
############# PART 4: POST-VACCINATION PERIOD PCV30 ##############################################################################################
#################################################################################################################################################################

vt = np.concatenate((vt, new_vt2))
vt_ind = np.concatenate((vt_ind, new_vt_ind2))

vacc3_stop_t = vacc2_stop_t + 3650 # current time + 10 years

# run IBM simulation
try:
    postPCV30 = pneumo_IBM_fun.pneumoIBMrun(t0 = vacc2_stop_t + 1, stop_t = vacc3_stop_t, carriers = carriers, 
                                            ind_lambdas = ind_lambdas, prev_col = prev_col, vt = vt)
except Exception: 
    error_message = traceback.print_exc()
    print('simulation '+ str(simnum))

carriers = postPCV30['carriers']
prev_col = postPCV30['prev_col']
ind_lambdas = postPCV30['ind_lambdas']
plotsy4 = postPCV30['plotsy']
sim_output4 = postPCV30['sim_output']

plotsy_melt4= pd.melt(plotsy4, id_vars = 'time')
plot4 = plotsy4.plot(kind = 'line')
#plot4 = plotsy_melt4.pivot_table('value', 'time', 'variable', aggfunc='mean').plot(kind='line')

phase = 'postPCV30'
endstring = "simulation" + str(simnum) + "_" + phase + "_" + country + "_"

#plotsy_melt4dest = filedest + endstring + 'carrdt.csv'
#plotsy_melt4.to_csv(path_or_buf= plotsy_meltdest, index = True)

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

disease_cas_lo4 = disinc4['carrprev']*disinc4['invasiveness.low']*indivs*(tprime/365)
disease_cas_hi4 = disinc4['carrprev']*disinc4['invasiveness.high']*indivs*(tprime/365)
disinc4['disease_cas_lo'] = disease_cas_lo4
disinc4['disease_cas_hi'] = disease_cas_hi4
disinc4['disease_inc_lo'] = (disease_cas_lo4/indivs)*100000 #dis per 100,000 ppl
disinc4['disease_inc_hi'] = (disease_cas_hi4/indivs)*100000 #dis per 100,000 ppl

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

disease_child4_lo = disinc4['carrprevchild']*disinc4['invasiveness.low']*(tprime/365)*pop_child
disease_adult4_lo = disinc4['carrprevadult']*disinc4['invasiveness.low']*(tprime/365)*pop_adult
disease_elder4_lo = disinc4['carrprevelder']*disinc4['invasiveness.low']*(tprime/365)*pop_elder
disease_child4_hi = disinc4['carrprevchild']*disinc4['invasiveness.high']*(tprime/365)*pop_child
disease_adult4_hi = disinc4['carrprevadult']*disinc4['invasiveness.high']*(tprime/365)*pop_adult
disease_elder4_hi = disinc4['carrprevelder']*disinc4['invasiveness.high']*(tprime/365)*pop_elder

# overall IPD incidence by age group
overalldis4 =  ((disease_child4 + disease_adult4 + disease_elder4).sum()/(indivs*(tprime/365)))*100000
overalldischild4 = (sum(disease_child4)/(pop_child*(tprime/365)))*100000
overalldisadult4 = (sum(disease_adult4)/(pop_adult*(tprime/365)))*100000
overalldiselder4 = (sum(disease_elder4)/(pop_elder*(tprime/365)))*100000

overalldis4_lo = ((disease_child4_lo + disease_adult4_lo + disease_elder4_lo).sum()/(indivs*(tprime/365)))*100000
overalldischild4_lo = (sum(disease_child4_lo)/(pop_child*(tprime/365)))*100000
overalldisadult4_lo = (sum(disease_adult4_lo)/(pop_adult*(tprime/365)))*100000
overalldiselder4_lo = (sum(disease_elder4_lo)/(pop_elder*(tprime/365)))*100000
overalldis4_hi = ((disease_child4_hi + disease_adult4_hi + disease_elder4_hi).sum()/(indivs*(tprime/365)))*100000
overalldischild4_hi = (sum(disease_child4_hi)/(pop_child*(tprime/365)))*100000
overalldisadult4_hi = (sum(disease_adult4_hi)/(pop_adult*(tprime/365)))*100000
overalldiselder4_hi = (sum(disease_elder4_hi)/(pop_elder*(tprime/365)))*100000

overallcarrprev4 = len(np.where(carriers.susc < max_carr_cap)[0])/indivs
overallcarrprevchild4 = len(np.where(carriers[carriers.agegroup == 0].susc < max_carr_cap)[0])/pop_child
overallcarrprevadult4 = len(np.where(carriers[carriers.agegroup == 1].susc < max_carr_cap)[0])/pop_adult
overallcarrprevelder4 = len(np.where(carriers[carriers.agegroup == 2].susc < max_carr_cap)[0])/pop_elder

agegrpstats4 = np.array([[overallcarrprev4, overallcarrprevchild4, overallcarrprevadult4, overallcarrprevelder4],
                        [overalldis4, overalldischild4, overalldisadult4, overalldiselder4],
                        [overalldis4_lo, overalldischild4_lo, overalldisadult4_lo, overalldiselder4_lo],
                        [overalldis4_hi, overalldischild4_hi, overalldisadult4_hi, overalldiselder4_hi]])

agegrpstats4 = pd.DataFrame(agegrpstats4, columns = ['overall','child', 'adult', 'elder'], 
                           index = ['carrprev', 'disease inc', 'disinc_lo', 'disinc_hi'])

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

# stdev
stdev_post13_child = math.sqrt((1/sum(disease_child)) + (1/sum(disease_child2)))
stdev_post20_child = math.sqrt((1/sum(disease_child3)) + (1/sum(disease_child2)))
stdev_post30_child = math.sqrt((1/sum(disease_child3)) + (1/sum(disease_child4)))
stdev_all_child = math.sqrt((1/sum(disease_child)) + (1/sum(disease_child4)))

# bounds 95% CI
post13_child_bounds = np.array([math.exp(math.log(post13_child)-(1.96*stdev_post13_child)),
                                math.exp(math.log(post13_child)+(1.96*stdev_post13_child))])
post20_child_bounds = np.array([math.exp(math.log(post20_child)-(1.96*stdev_post20_child)),
                                math.exp(math.log(post20_child)+(1.96*stdev_post20_child))])
post30_child_bounds = np.array([math.exp(math.log(post30_child)-(1.96*stdev_post30_child)),
                                math.exp(math.log(post30_child)+(1.96*stdev_post30_child))])
all_child_bounds = np.array([math.exp(math.log(all_child)-(1.96*stdev_all_child)),
                            math.exp(math.log(all_child)+(1.96*stdev_all_child))])

IRR_bounds_child = np.array((post13_child_bounds, post20_child_bounds, post30_child_bounds, all_child_bounds))
IRR_children = pd.DataFrame(IRR_bounds_child, columns = ['low', 'high'], 
             index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])
IRR_children['IRR'] = IRR_child
IRR_children['agegrp'] = 'children'

# adults
post13_adult = overalldisadult2/overalldisadult
post20_adult = overalldisadult3/overalldisadult2
post30_adult = overalldisadult4/overalldisadult3
all_adult = overalldisadult4/overalldisadult

IRR_adult = np.array((post13_adult, post20_adult, post30_adult, all_adult))

# stdev
stdev_post13_adult = math.sqrt((1/sum(disease_adult)) + (1/sum(disease_adult2)))
stdev_post20_adult = math.sqrt((1/sum(disease_adult3)) + (1/sum(disease_adult2)))
stdev_post30_adult = math.sqrt((1/sum(disease_adult3)) + (1/sum(disease_adult4)))
stdev_all_adult = math.sqrt((1/sum(disease_adult)) + (1/sum(disease_adult4)))

# bounds 95% CI
post13_adult_bounds = np.array([math.exp(math.log(post13_adult)-(1.96*stdev_post13_adult)),
                                math.exp(math.log(post13_adult)+(1.96*stdev_post13_adult))])
post20_adult_bounds = np.array([math.exp(math.log(post20_adult)-(1.96*stdev_post20_adult)),
                                math.exp(math.log(post20_adult)+(1.96*stdev_post20_adult))])
post30_adult_bounds = np.array([math.exp(math.log(post30_adult)-(1.96*stdev_post30_adult)),
                                math.exp(math.log(post30_adult)+(1.96*stdev_post30_adult))])
all_adult_bounds = np.array([math.exp(math.log(all_adult)-(1.96*stdev_all_adult)),
                                math.exp(math.log(all_adult)+(1.96*stdev_all_adult))])

IRR_bounds_adult = np.array((post13_adult_bounds, post20_adult_bounds, post30_adult_bounds, all_adult_bounds))
IRR_adults = pd.DataFrame(IRR_bounds_adult, columns = ['low', 'high'], 
             index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])
IRR_adults['IRR'] = IRR_adult
IRR_adults['agegrp'] = 'adults'

# elderly
post13_elder = overalldiselder2/overalldiselder
post20_elder = overalldiselder3/overalldiselder2
post30_elder = overalldiselder4/overalldiselder3
all_elder = overalldiselder4/overalldiselder

IRR_elderly = np.array((post13_elder, post20_elder, post30_elder, all_elder))

# stdev
stdev_post13_elder = math.sqrt((1/sum(disease_elder)) + (1/sum(disease_elder2)))
stdev_post20_elder = math.sqrt((1/sum(disease_elder3)) + (1/sum(disease_elder2)))
stdev_post30_elder = math.sqrt((1/sum(disease_elder3)) + (1/sum(disease_elder4)))
stdev_all_elder = math.sqrt((1/sum(disease_elder)) + (1/sum(disease_elder4)))

# bounds 95% CI
post13_elder_bounds = np.array([math.exp(math.log(post13_elder)-(1.96*stdev_post13_elder)),
                                math.exp(math.log(post13_elder)+(1.96*stdev_post13_elder))])
post20_elder_bounds = np.array([math.exp(math.log(post20_elder)-(1.96*stdev_post20_elder)),
                                math.exp(math.log(post20_elder)+(1.96*stdev_post20_elder))])
post30_elder_bounds = np.array([math.exp(math.log(post30_elder)-(1.96*stdev_post30_elder)),
                                math.exp(math.log(post30_elder)+(1.96*stdev_post30_elder))])
all_elder_bounds = np.array([math.exp(math.log(all_elder)-(1.96*stdev_all_elder)),
                                math.exp(math.log(all_elder)+(1.96*stdev_all_elder))])

IRR_bounds_elder = np.array((post13_elder_bounds, post20_elder_bounds, post30_elder_bounds, all_elder_bounds))
IRR_elders = pd.DataFrame(IRR_bounds_elder, columns = ['low', 'high'], 
             index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])
IRR_elders['IRR'] = IRR_elderly
IRR_elders['agegrp'] = 'elders'

# overall
post13_all = overalldis2/overalldis
post20_all = overalldis3/overalldis2
post30_all = overalldis4/overalldis3
all_all = overalldis4/overalldis

IRR_all = np.array((post13_all, post20_all, post30_all, all_all))

stdev_post13_all = math.sqrt((1/(disease_child + disease_adult + disease_elder).sum()) + (1/(disease_child2 + disease_adult2 + disease_elder2).sum()))
stdev_post20_all = math.sqrt((1/(disease_child3 + disease_adult3 + disease_elder3).sum()) + (1/(disease_child2 + disease_adult2 + disease_elder2).sum()))
stdev_post30_all = math.sqrt((1/(disease_child3 + disease_adult3 + disease_elder3).sum()) + (1/(disease_child4 + disease_adult4 + disease_elder4).sum()))
stdev_all_all = math.sqrt((1/(disease_child + disease_adult + disease_elder).sum()) + (1/(disease_child4 + disease_adult4 + disease_elder4).sum()))

# bounds 95% CI
post13_all_bounds = np.array([math.exp(math.log(post13_all)-(1.96*stdev_post13_all)),
                                math.exp(math.log(post13_all)+(1.96*stdev_post13_all))])
post20_all_bounds = np.array([math.exp(math.log(post20_all)-(1.96*stdev_post20_all)),
                                math.exp(math.log(post20_all)+(1.96*stdev_post20_all))])
post30_all_bounds = np.array([math.exp(math.log(post30_all)-(1.96*stdev_post30_all)),
                                math.exp(math.log(post30_all)+(1.96*stdev_post30_all))])
all_all_bounds = np.array([math.exp(math.log(all_all)-(1.96*stdev_all_all)),
                                math.exp(math.log(all_all)+(1.96*stdev_all_all))])

IRR_bounds_all = np.array((post13_all_bounds, post20_all_bounds, post30_all_bounds, all_all_bounds))
IRR_tot = pd.DataFrame(IRR_bounds_all, columns = ['low', 'high'], 
             index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])
IRR_tot['IRR'] = IRR_all
IRR_tot['agegrp'] = 'all'


### Make IRR full df
IRR_df = pd.concat([IRR_children, IRR_adults, IRR_elders, IRR_tot])
# IRR_df = pd.DataFrame({'children': IRR_child, 'adults': IRR_adult, 'elder': IRR_elder, 'overall': IRR_all}, 
#                       index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])

newfiledest5 = filedest + endstring + "IRR_allperiods" + ".csv" 
IRR_df.to_csv(path_or_buf = newfiledest5, index = True) # incidence rate ratio dataframe


### INCIDENCE RATES

pre = agegrpstats.loc['disease inc']
pre_lo = agegrpstats.loc['disinc_lo']
pre_hi = agegrpstats.loc['disinc_hi']
post13 = agegrpstats2.loc['disease inc']
post13_lo = agegrpstats2.loc['disinc_lo']
post13_hi = agegrpstats2.loc['disinc_hi']
post20 = agegrpstats3.loc['disease inc']
post20_lo = agegrpstats3.loc['disinc_lo']
post20_hi = agegrpstats3.loc['disinc_hi']
post30 = agegrpstats4.loc['disease inc']
post30_lo = agegrpstats4.loc['disinc_lo']
post30_hi = agegrpstats4.loc['disinc_hi']

agespec = pd.concat([pre, pre_lo, pre_hi, post13, post13_lo, post13_hi, post20, post20_lo, post20_hi, 
                     post30, post30_lo, post30_hi], axis=1)
agespec.columns = ['pre', 'pre_lo', 'pre_hi', 'post13', 'post13_lo', 'post13_hi', 'post20', 'post20_lo', 'post20_hi', 
                   'post30', 'post30_lo', 'post30_hi'] 

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