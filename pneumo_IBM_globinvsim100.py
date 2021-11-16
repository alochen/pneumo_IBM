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

country = sh.country
indivs = sh.indivs
globalinv = True
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

filedest = "//qdrive/homes/al817/Technical/Python/popsims/" + country + "_globinv/" 

phase = 'prevacc'
# mylist = [f for f in glob.glob(filedest + "*.csv")]
# mylist = [fname.rsplit('.', 1)[0] for fname in os.listdir(filedest)]
# filenums = [item.split('_',1)[0] for item in mylist]
# filenums = [item.split('simulation',1)[-1] for item in filenums]
# filenums = [int(i) for i in filenums] 
# simnum = max(filenums) + 1
string1 = str(simnum) + '_' + globalinv*("glob") + "inv" + "pop" + str(indivs)
endstring = "simulation" + string1 + "_" + phase + "_" + country + "_" #+ 'initFOIbeta_' #+ 'vaccdur15_'
#seroparams.to_csv(path_or_buf = filedest + endstring + 'dummy.csv', index = False)

# log file of main params
paramlogfile = dict(simnum = simnum, country = country, beta = beta, indivs = indivs, vacc_dur = vacc_dur, vacc_eff = vacc_eff,
                    sigma = sigma, epsilon = epsilon, theta = theta)
paramdf = pd.DataFrame(dict([(k,pd.Series(v)) for k,v in paramlogfile.items() ])) 
filedestparams = filedest + endstring + 'paramlogfile.csv'
paramdf.to_csv(path_or_buf = filedestparams, index = False)

# PCV13
# vt = np.array(country_inv.nlargest(13, 'invasiveness').Serotype)
# vt_ind = np.where(np.isin(colnames, vt))[0]
# # PCV20
# PCV20df = country_inv[country_inv['Serotype'].isin(vt) == False]
# new_vt = np.array(PCV20df.nlargest(7, 'invasiveness').Serotype)
# new_vt_ind = np.where(np.isin(colnames, new_vt))[0]
# # PCV30
# PCV30df = PCV20df[PCV20df['Serotype'].isin(new_vt) == False]
# new_vt2 = np.array(PCV30df.nlargest(10, 'invasiveness').Serotype)
# new_vt_ind2 = np.where(np.isin(colnames, new_vt2))[0]

# top 13 globally invasive serotypes get vaccinated against
# if globalinv == True:
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
############# PART 1: POST-VACCINATION PERIOD PCV30 #################################################################################################################
#################################################################################################################################################################
vt = np.concatenate((vt, new_vt))
vt_ind = np.concatenate((vt_ind, new_vt_ind))
vt = np.concatenate((vt, new_vt2))
vt_ind = np.concatenate((vt_ind, new_vt_ind2))

stop_t = 14600
postPCV20folder = "//qdrive/homes/al817/Technical/Python/popsims/" + country + "_inv/"
mylist2 = [f for f in glob.glob(postPCV20folder + "simulation" + str(simnum) + "_*")]
post20_carr_file = [s for s in mylist2 if "postPCV20_" + country + "_carrierstmax" in s]
carriers = pd.read_csv(post20_carr_file[0], index_col=0)
post20_prevcol_file = [s for s in mylist2 if "postPCV20_" + country + "_prevcol" in s]
prev_col = pd.read_csv(post20_prevcol_file[0], index_col=0)
ind_lambdas = carriers.iloc[:,0:num_sero]
disincfile = [s for s in mylist2 if "postPCV20_" + country + "_disinc" in s]
disinc = pd.read_csv(disincfile[0], index_col=0)
agegrpstatsfile = [s for s in mylist2 if "postPCV20_" + country + "_agegrpstats" in s]
agegrpstats = pd.read_csv(agegrpstatsfile[0], index_col=0)
overalldis = agegrpstats.loc['disease inc'][0]
overalldischild = agegrpstats.loc['disease inc'][1]
overalldisadult = agegrpstats.loc['disease inc'][2]
overalldiselder = agegrpstats.loc['disease inc'][3]

# index of children and sum of all children, adult and elderly carriers
index_children = np.where(carriers.age < age_breakdown[0])[0]
index_adults = np.where((carriers.age >= age_breakdown[0]) & (carriers.age < age_breakdown[1]))[0]
index_elderly = np.where(carriers.age >= age_breakdown[1])[0]

pop_child = len(index_children)
pop_adult = len(index_adults)
pop_elder = len(index_elderly)

# disease cases by serotype
disease_child = disinc['carrprevchild']*disinc['invasiveness']*(stop_t/365)*pop_child
disease_adult = disinc['carrprevadult']*disinc['invasiveness']*(stop_t/365)*pop_adult
disease_elder = disinc['carrprevelder']*disinc['invasiveness']*(stop_t/365)*pop_elder

totchildcarr = carriers.iloc[index_children,0:num_sero].sum(axis = 0)
totadultcarr = carriers.iloc[index_adults,0:num_sero].sum(axis = 0)
toteldercarr = carriers.iloc[index_elderly,0:num_sero].sum(axis = 0)

# calculation of initial rates of individuals
totFOIbysero_children = [0]*num_sero
totFOIbysero_adults = [0]*num_sero
totFOIbysero_elders = [0]*num_sero
    
for i in range(num_sero):
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

vacc_stop_t = stop_t + 3650#1825 # current time + 5 years

# run IBM simulation
try:
    postPCV30 = pneumo_IBM_fun.pneumoIBMrun(t0 = stop_t + 1, stop_t = vacc_stop_t, carriers = carriers, 
                                            ind_lambdas = ind_lambdas, prev_col = prev_col, vt = vt)
except Exception: 
    error_message = traceback.print_exc()
    print('globinv_simulation '+ str(simnum))

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
tprime = 3650 #vacc3_stop_t - vacc2_stop_t
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

# # children
# post13_child = overalldischild2/overalldischild
# post20_child = overalldischild3/overalldischild2
# post30_child = overalldischild4/overalldischild3
# all_child = overalldischild4/overalldischild
# 
# IRR_child = np.array((post13_child, post20_child, post30_child, all_child))
# 
# # stdev
# stdev_post13_child = math.sqrt((1/sum(disease_child)) + (1/sum(disease_child2)))
# stdev_post20_child = math.sqrt((1/sum(disease_child3)) + (1/sum(disease_child2)))
# stdev_post30_child = math.sqrt((1/sum(disease_child3)) + (1/sum(disease_child4)))
# stdev_all_child = math.sqrt((1/sum(disease_child)) + (1/sum(disease_child4)))
# 
# # bounds 95% CI
# post13_child_bounds = np.array([math.exp(math.log(post13_child)-(1.96*stdev_post13_child)),
#                                 math.exp(math.log(post13_child)+(1.96*stdev_post13_child))])
# post20_child_bounds = np.array([math.exp(math.log(post20_child)-(1.96*stdev_post20_child)),
#                                 math.exp(math.log(post20_child)+(1.96*stdev_post20_child))])
# post30_child_bounds = np.array([math.exp(math.log(post30_child)-(1.96*stdev_post30_child)),
#                                 math.exp(math.log(post30_child)+(1.96*stdev_post30_child))])
# all_child_bounds = np.array([math.exp(math.log(all_child)-(1.96*stdev_all_child)),
#                             math.exp(math.log(all_child)+(1.96*stdev_all_child))])
# 
# IRR_bounds_child = np.array((post13_child_bounds, post20_child_bounds, post30_child_bounds, all_child_bounds))
# IRR_children = pd.DataFrame(IRR_bounds_child, columns = ['low', 'high'], 
#              index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])
# IRR_children['IRR'] = IRR_child
# IRR_children['agegrp'] = 'children'
# 
# # adults
# post13_adult = overalldisadult2/overalldisadult
# post20_adult = overalldisadult3/overalldisadult2
# post30_adult = overalldisadult4/overalldisadult3
# all_adult = overalldisadult4/overalldisadult
# 
# IRR_adult = np.array((post13_adult, post20_adult, post30_adult, all_adult))
# 
# # stdev
# stdev_post13_adult = math.sqrt((1/sum(disease_adult)) + (1/sum(disease_adult2)))
# stdev_post20_adult = math.sqrt((1/sum(disease_adult3)) + (1/sum(disease_adult2)))
# stdev_post30_adult = math.sqrt((1/sum(disease_adult3)) + (1/sum(disease_adult4)))
# stdev_all_adult = math.sqrt((1/sum(disease_adult)) + (1/sum(disease_adult4)))
# 
# # bounds 95% CI
# post13_adult_bounds = np.array([math.exp(math.log(post13_adult)-(1.96*stdev_post13_adult)),
#                                 math.exp(math.log(post13_adult)+(1.96*stdev_post13_adult))])
# post20_adult_bounds = np.array([math.exp(math.log(post20_adult)-(1.96*stdev_post20_adult)),
#                                 math.exp(math.log(post20_adult)+(1.96*stdev_post20_adult))])
# post30_adult_bounds = np.array([math.exp(math.log(post30_adult)-(1.96*stdev_post30_adult)),
#                                 math.exp(math.log(post30_adult)+(1.96*stdev_post30_adult))])
# all_adult_bounds = np.array([math.exp(math.log(all_adult)-(1.96*stdev_all_adult)),
#                                 math.exp(math.log(all_adult)+(1.96*stdev_all_adult))])
# 
# IRR_bounds_adult = np.array((post13_adult_bounds, post20_adult_bounds, post30_adult_bounds, all_adult_bounds))
# IRR_adults = pd.DataFrame(IRR_bounds_adult, columns = ['low', 'high'], 
#              index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])
# IRR_adults['IRR'] = IRR_adult
# IRR_adults['agegrp'] = 'adults'
# 
# # elderly
# post13_elder = overalldiselder2/overalldiselder
# post20_elder = overalldiselder3/overalldiselder2
# post30_elder = overalldiselder4/overalldiselder3
# all_elder = overalldiselder4/overalldiselder
# 
# IRR_elderly = np.array((post13_elder, post20_elder, post30_elder, all_elder))
# 
# # stdev
# stdev_post13_elder = math.sqrt((1/sum(disease_elder)) + (1/sum(disease_elder2)))
# stdev_post20_elder = math.sqrt((1/sum(disease_elder3)) + (1/sum(disease_elder2)))
# stdev_post30_elder = math.sqrt((1/sum(disease_elder3)) + (1/sum(disease_elder4)))
# stdev_all_elder = math.sqrt((1/sum(disease_elder)) + (1/sum(disease_elder4)))
# 
# # bounds 95% CI
# post13_elder_bounds = np.array([math.exp(math.log(post13_elder)-(1.96*stdev_post13_elder)),
#                                 math.exp(math.log(post13_elder)+(1.96*stdev_post13_elder))])
# post20_elder_bounds = np.array([math.exp(math.log(post20_elder)-(1.96*stdev_post20_elder)),
#                                 math.exp(math.log(post20_elder)+(1.96*stdev_post20_elder))])
# post30_elder_bounds = np.array([math.exp(math.log(post30_elder)-(1.96*stdev_post30_elder)),
#                                 math.exp(math.log(post30_elder)+(1.96*stdev_post30_elder))])
# all_elder_bounds = np.array([math.exp(math.log(all_elder)-(1.96*stdev_all_elder)),
#                                 math.exp(math.log(all_elder)+(1.96*stdev_all_elder))])
# 
# IRR_bounds_elder = np.array((post13_elder_bounds, post20_elder_bounds, post30_elder_bounds, all_elder_bounds))
# IRR_elders = pd.DataFrame(IRR_bounds_elder, columns = ['low', 'high'], 
#              index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])
# IRR_elders['IRR'] = IRR_elderly
# IRR_elders['agegrp'] = 'elders'
# 
# # overall
# post13_all = overalldis2/overalldis
# post20_all = overalldis3/overalldis2
# post30_all = overalldis4/overalldis3
# all_all = overalldis4/overalldis
# 
# IRR_all = np.array((post13_all, post20_all, post30_all, all_all))
# 
# stdev_post13_all = math.sqrt((1/(disease_child + disease_adult + disease_elder).sum()) + (1/(disease_child2 + disease_adult2 + disease_elder2).sum()))
# stdev_post20_all = math.sqrt((1/(disease_child3 + disease_adult3 + disease_elder3).sum()) + (1/(disease_child2 + disease_adult2 + disease_elder2).sum()))
# stdev_post30_all = math.sqrt((1/(disease_child3 + disease_adult3 + disease_elder3).sum()) + (1/(disease_child4 + disease_adult4 + disease_elder4).sum()))
# stdev_all_all = math.sqrt((1/(disease_child + disease_adult + disease_elder).sum()) + (1/(disease_child4 + disease_adult4 + disease_elder4).sum()))
# 
# # bounds 95% CI
# post13_all_bounds = np.array([math.exp(math.log(post13_all)-(1.96*stdev_post13_all)),
#                                 math.exp(math.log(post13_all)+(1.96*stdev_post13_all))])
# post20_all_bounds = np.array([math.exp(math.log(post20_all)-(1.96*stdev_post20_all)),
#                                 math.exp(math.log(post20_all)+(1.96*stdev_post20_all))])
# post30_all_bounds = np.array([math.exp(math.log(post30_all)-(1.96*stdev_post30_all)),
#                                 math.exp(math.log(post30_all)+(1.96*stdev_post30_all))])
# all_all_bounds = np.array([math.exp(math.log(all_all)-(1.96*stdev_all_all)),
#                                 math.exp(math.log(all_all)+(1.96*stdev_all_all))])
# 
# IRR_bounds_all = np.array((post13_all_bounds, post20_all_bounds, post30_all_bounds, all_all_bounds))
# IRR_tot = pd.DataFrame(IRR_bounds_all, columns = ['low', 'high'], 
#              index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])
# IRR_tot['IRR'] = IRR_all
# IRR_tot['agegrp'] = 'all'
# 
# 
# ### Make IRR full df
# IRR_df = pd.concat([IRR_children, IRR_adults, IRR_elders, IRR_tot])
# # IRR_df = pd.DataFrame({'children': IRR_child, 'adults': IRR_adult, 'elder': IRR_elder, 'overall': IRR_all}, 
# #                       index = ['post13 vs pre', 'post20 vs post13', 'post30 vs post20', 'post30 vs pre'])
# 
# newfiledest5 = filedest + endstring + "IRR_allperiods" + ".csv" 
# IRR_df.to_csv(path_or_buf = newfiledest5, index = True) # incidence rate ratio dataframe


### INCIDENCE RATES
# 
# pre = agegrpstats.loc['disease inc']
# pre_lo = agegrpstats.loc['disinc_lo']
# pre_hi = agegrpstats.loc['disinc_hi']
# post13 = agegrpstats2.loc['disease inc']
# post13_lo = agegrpstats2.loc['disinc_lo']
# post13_hi = agegrpstats2.loc['disinc_hi']
# post20 = agegrpstats3.loc['disease inc']
# post20_lo = agegrpstats3.loc['disinc_lo']
# post20_hi = agegrpstats3.loc['disinc_hi']
# post30 = agegrpstats4.loc['disease inc']
# post30_lo = agegrpstats4.loc['disinc_lo']
# post30_hi = agegrpstats4.loc['disinc_hi']
# 
# agespec = pd.concat([pre, pre_lo, pre_hi, post13, post13_lo, post13_hi, post20, post20_lo, post20_hi, 
#                      post30, post30_lo, post30_hi], axis=1)
# agespec.columns = ['pre', 'pre_lo', 'pre_hi', 'post13', 'post13_lo', 'post13_hi', 'post20', 'post20_lo', 'post20_hi', 
#                    'post30', 'post30_lo', 'post30_hi'] 
# 
# newfiledest5 = filedest + endstring + "disinc_eachperiod" + ".csv" 
# agespec.to_csv(path_or_buf = newfiledest5, index = True) # incidence rate dataframe
# 
# ### CARRIAGE PREVALENCE
# 
# pre_cp = agegrpstats.loc['carrprev']
# post13_cp = agegrpstats2.loc['carrprev']
# post20_cp = agegrpstats3.loc['carrprev']
# post30_cp = agegrpstats4.loc['carrprev']
# 
# agespec_cp = pd.concat([pre_cp, post13_cp, post20_cp, post30_cp], axis=1)
# agespec_cp.columns = ['pre', 'post13', 'post20', 'post30']
# 
# newfiledest5 = filedest + endstring + "carrprev_eachperiod" + ".csv" 
# agespec_cp.to_csv(path_or_buf = newfiledest5, index = True) # carr prev dataframe

#################################################################################################################################################################

# import win32com.client as win32
# from win32 import win32api
# outlook = win32.Dispatch('outlook.application')
# mail = outlook.CreateItem(0)
# mail.To = 'a.lochen17@imperial.ac.uk'
# mail.Subject = country + 'simulation ' + str(simnum) + ' has finished running'
# mail.Body = 'Congrats!'
# mail.Send()