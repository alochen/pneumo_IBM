# pneumo_IBM
IBM for understanding the impact of increasing serotype valency in PCVs

File descriptions:

*Data import*

-fralocalinv.csv: Serotype-specific invasiveness in France post-PCV7, from https://github.com/alochen/pneumo-invasiveness
-usalocalinv.csv: Serotype-specific invasiveness in Massachussetts post-PCV7, from https://github.com/alochen/pneumo-invasiveness
-contactmatrices_balanced.xlsx: Balanced contact matrices for France and USA where France based on contacts from Beraud study https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0133203 and USA based on POLYMOD study in UK by Mossong et al https://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.0050074

*Parameter-fitting*

-IBM_parameterfitting.py: Python file that takes in carriage_episodes_AL.csv file (Maela refugee camp carriage data file) and fits recovery rate parameters using emcee package in Python; also runs naive population using these parameters for previous colonisations and shows qualitative comparison between observed and fitted distribution of previous colonisation.
-seroparams.csv: Output of IBM_parameterfitting.py that has the fitted parameters for the recovery rate distribution based on the previous colonisations.

*Simulations*

-simparams.csv: Input parameter file where each row represents the parameter set for a specific scenario
-pneumo_shared.py: Initialisation of population, previous colonisation distributions, rate assignment
-pneumo_IBM_fun.py: Simulation function file; includes functions to compute an individual's rates after they have been chosen for a switching event and function that runs IBM
-pneumo_IBM_localinv.py: File that runs pneumo_shared.py and pneumo_IBM_fun.py for the local invasiveness vaccine strategy. This file runs the pre-vaccination period which is shared by the disease incidence and global invasiveness vaccine strategies. Outputs files on the total number of carriers for each serotype at the end of each period, number of previous colonisations of each person, carriage prevalence of each serotype. Computes disease incidence from invasiveness files and the IRRs over the vaccination periods.
-pneumo_IBM_disinc.py: File that runs pneumo_shared.py and pneumo_IBM_fun.py for the disease incidence vaccine strategy. This file uses the pre-vaccination period from the corresponding local invasiveness scenario run. Outputs files on the total number of carriers for each serotype at the end of each period, number of previous colonisations of each person, carriage prevalence of each serotype. Computes disease incidence from invasiveness files and the IRRs over the vaccination periods.
-pneumo_IBM_globinv.py: File that runs pneumo_shared.py and pneumo_IBM_fun.py for the global invasiveness vaccine strategy. This file uses the pre-vaccination period from the corresponding local invasiveness scenario run. Outputs files on the total number of carriers for each serotype at the end of each period, number of previous colonisations of each person, carriage prevalence of each serotype. Computes disease incidence from invasiveness files and the IRRs over the vaccination periods.

*Post-simulation analysis*

-IBM_analysis.R: Post-simulation analysis that fixes IPD incidence rates and estimates fixed effects size of change in IPD incidence before and after vaccination

