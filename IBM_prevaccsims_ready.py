##########################################################
##### GET INVASIVE DATASET - COMMENTED OUT FOR FILE IMPORT, this only needs to be done one time
##########################################################

import pandas as pd
import numpy as np
import os
from os import path
from glob import glob

frainv_fol = "//qdrive/homes/al817/Technical/Python/popsims/fra_inv/"
usainv_fol = "//qdrive/homes/al817/Technical/Python/popsims/usa_inv/"
fraglobinv_fol = "//qdrive/homes/al817/Technical/Python/popsims/fra_globinv/"
usaglobinv_fol = "//qdrive/homes/al817/Technical/Python/popsims/usa_globinv/"
fradis_fol = "//qdrive/homes/al817/Technical/Python/popsims/fra_dis/"
usadis_fol = "//qdrive/homes/al817/Technical/Python/popsims/usa_dis/"


def find_ext(dr, ext):
    return glob(path.join(dr,"*.{}".format(ext)))

os.chdir(frainv_fol)

### Which fra locinv have finished running completely
pdflist = find_ext(".", "pdf")
newlist = [x.split('simulation')[1] for x in pdflist]
newlist2 = [int(x.split('_')[0]) for x in newlist]
newlist2.sort()
finito = np.unique(newlist2)

### Which fra locinv have not finished running completely (failed)
finishedlist = find_ext(".", "csv")
finishedmatch = [s for s in finishedlist if "postPCV30" in s]
finlist2 = [x.split('simulation')[1] for x in finishedmatch]
finlist2 = [int(x.split('_')[0]) for x in finlist2]
finlist2.sort()
finished = np.unique(finlist2)

os.chdir(usainv_fol)

### Which usa locinv have finished running completely
pdflist = find_ext(".", "pdf")
newlist = [x.split('simulation')[1] for x in pdflist]
newlist2 = [int(x.split('_')[0]) for x in newlist]
newlist2.sort()
np.unique(newlist2)

### Which usa locinv have not finished running completely (failed)
finishedlist_usa = find_ext(".", "csv")
finishedmatch_usa = [s for s in finishedlist_usa if "postPCV30" in s]
finlist2_usa = [x.split('simulation')[1] for x in finishedmatch_usa]
finlist2_usa = [int(x.split('_')[0]) for x in finlist2_usa]
finlist2_usa.sort()
finished_usa = np.unique(finlist2_usa)

bothfinished = [*finished, *finished_usa]
unfinished_full = list(range(1,444))
diff = np.setdiff1d(unfinished_full, bothfinished) # which simulations are done running but failed (and need to be re-run)

### Which disinc and globinv sims need to be re-run because they didn't run properly on the cluster?
os.chdir(fradis_fol)
import time
from datetime import datetime

ext = 'csv' 
start = '03/08/2021'
end = '25/08/2021'

def dateRange(createdDate, startDate, endDate):
    """determines if date is in range"""
    createdDate = datetime.strptime(createdDate, '%a %b %d %H:%M:%S %Y')
    startDate = datetime.strptime(startDate, '%d/%m/%Y')
    endDate = datetime.strptime(endDate, '%d/%m/%Y')
    return startDate < createdDate < endDate

# fra dis
src = fradis_fol
r = []
for filename in os.listdir(src):
    created = time.ctime(os.path.getmtime(src + filename))
    if filename.endswith('.' + ext) and dateRange(created, start, end):
        r.append(filename)

### which sims are finished?
rfinlist2 = [x.split('simulation')[1] for x in r]
# rfinlist2 = [s for s in r if "postPCV30" in s]
# rfinlist2 = [x.split('simulation')[1] for x in rfinlist2]

rfinlist2 = [int(x.split('_')[0]) for x in rfinlist2]
rfinlist2.sort()
rrunning = np.unique(rfinlist2)

# usa dis
src = usadis_fol
q = []
for filename in os.listdir(src):
    created = time.ctime(os.path.getmtime(src + filename))
    if filename.endswith('.' + ext) and dateRange(created, start, end):
        q.append(filename)


qfinlist2 = [x.split('simulation')[1] for x in q]
# qfinlist2 = [s for s in q if "postPCV30" in s]
# qfinlist2 = [x.split('simulation')[1] for x in qfinlist2]

qfinlist2 = [int(x.split('_')[0]) for x in qfinlist2]
qfinlist2.sort()
qrunning = np.unique(qfinlist2)

bothrunning = [*rrunning, *qrunning]
running_full = list(range(1,444))
diff_dis_run = np.setdiff1d(running_full, bothrunning) # which simulations are done but failed (and need to be re-run)


# fra globinv
src = fraglobinv_fol
r = []
for filename in os.listdir(src):
    created = time.ctime(os.path.getmtime(src + filename))
    if filename.endswith('.' + ext) and dateRange(created, start, end):
        r.append(filename)

rfinlist2 = [x.split('simulation')[1] for x in r]
rfinlist2 = [int(x.split('_')[0]) for x in rfinlist2]
rfinlist2.sort()
rrunning = np.unique(rfinlist2)

# usa globinv
src = usaglobinv_fol
q = []
for filename in os.listdir(src):
    created = time.ctime(os.path.getmtime(src + filename))
    if filename.endswith('.' + ext) and dateRange(created, start, end):
        q.append(filename)

qfinlist2 = [x.split('simulation')[1] for x in q]
qfinlist2 = [int(x.split('_')[0]) for x in qfinlist2]
qfinlist2.sort()
qrunning = np.unique(qfinlist2)

gl_bothrunning = [*rrunning, *qrunning]
gl_running_full = list(range(1,444))
gl_diff_dis_run = np.setdiff1d(gl_running_full, gl_bothrunning) # which simulations are done but failed (and need to be re-run)

