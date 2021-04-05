# network is randomized to test the significance of hits
import getVarGene
import random
import matplotlib.pyplot as mplt
import numpy
#import seaborn as splt

iters = 1000

pnet = getVarGene.biogridInts()
pl = list(pnet.keys())
vad = getVarGene.ebigwas()
vl = list(vad.keys())
vv = list(set(sum(list(vad.values()),[])))
vvp = list(set(sum(list(vad.values()),[])) & set(pl))

# reported vis
impvis = []
realnenos = {}
randomnenos = {}

nsizes = []
for vi in vl:
    tls = []
    tv = vad[vi]
    for ti in tv:
        if ti in pnet.keys():
            tls.extend(pnet[ti])
    tls = list(set(tls))
    tno = len(tls)
    nsizes.append(tno)
    realnenos[vi] = tno
    randomnenos[vi] = []
    if tno > 200:
        impvis.append(vi)
        #print(vi)
        #print(vad[vi])
        #print(tno)

# randomize pnet and repeat the same for randoms
sumnsizes = []
mednsizes = []
onepnsizes = []
tenpnsizes = []
meanrns = []
stdrns = []
rnall = []
for ni in nsizes:
    rnall.append([])

for i in range(iters):
    print(i)
    rpnet = getVarGene.netShuffle(pnet)
    rnsizes = []
    for vi in vl:
        tls = []
        tv = vad[vi]
        for ti in tv:
            if ti in rpnet.keys():
                tls.extend(rpnet[ti])
        tls = list(set(tls))
        tno = len(tls)
        rnsizes.append(tno)
        randomnenos[vi].append(tno)
            
    rnsizes.sort()
    sumnsizes.append(sum(rnsizes))
    mednsizes.append(numpy.median(rnsizes))
    onepnsizes.append(rnsizes[len(rnsizes)-1-int(len(rnsizes)/100)])
    tenpnsizes.append(rnsizes[len(rnsizes)-1-int(len(rnsizes)/10)])
    for rr in range(len(rnsizes)):
        rnall[rr].append(rnsizes[rr])

for rr in rnall:
    meanrns.append(numpy.mean(rr))
    stdrns.append(numpy.std(rr))
    
    
hpvs = {}

lpvs = {}

for i in list(realnenos.keys()):
    hco = 0
    lco = 0
    for j in randomnenos[i]:
        if j >= realnenos[i]:
            hco = hco +1
        if j <= realnenos[i]:
            lco = lco +1
    hpvs[i] = hco/iters
    lpvs[i] = lco/iters
    

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)

mplt.hist(nsizes,bins=100) 
#mplt.xticks(range(max(hitnos)+1))
#mplt.xlim([-1, max(hitnos)+1])  
mplt.xlabel('Size of neighborhood')
mplt.ylabel('Number of SNPs')
mplt.savefig('T2DMSNPHist.tiff',dpi = 600)
mplt.show()

# sample from pnet the same size and compare degree dists
# plot deg vals together with randoms
# calculate p vals for min median and max values
print('real values')
minpv = 0
medpv = 0
maxpv = 0
rds = []
for vi in vvp:
    rds.append(len(pnet[vi]))
rds.sort()

rands = [0]*len(vvp)
for i in range(iters):
    randgs = random.sample(pl,len(vvp))
    thisrs = []
    for ri in randgs:
        thisrs.append(len(pnet[ri]))
    thisrs.sort()
    for j in range(len(thisrs)):
        rands[j] = rands[j] + thisrs[j]
        
    if min(thisrs) <= min(rds):
        minpv = minpv + 1
    if numpy.median(thisrs) >= numpy.median(rds):
        medpv = medpv + 1
    if max(thisrs) >= max(rds):
        maxpv = maxpv + 1

minpv = minpv/iters
medpv = medpv/iters
maxpv = maxpv/iters        

for i in range(len(rands)):
    rands[i] = rands[i]/iters
    
mplt.hist([rds,rands],bins=50,color=['y','c'])
mplt.xlabel('Degree value')
mplt.ylabel('Frequency')
mplt.legend(['T2DM','Random'])
mplt.savefig('T2DMGenes_vs_Rands.tiff',dpi = 600)
mplt.show()



