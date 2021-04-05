#simulate a random individual for SNP events
import getVarGene
import random
import matplotlib.pyplot as mplt

noSNPs = 4000000
itersno = 50000
chd = getVarGene.chrsize()
pnet = getVarGene.biogridInts()
impgl = getVarGene.keggMembers()
chl = list(chd.keys())
tot = 0
for i in chl:
    tot = tot + chd[i]

#no of snps per chr
snd = {}
snpratio = noSNPs/tot
snptot = 0
for i in chl:
    snd[i] = int(chd[i]*snpratio)
    snptot = snptot + snd[i]

while snptot > noSNPs:
    j = random.sample(chl,1)
    snd[j[0]] = snd[j[0]] - 1
    snptot = snptot - 1
     
while snptot < noSNPs:
    j = random.sample(chl,1)
    snd[j[0]] = snd[j[0]] + 1
    snptot = snptot + 1

#get real SNP locations (clinically associated)
#keep locs per chr in a dic
vad = getVarGene.ebigwas()
vl = list(vad.keys())
scloc = {}
for i in chl:
    scloc[i] = []
for i in vl:
    ii = i.split(' ')
    if int(ii[1]) not in scloc[ii[0]]:
        scloc[ii[0]].append(int(ii[1]))

#for a number of iterations get random locs and check intersection
#keep number of intersections
#hitgenes and hitneighs is repetitive (gene may occur more than once)
#hitimps are unique (only once per gene)
hitnos = []
hitlocs = []
hitgenes = []
hitimpgenes = []
hitneighs = []
hitimpneighs = []

for i in range(itersno):
    #get random SNP locations
    print(i)
    thisno = 0
    thisloc = []
    thisgene = []
    thisneigh = []
    
    for ci in chl:
        rc = random.sample(range(1,chd[ci]),snd[ci])
        comlocs = list(set(rc) & set(scloc[ci]))
        
        if comlocs != []:
            thisno = thisno + len(comlocs)
            for li in comlocs:
                chrloc = str(ci) + ' ' + str(li)
                thisloc.append(chrloc)
                thisgene.extend(vad[chrloc])
        
    hitnos.append(thisno)
    hitlocs.append(thisloc)
    hitgenes.append(thisgene)

    thisimpgene = list(set(impgl) & set(thisgene))
    hitimpgenes.append(thisimpgene)
    
    gs = list(set(thisgene))
    for ni in gs:
        if ni in pnet.keys():
            thisneigh.extend(pnet[ni])
    hitneighs.append(thisneigh)
        
    thisimpneigh = list(set(impgl) & set(thisneigh))
    hitimpneighs.append(thisimpneigh)    
    
import numpy    
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'xx-large',
          'figure.figsize': (10, 8),
         'axes.labelsize': 'xx-large',
         'axes.titlesize':'xx-large',
         'xtick.labelsize':'xx-large',
         'ytick.labelsize':'xx-large'}
pylab.rcParams.update(params)

bins = numpy.arange(max(hitnos)+2) - 0.5
mplt.hist(hitnos,bins,weights=numpy.ones(len(hitnos)) / len(hitnos)) 
mplt.xticks(range(max(hitnos)+1))
mplt.xlim([-1, max(hitnos)+1])  
mplt.xlabel('Number of hits')
mplt.ylabel('Fraction of random genomes')
mplt.savefig('Hitsfig.tiff',dpi = 300)
mplt.show()

hitgc = []
hitgenesunq = []
hitgcunq = []
for i in hitgenes:
    hitgc.append(len(i))
    hitgenesunq.append(list(set(i)))
    hitgcunq.append(len(set(i)))

bins = numpy.arange(max(hitgc)+2) - 0.5
mplt.hist(hitgc,bins, weights=numpy.ones(len(hitgc)) / len(hitgc))  
mplt.xticks(range(max(hitgc)+1))
mplt.xlim([-1, max(hitgc)+1])   
mplt.xlabel('Number of mapped genes')
mplt.ylabel('Fraction of random genomes')
mplt.savefig('Hitsfig_gene.tiff',dpi = 300)
mplt.show()


bins = numpy.arange(max(hitgcunq)+2) - 0.5
mplt.hist(hitgcunq,bins, weights=numpy.ones(len(hitgcunq)) / len(hitgcunq))  
mplt.xticks(range(max(hitgcunq)+1))
mplt.xlim([-1, max(hitgcunq)+1])   
mplt.xlabel('Number of mapped unique genes')
mplt.ylabel('Fraction of random genomes')
mplt.savefig('Hitsfig_geneunique.tiff',dpi = 300)
mplt.show()


hitnc = []
hitneighsunq = []
hitncunq = []
for i in hitneighs:
    hitnc.append(len(i))
    hitneighsunq.append(list(set(i)))
    hitncunq.append(len(set(i)))

bins = numpy.arange(0,max(hitnc)+2,10) - 0.5
mplt.hist(hitnc,bins, weights=numpy.ones(len(hitnc)) / len(hitnc)) 
mplt.xticks(range(0,max(hitnc)+1,150))
mplt.xlim([-1, max(hitnc)+1])    
mplt.xlabel('Size of gene neighborhood')
mplt.ylabel('Fraction of random genomes')
mplt.savefig('Hitsfig_neigh.tiff',dpi = 300)
mplt.show()


bins = numpy.arange(0,max(hitncunq)+2,10) - 0.5
mplt.hist(hitncunq,bins, weights=numpy.ones(len(hitncunq)) / len(hitncunq)) 
mplt.xticks(range(0,max(hitncunq)+1,150))
mplt.xlim([-1, max(hitncunq)+1])    
mplt.xlabel('Size of unique gene neighborhood')
mplt.ylabel('Fraction of random genomes')
mplt.savefig('Hitsfig_neighunique.tiff',dpi = 300)
mplt.show()


hitic = []
for i in hitimpgenes:
    hitic.append(len(i))

bins = numpy.arange(max(hitic)+2) - 0.5
mplt.hist(hitic,bins, weights=numpy.ones(len(hitic)) / len(hitic))  
mplt.xticks(range(max(hitic)+1))
mplt.xlim([-1, max(hitic)+1])   
mplt.xlabel('Number of mapped IS Pathway genes')
mplt.ylabel('Fraction of random genomes')
mplt.savefig('Hitsfig_ISPgenes.tiff',dpi = 300)
mplt.show()

hitinc = []
hitimpneighsunq = []
hitincunq = []
for i in hitimpneighs:
    hitinc.append(len(i))
    hitimpneighsunq.append(list(set(i)))
    hitincunq.append(len(set(i)))

bins = numpy.arange(max(hitinc)+2) - 0.5
mplt.hist(hitinc,bins, weights=numpy.ones(len(hitinc)) / len(hitinc)) 
mplt.xticks(range(0,max(hitinc)+1,2))
mplt.xlim([-1, max(hitinc)+1])    
mplt.xlabel('Size of IS Pathway neighborhood')
mplt.ylabel('Fraction of random genomes')
mplt.savefig('Hitsfig_ISPgenesneigh.tiff',dpi = 300)
mplt.show()


bins = numpy.arange(max(hitincunq)+2) - 0.5
mplt.hist(hitincunq,bins, weights=numpy.ones(len(hitincunq)) / len(hitincunq)) 
mplt.xticks(range(0,max(hitincunq)+1,2))
mplt.xlim([-1, max(hitincunq)+1])    
mplt.xlabel('Size of unique IS Pathway neighborhood')
mplt.ylabel('Fraction of random genomes')
mplt.savefig('Hitsfig_ISPgenesneighunique.tiff',dpi = 300)
mplt.show()



allhitgenes = sum(hitgenes,[])
allhitgenesunq = sum(hitgenesunq,[])
allhitimpgenes = sum(hitimpgenes,[])
allhitneighs = sum(hitneighs,[])
allhitimpneighs = sum(hitimpneighs,[])
allhitneighsunq = sum(hitneighsunq,[])
allhitimpneighsunq = sum(hitimpneighsunq,[])

# get frequently hit genes
print('hit genes: ')
for i in set(allhitgenes):
    if allhitgenes.count(i)>(itersno/100):
        print(i)
        print(allhitgenes.count(i))
        
print('hit unique genes: ')
for i in set(allhitgenesunq):
    if allhitgenesunq.count(i)>(itersno/100):
        print(i)
        print(allhitgenesunq.count(i))

print('hit ISP genes more than 1: ')        
for i in set(allhitimpgenes):
    if allhitimpgenes.count(i)>1:
        print(i)
        print(allhitimpgenes.count(i))

print('All hit ISP genes: ')        
print(set(allhitimpgenes))
        
print('hit gene neighbors: ')          
for i in set(allhitneighs):
    if allhitneighs.count(i)>(itersno/10):
        print(i)
        print(allhitneighs.count(i))
        
print('hit unique gene neighbors: ')          
for i in set(allhitneighsunq):
    if allhitneighsunq.count(i)>(itersno/10):
        print(i)
        print(allhitneighsunq.count(i))        
    

print('hit ISP gene neighbors: ')         
for i in set(allhitimpneighs):
    if allhitimpneighs.count(i)>(itersno/100):
        print(i)
        print(allhitimpneighs.count(i))
        
print('hit unique ISP gene neighbors: ')         
for i in set(allhitimpneighsunq):
    if allhitimpneighsunq.count(i)>(itersno/100):
        print(i)
        print(allhitimpneighsunq.count(i))

