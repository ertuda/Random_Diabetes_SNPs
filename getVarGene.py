#Ertugrul Dalgic,2019
#import os
#os.chdir('C:/Users/Dell/Documents/Projects/DiabPatGene/')
import random
import copy

def ebigwas(cvfile='DM_GWAS_2019-03-21-EFO_0001360-withChildTraits.tsv'):
    # collect lists of genes from ClinVar
    vardic = {}
    
    a = open(cvfile)
    al = a.readlines()
    a.close()
    # first (0) line contains colnames
    # get indiced for columns
    nameslist = al[0].split('\t')
    gind = nameslist.index('MAPPED_GENE')
    c1ind = nameslist.index('CHR_ID')
    c2ind = nameslist.index('CHR_POS')
    sind = nameslist.index('SNPS')
    
    for i in range(1,len(al)):
        plist = al[i].split('\t')
        
        if plist[c1ind] != '':
            cinf = plist[c1ind].split(';')
            if len(cinf) > 1:
                pinf = plist[c2ind].split(';')
                ginf = plist[gind].split(';')
                for ci in range(len(cinf)):
                    gg = ginf[ci].split(' - ')
                    for ggi in gg:
                        theloc = cinf[ci].strip()+' '+pinf[ci].strip()
                        if theloc not in vardic.keys():
                            vardic[theloc] = []
                        vardic[theloc].append(ggi.strip())
                                                    
            else:
                cinf2 = plist[c1ind].split(' x ')
                if len(cinf2) > 1:
                    pinf = plist[c2ind].split(' x ')
                    ginf = plist[gind].split(' x ')
                    for ci in range(len(cinf2)):
                        gg = ginf[ci].split(' - ')
                        for ggi in gg:
                            theloc = cinf2[ci].strip()+' '+pinf[ci].strip()
                            if theloc not in vardic.keys():
                                vardic[theloc] = []
                            vardic[theloc].append(ggi.strip())
                    
                    
                else:
                    theloc = plist[c1ind].strip()+' '+plist[c2ind].strip()
                    if theloc not in vardic.keys():
                        vardic[theloc] = []
                        
                    if ',' in plist[gind]:
                        gg = plist[gind].split(',')
                        for ggi in gg:
                            vardic[theloc].append(ggi.strip())
                            
                    elif ' - ' in plist[gind]:
                        gg = plist[gind].split(' - ')
                        for ggi in gg:
                            vardic[theloc].append(ggi.strip())
                    
                    else:
                        vardic[theloc].append(plist[gind].strip())
                    
    
    # make unique entries
    vkeys = list(vardic.keys())
    for vi in vkeys:
        vardic[vi] = list(set(vardic[vi]))                    
                               
    return vardic


def keggMembers(kfile='KEGG_INSRPathwayMembers.txt',idfile = 'HGNC_03282019.txt'):   
    # id to symbol dic
    idname = {}
    a = open(idfile,encoding="utf8")
    al = a.readlines()
    a.close()
    # first (0) line contains colnames
    # get indiced for columns
    nameslist = al[0].split('\t')
    eind = nameslist.index('entrez_id')
    sind = nameslist.index('symbol')
    
    for i in range(1,len(al)):
        plist = al[i].split('\t')
        if plist[eind] != '' and plist[sind] != '':
            idname[plist[eind]] = plist[sind]
        
    # get kegg pathway member genes
    # there is no name row
    keggpathlist = []    
    a = open(kfile)
    al = a.readlines()
    a.close()
    
    for i in range(1,len(al)):
        plist = al[i].split()
        if plist[0] in idname.keys():
            keggpathlist.append(idname[plist[0]])
        
    return keggpathlist


def chrsize(chrsizefile='GRCh38p12_Chr_Sizes_noY.txt'):
    # get chr sizes
    chrsi = {}
    
    a = open(chrsizefile)
    al = a.readlines()
    a.close()
    # first (0) line contains colnames
    # get indiced for columns
    nameslist = al[0].split('\t')
    cind = nameslist.index('Chr')
    sind = nameslist.index('Size(bp)\n')
    
    for i in range(1,len(al)):
        plist = al[i].split('\t')
        chrsi[plist[cind]] = int(plist[sind].strip())
    return chrsi

def biogridInts(intfile='BIOGRID_3.5.170_Human.txt',idfile = 'HGNC_03282019.txt'):
    # collect neighborhoods of genes from undirected network
    hnet = {}
    
    # id to symbol dic
    idname = {}
    a = open(idfile,encoding="utf8")
    al = a.readlines()
    a.close()
    # first (0) line contains colnames
    # get indiced for columns
    nameslist = al[0].split('\t')
    eind = nameslist.index('entrez_id')
    sind = nameslist.index('symbol')
    
    for i in range(1,len(al)):
        plist = al[i].split('\t')
        if plist[eind] != '' and plist[sind] != '':
            idname[plist[eind]] = plist[sind]
        
    
    a = open(intfile)
    al = a.readlines()
    a.close()
    # first (0) line contains colnames
    # get indiced for columns
    nameslist = al[0].split('\t')
    gAind = nameslist.index('Entrez.Gene.Interactor.A') + 1
    gBind = nameslist.index('Entrez.Gene.Interactor.B') + 1
    
    for i in range(1,len(al)):
        plist = al[i].split('\t')
        if plist[gAind] in idname.keys() and plist[gBind] in idname.keys():
            gA = idname[plist[gAind]]
            gB = idname[plist[gBind]]
 
            if gA not in hnet.keys():
                hnet[gA] = []
            if gB not in hnet.keys():
                hnet[gB] = []
            if gB not in hnet[gA]:
                hnet[gA].append(gB)
            if gA not in hnet[gB]:
                hnet[gB].append(gA)
            
    return hnet



def netShuffle(givnet):
    # randomize network which is given as dic of neighbors
    # self interactions are allowed
    # keep the number of self interactions same
    # keep degree distribution same 
    # a sorted list from higher degs to lower 
    
    shufnet = {}
    degdic = {}
    gl = list(givnet.keys())
    
    # count self ints and redistribute
    selfcount = 0
    for i in gl:
        if i in givnet[i]:
            selfcount = selfcount + 1
    newselfints = random.sample(range(len(gl)),selfcount)
    
    for i in range(len(gl)):
        gn = gl[i]
        shufnet[gn] = []
        if i in newselfints:
            shufnet[gn].append(gn)
            
        degdic[gn] = len(givnet[gn])

    sortdegs = list(degdic.values())
    sortdegs.sort(reverse=1)

    # randomize gene list and reassign degree vals
    random.shuffle(gl)
    glsc = copy.copy(gl)
    
    # use onedegs lastly !
    onedegs = []
    
    for i in range(len(gl)):
        gn = gl[i]
        dg = sortdegs[i]
        if gn in shufnet[gn]:
            dg = dg-1
        degdic[gn] = dg
        if dg == 0:
            glsc.remove(gn)
        if dg == 1:
            onedegs.append(gn)
    notonedegs = copy.copy(list(set(glsc)-set(onedegs)))
          
    # randomly assign neighbors (start from higher connected ones)
    for i in range(len(gl)):
        #print(glsc)
        #print(onedegs)
        #print(notonedegs)
        gn = gl[i]
        #print(i)
        #print(gn)
        if degdic[gn] > 0:
            glsc.remove(gn)
            #print(glsc)
            if gn in onedegs:
                onedegs.remove(gn)
            if gn in notonedegs:
                notonedegs.remove(gn)
            
            subst  = degdic[gn]-len(notonedegs)
            #print(subst)
            if subst <= 0:
                rg = random.sample(notonedegs,degdic[gn])
            if subst > 0:
                rg = notonedegs
                rgx = random.sample(onedegs,subst)
                rg.extend(rgx)
            
            rg = copy.copy(rg)
            shufnet[gn].extend(rg)
            degdic[gn] = 0
            for ri in rg:
                shufnet[ri].append(gn)
                degdic[ri] = degdic[ri]-1
                if degdic[ri] == 0:
                    glsc.remove(ri)
                    if ri in onedegs:
                        onedegs.remove(ri)
                    if ri in notonedegs:
                        notonedegs.remove(ri)
        if glsc == []:
            break
       
    return shufnet         
               

def vargs(cvfile='DM_clinvar_03202019.txt'):
    # collect lists of genes from ClinVar
    vargen = []
    
    a = open(cvfile)
    al = a.readlines()
    a.close()
    # first (0) line contains colnames
    # get indiced for columns
    nameslist = al[0].split('\t')
    gind = nameslist.index('Gene(s)')
    
    for i in range(1,len(al)):
        plist = al[i].split('\t')
        vargen.extend(plist[gind].split('|'))
    vargen = list(set(vargen))
    vargen.remove('')
                               
    return vargen