import sys
import os
import subprocess
import numpy as np
import math
import time
import random
from copy import deepcopy
import argparse
import multiprocessing as mp
from functools import partial
import gzip
import mmap

#*****************************************************************
#*****************************************************************
def emissionprob(scaledmutrate):
    emissionmat = np.array([[0.0 for j in range(3)] for i in range(3)])
    #emissionmat[observed genotype][combined reference genotype]
    emissionmat[0,0] = (1-scaledmutrate)**2
    emissionmat[2,2] = (1-scaledmutrate)**2
    emissionmat[0,1] = scaledmutrate*(1-scaledmutrate)
    emissionmat[2,1] = scaledmutrate*(1-scaledmutrate)
    emissionmat[0,2] = scaledmutrate**2
    emissionmat[2,0] = scaledmutrate**2
    emissionmat[1,0] = 2*scaledmutrate*(1-scaledmutrate)
    emissionmat[1,1] = scaledmutrate**2 + (1-scaledmutrate)**2
    emissionmat[1,2] = 2*scaledmutrate*(1-scaledmutrate)
    return(emissionmat)
#*****************************************************************
#*****************************************************************
def transprob_unbiased(happop, localrate):
    return((np.exp(-localrate/happop) + ((1-np.exp(-localrate/happop))/happop), ((1-np.exp(-localrate/happop))/happop)))
#*****************************************************************
#*****************************************************************
def row1Dto2D(k, totvals):
    return(int(-0.5+0.5*math.sqrt(1+8*k)))
#*****************************************************************
#*****************************************************************
def map1Dto2D(k, totvals):
    rowID = row1Dto2D(k,totvals)
    colID = k - rowID*(rowID+1)/2
    return((rowID,colID))
#*****************************************************************
#*****************************************************************
def generate_recomb(prefix, effpop, chromID, recombrate, obssamplefile):
    infile = open(obssamplefile)
    snplist = []
    for line in infile:
        terms = line.strip().split("\t")
        snplist.append(int(terms[1]))
    snplist = sorted(snplist)
    infile.close()
    distlist = [(snplist[snpindex] - snplist[snpindex-1])/(10.0**6) for snpindex in range(1,len(snplist),1)]
    recombfile = prefix+"_"+chromID+"_Linear_distance_recombination.tsv"
    outfile = open(recombfile,'w')
    for dist in distlist:
        outfile.write(str(4*effpop*dist*recombrate)+"\n")
    outfile.close()
    return(recombfile)
#*****************************************************************
#*****************************************************************
def update_func(pmat, lhaps, delta, si, nsnps, tolerance, index_tuple):
    localp = deepcopy(pmat)
    i = index_tuple[0]
    j = index_tuple[1]
    localp[int(i*(i+1)/2):int((i+1)*(i+2)/2)] += delta
    colargs = [(k+1)*k/2 + j for k in range(j,lhaps,1)]
    for k in colargs:
        localp[int(k)] += delta
    maxval = np.amax(localp)
    return (maxval, np.argwhere(localp >= maxval+tolerance*maxval))
#*****************************************************************
#*****************************************************************
def process_SNPs_position_specific(dsfilename, filename, observedsample, refpop):
    #Extract the SNP IDs from the observed sample
    prefilterds = ".".join(dsfilename.split(".")[:-1])+".prefiltered.txt"
    bcfcommand = "bcftools view -Ov -R "+observedsample+" "+filename+" | awk '$1 ~ /CHROM/ || $1 !~ /##/' > "+prefilterds
    subprocess.call(bcfcommand,shell=True)
    obsmatchdict = {}
    mutratedict = {}
    obsfile = open(observedsample,'r')
    for line in obsfile:
        terms = line.strip().split("\t")
        finmatch = "_".join([terms[0],terms[1],terms[3]])
        obsmatchdict[finmatch] = terms[4]
        mutratedict[finmatch] = (terms[5],terms[6])
    obsfile.close()
    ####Remove Multi-Allelic SNPs and deletions/insertions
    obsgtypelist = []
    mutratelist = []
    snplist = []
    infile = open(prefilterds)
    outfile = open(dsfilename,'w')
    count=0
    nsnps = 0
    for line in infile:
        terms = line.strip().split("\t")
        finmatch = "_".join([terms[0],terms[1],terms[4]])
        if count == 0:
            outfile.write(line)
            terms = line.strip().split("\t")
            samplelistprime = terms[9:]
        elif (len(terms[4].split(","))==1) and ("VT=SNP" in terms[7]) and (finmatch in obsmatchdict.keys()):
            outfile.write(line)
            snplist.append("chr"+finmatch)
            nsnps += 1
            obsgtypelist.append(obsmatchdict[finmatch])
            mutratelist.append(mutratedict[finmatch])
        count += 1
    infile.close()
    outfile.close()
    ###Clean-up for specific pathology in Sample IDs for 1kGenomes vcf files
    sampletemp = []
    for item in samplelistprime:
        if len(item)==14:
            sampletemp.append(item[:7])
            sampletemp.append(item[7:])
        else:
            sampletemp.append(item)

    samplelist = [sampletemp[int(i/2)]+"_A" if i%2 ==0 else sampletemp[int((i-1)/2)]+"_B" for i in range(2*len(sampletemp))][:2*refpop]
    return((samplelist,obsgtypelist,mutratelist,nsnps,snplist))
#*****************************************************************
#*****************************************************************
def process_SNPs(dsfilename, filename, observedsample, refpop):
    #Extract the SNP IDs from the observed sample
    prefilterds = ".".join(dsfilename.split(".")[:-1])+".prefiltered.txt"
    bcfcommand = "bcftools view -Ov -R "+observedsample+" "+filename+" | awk '$1 ~ /CHROM/ || $1 !~ /##/' > "+prefilterds
    subprocess.call(bcfcommand,shell=True)

    obsmatchdict = {}
    obsfile = open(observedsample,'r')
    for line in obsfile:
        terms = line.strip().split("\t")
        obsmatchdict["_".join([terms[0],terms[1],terms[3]])] = terms[4]
    obsfile.close()
    ####Remove Multi-Allelic SNPs and deletions/insertions
    obsgtypelist = []
    snplist = []
    infile = open(prefilterds)
    outfile = open(dsfilename,'w')
    count=0
    nsnps = 0
    for line in infile:
        terms = line.strip().split("\t")
        finmatch = "_".join([terms[0],terms[1],terms[4]])
        if count == 0:
            outfile.write(line)
            terms = line.strip().split("\t")
            samplelistprime = terms[9:]
        elif (len(terms[4].split(","))==1) and ("VT=SNP" in terms[7]) and (finmatch in obsmatchdict.keys()):
            outfile.write(line)
            snplist.append("chr"+finmatch)
            nsnps += 1
            obsgtypelist.append(obsmatchdict[finmatch])
        count += 1
    infile.close()
    outfile.close()
    ###Clean-up for specific pathology in Sample IDs for 1kGenomes vcf files
    sampletemp = []
    for item in samplelistprime:
        if len(item)==14:
            sampletemp.append(item[:7])
            sampletemp.append(item[7:])
        else:
            sampletemp.append(item)

    samplelist = [sampletemp[int(i/2)]+"_A" if i%2 ==0 else sampletemp[int((i-1)/2)]+"_B" for i in range(2*len(sampletemp))][:2*refpop]
    return((samplelist,obsgtypelist,nsnps,snplist))
#*****************************************************************
#*****************************************************************
def selection_func(samplelist, subgroup):
        sampleindices = list(range(len(samplelist)))
        random.shuffle(sampleindices)
        select_sample = [sampleindices[i*subgroup:(i+1)*subgroup] for i in range((len(sampleindices)+subgroup-1) // subgroup)]
        return(select_sample)
#*****************************************************************
#*****************************************************************
def read_backtrace_in_reference(prefix, btfilename, termbtcollapse, lhaps, samplelist, snplist, chromID):
    infile = open(btfilename,'rb')
    mm = mmap.mmap(infile.fileno(), 0, access = mmap.ACCESS_READ)
    gzfile = gzip.GzipFile(mode="r", fileobj=mm)
    trajfile = prefix+"_"+chromID + "_Best_trajectories.tsv"

    seqindex = termbtcollapse.split(";")


    numstates = len(seqindex)
    seq2D = [int(item) for item in seqindex]
    finsamples = [samplelist[2*item]+"_"+samplelist[2*item+1] for item in seq2D]
    with open(trajfile,'w') as outfile:
        outfile.write(snplist[-1]+"\t"+"\t".join([samplelist[2*item]+"_"+samplelist[2*item+1] for item in seq2D])+"\n")
        for snpindex in range(len(snplist)-1,0,-1):
          seqindexdict = {}

          outfile.write(snplist[snpindex-1]+"\t")
          gzfile.seek(0)
          for lindex in range(2*(snpindex-1)+1):
              gzfile.readline()
          btflat = gzfile.readline().strip().split(b",")
          btflat = np.array([bitem.decode('utf-8') for bitem in btflat])

          seqlist = []
          for poss in seqindex:

              posskey = samplelist[2*int(poss)]+"_"+samplelist[2*int(poss)+1]

              seq = btflat[int(poss)]

              seqkey = samplelist[2*int(seq)]+"_"+samplelist[2*int(seq)+1]
              seqindexdict[posskey] = seqkey

              seqlist += seq.split(";")
              seqindex = list(set(seqlist))
          for seqpair in seqindexdict.items():
              outfile.write(":".join(seqpair)+"\t")
          outfile.write("\n")
    infile.close()
    rmcall = f"rm {btfilename}"
    subprocess.call(rmcall,shell=True)
    return(trajfile)
#*****************************************************************
#*****************************************************************
def read_backtrace_exact(prefix, btfilename, termbtcollapse, lhaps, samplelist, snplist, chromID):

    infile = open(btfilename,'rb')
    mm = mmap.mmap(infile.fileno(), 0, access = mmap.ACCESS_READ)
    gzfile = gzip.GzipFile(mode="r", fileobj=mm)
    trajfile = prefix+"_"+chromID + "_Best_trajectories.tsv"
    seqindex = termbtcollapse.split(";")
    seq2D = [map1Dto2D(int(item),lhaps) for item in seqindex]
    with open(trajfile,'w') as outfile:
        outfile.write(snplist[-1]+"\t"+"\t".join([samplelist[int(item[0])]+"_"+samplelist[int(item[1])] for item in seq2D])+"\n")
        for snpindex in range(len(snplist)-1,0,-1):
          seqindexdict = {}

          outfile.write(snplist[snpindex-1]+"\t")
          gzfile.seek(0)
          for lindex in range(2*(snpindex-1)+1):
              gzfile.readline()
          btflat = gzfile.readline().strip().split(b",")
          btflat = np.array([bitem.decode('utf-8') for bitem in btflat])


          seqlist = []
          for poss in seqindex:
              poss2D = map1Dto2D(int(poss),lhaps)
              posskey = samplelist[int(poss2D[0])]+"_"+samplelist[int(poss2D[1])]

              seq = btflat[int(poss)]
              seq2D = [map1Dto2D(int(item.split(".")[0]),lhaps) for item in seq.split(";")]
              seqkey = ";".join([samplelist[int(item[0])]+"_"+samplelist[int(item[1])] for item in seq2D])
              seqindexdict[posskey] = seqkey

              seqlist += seq.split(";")
              seqindex = list(set(seqlist))
          for seqpair in seqindexdict.items():
              outfile.write(":".join(seqpair)+"\t")
          outfile.write("\n")


    infile.close()
    rmcall = f"rm {btfilename}"
    subprocess.call(rmcall,shell=True)
    return(trajfile)
#*****************************************************************
#*****************************************************************
def read_backtrace_iterative(prefix, btfilename, termbtcollapse, lhaps, samplelist, samplelisttot, snplist, chromID, writeflag):
    besthaps = []
    infile = open(btfilename,'rb')
    mm = mmap.mmap(infile.fileno(), 0, access = mmap.ACCESS_READ)
    gzfile = gzip.GzipFile(mode="r", fileobj=mm)
    trajfile = prefix+"_"+chromID + "_Best_trajectories.tsv"

    seqindex = termbtcollapse.split(";")

    seq2D = [map1Dto2D(int(item),lhaps) for item in seqindex]
    besthaps.extend(seq2D)

    if(writeflag=="final"):
        outfile = open(trajfile,'w')
        outfile.write(snplist[-1]+"\t"+"\t".join([samplelist[int(item[0])]+"_"+samplelist[int(item[1])] for item in seq2D])+"\n")

    for snpindex in range(len(snplist)-1,0,-1):
      seqindexdict = {}

      if(writeflag=="final"):
          outfile.write(snplist[snpindex-1]+"\t")

      gzfile.seek(0)
      for lindex in range(2*(snpindex-1)+1):
          gzfile.readline()
      btflat = gzfile.readline().strip().split(b",")
      btflat = np.array([bitem.decode('utf-8') for bitem in btflat])

      seqlist = []
      for poss in seqindex:
          poss2D = map1Dto2D(int(poss),lhaps)
          posskey = samplelist[int(poss2D[0])]+"_"+samplelist[int(poss2D[1])]

          seq = btflat[int(poss)]
          seq2D = [map1Dto2D(int(item.split(".")[0]),lhaps) for item in seq.split(";")]
          besthaps.extend(seq2D)
          besthaps = list(set(besthaps))
          seqkey = ";".join([samplelist[int(item[0])]+"_"+samplelist[int(item[1])] for item in seq2D])
          seqindexdict[posskey] = seqkey

          seqlist += seq.split(";")
          seqindex = list(set(seqlist))

      if(writeflag=="final"):
          for seqpair in seqindexdict.items():
              outfile.write(":".join(seqpair)+"\t")
          outfile.write("\n")

    if(writeflag=="final"):
        outfile.close()

    infile.close()
    besthaps = list(set([samplelist[int(item)] for item2D in besthaps for item in item2D]))
    besthapindices = [i for i, val in enumerate(samplelisttot) if val in besthaps]

    rmcall = f"rm {btfilename}"
    subprocess.call(rmcall,shell=True)
    return(besthapindices)
#*****************************************************************
#*****************************************************************
def run_hmm_exact(prefix, dsfilename, filename, observedsample, chromID, refpop, recombfile, scaledmutrate, tolerance, numproc, posspecific):
    indep_p = 0
    emissionmat = emissionprob(scaledmutrate)
    lowerlimit = sys.float_info.min

    prefilterds = ".".join(dsfilename.split(".")[:-1])+".prefiltered.txt"
    bcfcommand = "bcftools view -Ov -R "+observedsample+" "+filename+" | awk '$1 ~ /CHROM/ || $1 !~ /##/' > "+prefilterds
    subprocess.call(bcfcommand,shell=True)
    obsmatchdict = {}
    if posspecific=="True":
        mutratedict = {}
    obsfile = open(observedsample,'r')
    for line in obsfile:
        terms = line.strip().split("\t")
        finmatch = "_".join([terms[0],terms[1],terms[3]])
        obsmatchdict[finmatch] = terms[4]
        if posspecific=="True":
            mutratedict[finmatch] = (terms[5],terms[6])
    obsfile.close()
    ####Remove Multi-Allelic SNPs and deletions/insertions
    obsgtypelist = []
    if posspecific=="True":
        mutratelist = []
    snplist = []
    infile = open(prefilterds)
    outfile = open(dsfilename,'w')
    count=0
    nsnps = 0
    for line in infile:
        terms = line.strip().split("\t")
        finmatch = "_".join([terms[0],terms[1],terms[4]])
        if count == 0:
            outfile.write(line)
            terms = line.strip().split("\t")
            samplelistprime = terms[9:]
        elif (len(terms[4].split(","))==1) and ("VT=SNP" in terms[7]) and (finmatch in obsmatchdict.keys()):
            outfile.write(line)
            snplist.append("chr"+finmatch)
            nsnps += 1
            obsgtypelist.append(obsmatchdict[finmatch])
            if posspecific=="True":
                mutratelist.append(mutratedict[finmatch])
        count += 1
    infile.close()
    outfile.close()
    if posspecific=="False":
        mutratelist = None
    genposs = ["0","1"]
    #EXTRACT the haplotypes of all the individuals from 1kGenomes
    btfilename = prefix+"_"+chromID+"_BackTrace_Sequences.txt.gz"
    btfile = gzip.open(btfilename,'wb')
    infile = open(dsfilename,'r')
    inrecomb = open(recombfile, 'r')
    mm = mmap.mmap(infile.fileno(), 0, access = mmap.ACCESS_READ)
    mm.seek(0)

    for si in range(-1,nsnps,1):
        print(si)
        flag = 0
        terms = mm.readline().strip().split(b"\t")
        terms = [bitem.decode('utf-8') for bitem in terms]
        if(terms[0][:6]=="#CHROM"):

            samplelistprime = terms[9:]
            ###Clean-up for specific pathology in Sample IDs for 1kGenomes vcf files
            sampletemp = []
            for item in samplelistprime:
                if len(item)==14:
                    sampletemp.append(item[:7])
                    sampletemp.append(item[7:])
                else:
                    sampletemp.append(item)


            samplelistprime = sampletemp[:refpop]
            samplelist = [samplelistprime[int(i/2)]+"_A" if i%2 ==0 else samplelistprime[int((i-1)/2)]+"_B" for i in range(2*len(samplelistprime))]
            #********************
            indexlist = [(i,j) for i in range(len(samplelist)) for j in range(0,i+1,1)]
        elif terms[0][0]!="#":

            if (si>0):
                localrate = float(inrecomb.readline())
                transmat = transprob_unbiased(2*refpop,localrate)

                transmat = [math.log(lowerlimit) if item == 0 else math.log(item) for item in transmat]
                ton = transmat[0]
                toff = transmat[1]



            obsgtype = obsgtypelist[si]

            if mutratelist == None:
                emslice = [lowerlimit if item==0 else item for item in emissionmat[int(obsgtype)]]
                em0 = emslice[0]
                em1 = emslice[1]
                em2 = emslice[2]
                print("Non-position-specific")
            else:
                mutrate = mutratelist[si]
                mutation_tuples = [(int(item.split(":")[0]),float(item.split(":")[1])) for item in mutrate]
                sumprob = sum([tup[1] for tup in mutation_tuples])
                mutation_tuples.append((int(obsgtype), 1-sumprob))
                mutation_tuples.sort(key=lambda x: x[0])
                emslice = [item[1] for item in mutation_tuples]
                emslice = [lowerlimit if item==0 else item for item in emslice]
                em0 = emslice[0]
                em1 = emslice[1]
                em2 = emslice[2]


            freqs = terms[7].split(";")

            allelefreqprime = [freq.split("=")[1] for freq in freqs if (freq[:3]=="AF=")][0]
            if len(allelefreqprime.split(","))>1:
                allelefreq = float(allelefreqprime.split(",")[0])
            else:
                allelefreq = float(allelefreqprime)
            p1 = allelefreq
            p0 = 1.0-p1
            haps = terms[9:refpop+9]
            haps = [indhap+"|"+indhap if len(indhap.split("|")) == 1 else indhap for indhap in haps]
            haps = [item for indhap in haps for item in indhap.split("|")]

            if int(obsgtype) == 2:
                SNP_prob = np.log(p1*p1)
            elif int(obsgtype) == 1:
                SNP_prob = np.log(2*p1*p0)
            elif int(obsgtype) == 0:
                SNP_prob = np.log(p0*p0)
            indep_p += SNP_prob
            lhaps = len(haps)

            #Explicit calculation of emission probabilities for all possibilities
            haps0 = np.array([i for i in range(lhaps) if haps[i]=="0"], dtype = np.uint32)
            haps1 = np.array([i for i in range(lhaps) if haps[i]=="1"], dtype = np.uint32)
            hapsm = np.array([i for i in range(lhaps) if (haps[i] != "0") and (haps[i] != "1")], dtype = np.uint32)
            if len(hapsm) == 0:
                flag = 1


            if flag != 1:
                emmm = em0*(p0**2) + em2*(p1**2) + em1*2*(p0*p1)
                emm0 = em0*p0 + em1*p1
                emm1 = em1*p0 + em2*p1
            pvec = np.array([math.log(lowerlimit) for i in range(lhaps) for j in range(0,i+1,1)], dtype = np.float32)

            for i in haps0:
                for j in haps0:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(em0)
                for j in haps1:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(em1)
                for j in hapsm:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emm0)
            for i in haps1:
                for j in haps0:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(em1)
                for j in haps1:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(em2)
                for j in hapsm:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emm1)
            for i in hapsm:
                for j in haps0:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emm0)
                for j in haps1:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emm1)
                for j in hapsm:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emmm)


            #pvec = lemission

            if si == 0:
                reflog = -2*math.log(lhaps)
                pvec += reflog

            else:
                pmat = prevpvec + 2*toff

                delta = ton-toff
                pool = mp.Pool(processes=numproc)

                fix_update=partial(update_func, pmat, lhaps, delta, si, nsnps, tolerance)

                returnlist = pool.map(fix_update, indexlist)

                pool.close()
                pool.join()

                totpvec = np.array([item[0] for item in returnlist],dtype = np.float32)
                btvec = [item[1] for item in returnlist]

                returnlist = []

                pvec += totpvec
                #********************

                btflat = [";".join(np.ndarray.flatten(item).astype(str)) for item in btvec]

                str1 = "SNP_"+str(si)+"\n"
                btfile.write(str1.encode('utf-8'))
                str1 = ",".join(btflat) + "\n"
                btfile.write(str1.encode('utf-8'))

            prevpvec = pvec


    #Termination step
    pmat = np.array(pvec)
    termp = np.amax(pmat)
    termbt = np.argwhere(pmat >= termp+tolerance*termp)

    termbtcollapse = ";".join(np.ndarray.flatten(termbt).astype(str))

    infile.close()
    inrecomb.close()
    btfile.close()

    sum_p = termp + np.log(np.sum(np.exp(-termp + pmat)))
    outfile = open(prefix+"_"+str(chromID)+"_Probability_value.txt",'w')
    outfile.write("Log-Probability value of best-fit trajectories = "+str(termp))
    outfile.write("log(Joint Probability of SNPs) - log(Product of Independent HWE Probabilities of SNP Genotypes) = "+str(sum_p - indep_p))
    outfile.close()
    #********************

    return btfilename, termbtcollapse, lhaps, samplelist, snplist, chromID
#*****************************************************************
#*****************************************************************
def run_hmm_in_reference(prefix, dsfilename, filename, observedsample, chromID, refpop, recombfile, scaledmutrate, tolerance, numproc, posspecific):
    indep_p = 0
    emissionmat = emissionprob(scaledmutrate)
    lowerlimit = sys.float_info.min

    prefilterds = ".".join(dsfilename.split(".")[:-1])+".prefiltered.txt"
    bcfcommand = "bcftools view -Ov -R "+observedsample+" "+filename+" | awk '$1 ~ /CHROM/ || $1 !~ /##/' > "+prefilterds
    subprocess.call(bcfcommand,shell=True)
    obsmatchdict = {}
    if posspecific=="True":
        mutratedict = {}
    obsfile = open(observedsample,'r')
    for line in obsfile:
        terms = line.strip().split("\t")
        finmatch = "_".join([terms[0],terms[1],terms[3]])
        obsmatchdict[finmatch] = terms[4]
        if posspecific=="True":
            mutratedict[finmatch] = (terms[5],terms[6])
    obsfile.close()
    ####Remove Multi-Allelic SNPs and deletions/insertions
    obsgtypelist = []
    if posspecific=="True":
        mutratelist = []
    snplist = []
    infile = open(prefilterds)
    outfile = open(dsfilename,'w')
    count=0
    nsnps = 0
    for line in infile:
        terms = line.strip().split("\t")
        finmatch = "_".join([terms[0],terms[1],terms[4]])
        if count == 0:
            outfile.write(line)
            terms = line.strip().split("\t")
            samplelistprime = terms[9:]
        elif (len(terms[4].split(","))==1) and ("VT=SNP" in terms[7]) and (finmatch in obsmatchdict.keys()):
            outfile.write(line)
            snplist.append("chr"+finmatch)
            nsnps += 1
            obsgtypelist.append(obsmatchdict[finmatch])
            if posspecific=="True":
                mutratelist.append(mutratedict[finmatch])
        count += 1
    infile.close()
    outfile.close()
    if posspecific=="False":
        mutratelist = None
    genposs = ["0","1"]
    #EXTRACT the haplotypes of all the individuals from 1kGenomes
    btfilename = prefix+"_"+chromID+"_BackTrace_Sequences.txt.gz"
    btfile = gzip.open(btfilename,'wb')
    infile = open(dsfilename,'r')
    inrecomb = open(recombfile, 'r')
    mm = mmap.mmap(infile.fileno(), 0, access = mmap.ACCESS_READ)
    mm.seek(0)
    freqlist = []
    #si = 0
    for si in range(-1,nsnps,1):
        flag = 0
        terms = mm.readline().strip().split(b"\t")
        terms = [bitem.decode('utf-8') for bitem in terms]
        if(terms[0][:6]=="#CHROM"):

            samplelistprime = terms[9:]
            ###Clean-up for specific pathology in Sample IDs for 1kGenomes vcf files
            sampletemp = []
            for item in samplelistprime:
                if len(item)==14:
                    sampletemp.append(item[:7])
                    sampletemp.append(item[7:])
                else:
                    sampletemp.append(item)


            samplelistprime = sampletemp[:refpop]
            samplelist = [samplelistprime[int(i/2)]+"_A" if i%2 ==0 else samplelistprime[int((i-1)/2)]+"_B" for i in range(2*len(samplelistprime))]
            indexlist = [i for i in range(len(samplelist))]
        elif terms[0][0]!="#":

            obsgtype = obsgtypelist[si]
            if mutratelist == None:
                emslice = [lowerlimit if item==0 else item for item in emissionmat[int(obsgtype)]]
                em0 = emslice[0]
                em1 = emslice[1]
                em2 = emslice[2]
            else:
                mutrate = mutratelist[si]
                mutation_tuples = [(int(item.split(":")[0]),float(item.split(":")[1])) for item in mutrate]
                sumprob = sum([tup[1] for tup in mutation_tuples])
                mutation_tuples.append((int(obsgtype), 1-sumprob))
                mutation_tuples.sort(key=lambda x: x[0])
                emslice = [item[1] for item in mutation_tuples]
                emslice = [lowerlimit if item==0 else item for item in emslice]
                em0 = emslice[0]
                em1 = emslice[1]
                em2 = emslice[2]
            freqs = terms[7].split(";")

            allelefreqprime = [freq.split("=")[1] for freq in freqs if (freq[:3]=="AF=")][0]
            if len(allelefreqprime.split(","))>1:
                allelefreq = float(allelefreqprime.split(",")[0])
            else:
                allelefreq = float(allelefreqprime)
            freqlist.append(allelefreq)
            p1 = allelefreq
            p0 = 1.0-p1
            haps = terms[9:refpop+9]

            haps = [f"{indhap.split(':')[0].split('/')[0]}|{indhap.split(':')[0].split('/')[1]}" if "/" in indhap.split(':')[0] else indhap.split(':')[0] for indhap in haps]
            if int(obsgtype) == 2:
                SNP_prob = np.log(p1*p1)
            elif int(obsgtype) == 1:
                SNP_prob = np.log(2*p1*p0)
            elif int(obsgtype) == 0:
                SNP_prob = np.log(p0*p0)
            indep_p += SNP_prob
            lhaps = len(haps)

            pvec = np.array([math.log(lowerlimit) for i in range(lhaps)], dtype = np.float32)
            for i,indhap in enumerate(haps):
                sp = indhap.split("|")
                if (len(sp) == 1):
                    if indhap in ["0","1","2"]:
                        pvec[i] = math.log(emslice[int(indhap)])
                    else:
                        pvec[i] = math.log(em0*(p0**2) + em2*(p1**2) + em1*2*(p0*p1))
                else:
                    if (sp[0] in genposs) and (sp[1] in genposs):
                        pvec[i] = math.log(emslice[int(sp[0])+int(sp[1])])
                    elif (sp[0]=="0") or (sp[1]=="0"):
                        pvec[i] = math.log(em0*p0 + em1*p1)
                    elif (sp[0]=="1") or (sp[1]=="1"):
                        pvec[i] = math.log(em1*p0 + em2*p1)
                    else:
                        pvec[i] = math.log(em0*(p0**2) + em2*(p1**2) + em1*2*(p0*p1))



            if si == 0:
                reflog = -2*math.log(lhaps)
                pvec += reflog

            else:

                pvec += prevpvec

                btflat = [str(item) for item in range(lhaps)]

                str1 = "SNP_"+str(si)+"\n"
                btfile.write(str1.encode('utf-8'))
                str1 = ",".join(btflat) + "\n"
                btfile.write(str1.encode('utf-8'))

            prevpvec = pvec



    #Termination step
    pmat = np.array(pvec)
    termp = np.amax(pmat)
    termbt = np.argwhere(pmat >= termp+tolerance*termp)

    termbtcollapse = ";".join(np.ndarray.flatten(termbt).astype(str))

    infile.close()
    inrecomb.close()
    btfile.close()
    sum_p = termp + np.log(np.sum(np.exp(-termp + pmat)))
    outfile = open(prefix+"_"+str(chromID)+"_Probability_value.txt",'w')
    outfile.write("Log-Probability value of best-fit trajectories = "+str(termp))
    outfile.write("log(Joint Probability of SNPs) - log(Product of Independent HWE Probabilities of SNP Genotypes) = "+str(sum_p - indep_p))
    outfile.close()
    return btfilename, termbtcollapse, lhaps, samplelist, snplist, chromID
#*****************************************************************
#*****************************************************************
def run_hmm_iterative(prefix, dsfilename, chromID, refpop, recombfile, scaledmutrate, tolerance, nsnps, obsgtypelist, mutratelist, numproc, en_select_sample):
    indep_p = 0
    select_index = en_select_sample[0]
    select_sample = en_select_sample[1]

    emissionmat = emissionprob(scaledmutrate)
    lowerlimit = sys.float_info.min
    genposs = ["0","1"]
    #EXTRACT the haplotypes of all the individuals from 1kGenomes
    btfilename = prefix+"_"+str(select_index)+"_"+chromID+"_BackTrace_Sequences.txt.gz"
    btfile = gzip.open(btfilename,'wb')
    infile = open(dsfilename,'r')
    inrecomb = open(recombfile, 'r')
    mm = mmap.mmap(infile.fileno(), 0, access = mmap.ACCESS_READ)
    mm.seek(0)


    for si in range(-1,nsnps,1):
        flag = 0
        terms = mm.readline().strip().split(b"\t")
        terms = [bitem.decode('utf-8') for bitem in terms]
        if(terms[0][:6]=="#CHROM"):

            samplelistprime = terms[9:]
            ###Clean-up for specific pathology in Sample IDs for 1kGenomes vcf files
            sampletemp = []
            for item in samplelistprime:
                if len(item)==14:
                    sampletemp.append(item[:7])
                    sampletemp.append(item[7:])
                else:
                    sampletemp.append(item)

            samplelisttot = [sampletemp[int(i/2)]+"_A" if i%2 ==0 else sampletemp[int((i-1)/2)]+"_B" for i in range(2*len(sampletemp))][:2*refpop]

            samplelist = [samplelisttot[ind] for ind in select_sample]
            indexlist = [(i,j) for i in range(len(samplelist)) for j in range(0,i+1,1)]
        elif terms[0][0]!="#":

            if (si>0):
                localrate = float(inrecomb.readline())
                transmat = transprob_unbiased(2*refpop,localrate)

                transmat = [math.log(lowerlimit) if item == 0 else math.log(item) for item in transmat]
                ton = transmat[0]
                toff = transmat[1]



            obsgtype = obsgtypelist[si]
            if mutratelist == None:
                emslice = [lowerlimit if item==0 else item for item in emissionmat[int(obsgtype)]]
                em0 = emslice[0]
                em1 = emslice[1]
                em2 = emslice[2]
            else:
                mutrate = mutratelist[si]
                mutation_tuples = [(int(item.split(":")[0]),float(item.split(":")[1])) for item in mutrate]
                sumprob = sum([tup[1] for tup in mutation_tuples])
                mutation_tuples.append((int(obsgtype), 1-sumprob))
                mutation_tuples.sort(key=lambda x: x[0])
                emslice = [item[1] for item in mutation_tuples]
                emslice = [lowerlimit if item==0 else item for item in emslice]
                em0 = emslice[0]
                em1 = emslice[1]
                em2 = emslice[2]

            freqs = terms[7].split(";")

            allelefreqprime = [freq.split("=")[1] for freq in freqs if (freq[:3]=="AF=")][0]
            if len(allelefreqprime.split(","))>1:
                allelefreq = float(allelefreqprime.split(",")[0])
            else:
                allelefreq = float(allelefreqprime)
            p1 = allelefreq
            p0 = 1.0-p1
            haps = terms[9:refpop+9]
            haps = [indhap+"|"+indhap if len(indhap.split("|")) == 1 else indhap for indhap in haps]
            haps = [item2 for indhap in haps for item in indhap.split(":")[0].split("|") for item2 in item.split("/")]

            haps = [haps[ind] for ind in select_sample]
            if int(obsgtype) == 2:
                SNP_prob = np.log(p1*p1)
            elif int(obsgtype) == 1:
                SNP_prob = np.log(2*p1*p0)
            elif int(obsgtype) == 0:
                SNP_prob = np.log(p0*p0)
            indep_p += SNP_prob
            lhaps = len(haps)

            #Explicit calculation of emission probabilities for all possibilities
            haps0 = np.array([i for i in range(lhaps) if haps[i]=="0"], dtype = np.uint32)
            haps1 = np.array([i for i in range(lhaps) if haps[i]=="1"], dtype = np.uint32)
            hapsm = np.array([i for i in range(lhaps) if (haps[i] != "0") and (haps[i] != "1")], dtype = np.uint32)
            if len(hapsm) == 0:
                flag = 1


            if flag != 1:
                emmm = em0*(p0**2) + em2*(p1**2) + em1*2*(p0*p1)
                emm0 = em0*p0 + em1*p1
                emm1 = em1*p0 + em2*p1
            pvec = np.array([math.log(lowerlimit) for i in range(lhaps) for j in range(0,i+1,1)], dtype = np.float32)

            for i in haps0:
                for j in haps0:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(em0)
                for j in haps1:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(em1)
                for j in hapsm:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emm0)
            for i in haps1:
                for j in haps0:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(em1)
                for j in haps1:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(em2)
                for j in hapsm:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emm1)
            for i in hapsm:
                for j in haps0:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emm0)
                for j in haps1:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emm1)
                for j in hapsm:
                    if(j<=i):
                        pvec[int(i*(i+1)/2+j)] = math.log(emmm)




            if si == 0:
                reflog = -2*math.log(lhaps)
                pvec += reflog

            else:
                pmat = prevpvec + 2*toff

                delta = ton-toff
                pool = mp.Pool(processes=numproc)

                fix_update=partial(update_func, pmat, lhaps, delta, tolerance)

                returnlist = pool.map(fix_update, indexlist)
                pool.close()
                pool.join()

                totpvec = np.array([item[0] for item in returnlist],dtype = np.float32)
                btvec = [item[1] for item in returnlist]
                pvec += totpvec
                btflat = [";".join(np.ndarray.flatten(item).astype(str)) for item in btvec]

                str1 = "SNP_"+str(si)+"\n"
                btfile.write(str1.encode('utf-8'))
                str1 = ",".join(btflat) + "\n"
                btfile.write(str1.encode('utf-8'))

            prevpvec = pvec

    #Termination step
    pmat = np.array(pvec)
    termp = np.amax(pmat)
    termbt = np.argwhere(pmat >= termp+tolerance*termp)

    termbtcollapse = ";".join(np.ndarray.flatten(termbt).astype(str))

    infile.close()
    inrecomb.close()
    btfile.close()

    sum_p = termp + np.log(np.sum(np.exp(-termp + pmat)))

    return sum_p, indep_p, termp, btfilename, termbtcollapse, lhaps, samplelist, samplelisttot, chromID

#*****************************************************************
#*****************************************************************
def iterative_process(prefix, dsfilename, filename, observedsample, refpop, subgroup, chromID, recombfile, scaledmutrate, tolerance, numproc, niter, posspecific):
    if posspecific=="True":
        (samplelist,obsgtypelist,mutratelist,nsnps,snplist) = process_SNPs_position_specific(dsfilename, filename, observedsample, refpop)
    else:
        (samplelist,obsgtypelist,nsnps,snplist) = process_SNPs(dsfilename, filename, observedsample, refpop)
        mutratelist = None
    writeflag = "init"
    count = 0
    bestlist = samplelist[:2*refpop]
    prevbest = samplelist
    while len(bestlist)>subgroup:
        toplist = []
        for iter in range(niter):
            if count == 0:
                select_sample = selection_func(samplelist, subgroup)
            else:
                random.shuffle(bestlist)
                select_sample = [bestlist[i*subgroup:(i+1)*subgroup] for i in range((len(bestlist)+subgroup-1) // subgroup)]
            for en_select_sample in enumerate(select_sample):
                sum_p, indep_p, termp, btfilename, termbtcollapse, lhaps, samplelist, samplelisttot, chromID = run_hmm_iterative(prefix, dsfilename, chromID, refpop, recombfile, scaledmutrate, tolerance, nsnps, obsgtypelist, mutratelist, numproc, en_select_sample)
                returnlist = read_backtrace_iterative(prefix, btfilename, termbtcollapse, lhaps, samplelist, samplelisttot, snplist, chromID, writeflag)
                toplist.extend(returnlist)
        bestlist = list(set(toplist))
        if (set(prevbest) == set(bestlist)) or (len(bestlist) > len(prevbest)):
            break
        prevbest = bestlist
        count += 1
    writeflag = "final"
    sum_p, indep_p, termp, btfilename, termbtcollapse, lhaps, samplelist, samplelisttot, chromID = run_hmm_iterative(prefix, dsfilename, chromID, refpop, recombfile, scaledmutrate, tolerance, nsnps, obsgtypelist, mutratelist, numproc, (0,bestlist))
    read_backtrace_iterative(prefix, btfilename, termbtcollapse, lhaps, samplelist, samplelisttot, snplist, chromID, writeflag)

    outfile = open(prefix+"_"+str(chromID)+"_Probability_value.txt",'w')
    outfile.write("Log-Probability value of best-fit trajectories = "+str(termp))
    outfile.write("log(Joint Probability of SNPs) - log(Product of Independent HWE Probabilities of SNP Genotypes) = "+str(sum_p - indep_p))
    outfile.close()
    return

#*****************************************************************
#*****************************************************************
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Remove most informative SNP and rerun PLIGHt algorithm of choice')
    parser.add_argument('--bftfile', required=True, help='Best-fit trajectory file from previous run of PLIGHT')
    parser.add_argument('--algorithm_type', required=True, choices = ['Exact','Iterative','InReference'], help='Version of PLIGHT to run the post-sanitization protocol on (Truncated option not included)')
    parser.add_argument('-c','--chromosomefile', required=True, help='Chromosome file name')
    parser.add_argument('-O','--observedsample', required=True, help='Observed Sample Genotype File')
    parser.add_argument('-I','--chromosomeID', required=False, help='Chromosome ID')
    parser.add_argument('-m','--metadata', required=False, help='Metadata file with ancestry information', default='integrated_call_samples_v3.20130502.ALL.panel')
    parser.add_argument('-F','--genfolder', required=False, help='Genotype folder', default='Genotypes/')
    parser.add_argument('-M','--thetamutationrate', required=False, type = float, help='Theta (Coalescent) Mutation rate')
    parser.add_argument('-L','--lambdamutationrate', required=False, type = float, help='Lambda (direct error) Mutation rate')
    parser.add_argument('-e','--effpop', required=False, type = int, help='Effective population', default=11418)
    parser.add_argument('-r','--refpop', required=False, type = int, help='Reference population', default=2504)
    parser.add_argument('-b','--recombrate', required=False, type = float, help='Recombination rate in cM/Mb (if /distance/ option is chosen)', default=0.5)
    parser.add_argument('-s','--recombswitch', required=False, choices = ['distance', 'custom'], help='Recombination Model Switch (distance = simple linear increase of recombination rate with distance, custom  = alternative model)', default='distance')
    parser.add_argument('-f','--recombfile', required=False, help='Custom Recombination Model File')
    parser.add_argument('-t','--tolerance', required=False, type = float, help='Fraction of maximum value used as allowance for inclusion of a path', default=0.01)
    parser.add_argument('-C','--currdir', required=False, help='Current working directory', default='./')
    parser.add_argument('--numproc', required=False, type = int, help='Number of processes for parallelization', default=1)
    parser.add_argument('--subgroup', required=False, help='Number of haplotypes in each iterative scheme subgroup', default=200, type=int)
    parser.add_argument('--niter', required=False, help='Number of iterations of bootstrapping process', default=20, type=int)
    parser.add_argument('--posspecific', required=False, help='Position-specific mutation rates included in observation file? (True/False)', default="False")
    parser.add_argument('--prefix', required=False, help='String prefix to append to output Best_trajectories file, in addition to chromosome number', default="")
    args = parser.parse_args()

    bftfile = args.bftfile
    algorithm_type = args.algorithm_type
    chromfile = args.chromosomefile
    chromID = args.chromosomeID
    metadata = args.metadata
    genfolder = args.genfolder
    currdir = args.currdir
    tolerance = args.tolerance
    effpop = args.effpop
    refpop = args.refpop
    numproc = args.numproc
    recombrate = 0.01*args.recombrate
    filename = genfolder + chromfile

    observedsample = args.observedsample
    subgroup = args.subgroup
    niter = args.niter
    posspecific = args.posspecific
    prefix = args.prefix
    dsfilename = prefix+"_"+chromID+"_SNP_Haplotypes.txt"
    if (args.thetamutationrate is None) and (args.lambdamutationrate is None):
        mutrate = (sum(1.0/i for i in range(1,2*refpop-1)))**(-1)
        scaledmutrate = 0.5*mutrate/(mutrate + 2*refpop)
    elif (args.thetamutationrate is None):
        scaledmutrate = args.lambdamutationrate
    elif (args.lambdamutationrate is None):
        mutrate = args.thetamutationrate
        scaledmutrate = 0.5*mutrate/(mutrate + 2*refpop)

    #*****************************************************************
    #*****************************************************************
    #Read in the current Best-fit trajectories file and figure out the SNPs with the maximum leakage
    infile = open(bftfile,'r')

    count = 0
    snplist = []
    haplength = []
    for line in infile:
        terms = line.strip().split("\t")
        SNP_ID = terms[0]
        snplist.append(SNP_ID)
        gens = terms[1:]
        if count==0:
            haps = list(set([gen.split("_")[i]+"_"+gen.split("_")[i+1] for gen in gens for i in range(0,3,2)]))
            haplength.append(len(haps))
        else:
            haps = list(set([nextgen.split("_")[i]+"_"+nextgen.split("_")[i+1] for gen in gens for nextgen in gen.split(":")[1].split(";") for i in range(0,3,2)]))
            haplength.append(len(haps))

        count += 1
    infile.close()
    min_hl = np.min(np.array(haplength))
    min_SNPs = np.array(snplist)[np.where(np.array(haplength) == min_hl)]
    print(min_SNPs)
    #For the set of maximally leaking SNPs, remove the one with the lowest MAF
    af_file = f"{chromID}_Observed_AlleleFreq.txt"
    bcfcommand = f"bcftools view -m2 -M2 -v snps -Ov -R {observedsample} {filename} | bcftools query -f '%CHROM\t%POS\t%AF\n' > {af_file}"
    subprocess.call(bcfcommand,shell=True)

    infile = open(af_file,'r')
    af_dict = {}
    for line in infile:
        terms = line.strip().split("\t")
        SNP_ID = f"chr{terms[0].strip('chr')}_{terms[1]}"
        AF = float(terms[2])
        af_dict[SNP_ID] = AF
    infile.close()
    min_AFs = np.array([af_dict["_".join(SNP.split("_")[:-1])] for SNP in min_SNPs])
    most_informative_SNP = min_SNPs[np.where(min_AFs == np.min(min_AFs))][0]
    print(most_informative_SNP)

    newname = ".".join(observedsample.split(".")[:-1])+".Sanitized.txt"
    infile = open(observedsample,'r')
    outfile = open(newname,'w')
    for line in infile:
        terms = line.strip().split("\t")
        SNP_ID = "chr"+terms[0].strip("chr")+"_"+terms[1]
        if SNP_ID not in most_informative_SNP:
            outfile.write(line)
    outfile.close()
    infile.close()
    observedsample = newname
    #*****************************************************************
    #*****************************************************************
    recombswitch = args.recombswitch
    if recombswitch == 'custom':
        recombfile = args.recombfile
    else:
        recombfile = generate_recomb(prefix, effpop, chromID, recombrate, observedsample)

    if algorithm_type == 'Exact':
        btfilename, termbtcollapse, lhaps, samplelist, snplist, chromID = run_hmm_exact(prefix, dsfilename, filename, observedsample, chromID, refpop, recombfile, scaledmutrate, tolerance, numproc, posspecific)
        trajfile = read_backtrace_exact(prefix, btfilename, termbtcollapse, lhaps, samplelist, snplist, chromID)
    elif algorithm_type == 'Iterative':
        iterative_process(prefix, dsfilename, filename, observedsample, refpop, subgroup, chromID, recombfile, scaledmutrate, tolerance, numproc, niter, posspecific)
    elif algorithm_type == 'InReference':
        btfilename, termbtcollapse, lhaps, samplelist, snplist, chromID = run_hmm_in_reference(prefix, dsfilename, filename, observedsample, chromID, refpop, recombfile, scaledmutrate, tolerance, numproc, posspecific)
        trajfile = read_backtrace_in_reference(prefix, btfilename, termbtcollapse, lhaps, samplelist, snplist, chromID)
