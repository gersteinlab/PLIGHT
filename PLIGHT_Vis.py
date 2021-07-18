import sys
import os
import subprocess
import argparse
import numpy as np
import collections
from collections import Counter
import matplotlib.pyplot as plt
class Node:
    def __init__(self,name = None,lyr = None,*links):
        self.name = name
        self.layer = lyr

        for idx, item in enumerate(links):
            setattr(self, "link{}".format(idx), item)
    def getName(self):
        return self.name
    def getLayer(self):
        return self.layer
    def getLinks(self):
        return [getattr(self,attr) for attr in dir(self) if attr.startswith("link")]
    def setname(self,name):
        self.name = name
    def setLayer(self,lyr):
        self.layer = lyr
    def setLinks(self, *linkvals):
        for idx, item in enumerate(linkvals):
            setattr(self, "link{}".format(idx), item)
class Trajectory:
    def __init__(self):
        self.head = [Node("Tail",0,[None])]
    def addLayer(self, genotypes, lyr, membership):

        temp = [Node(gen, lyr, [prevnode for j,prevnode in enumerate(self.head) if prevnode.getName() in membership[gen] ] ) for gen in genotypes]

        self.head = temp
def generate_trajset(trajfolder, chromID, trajfile):
    infile = open(trajfile)

    count = 0
    for line in infile:
        nodes = line.strip().split("\t")[1:]

        if count==0:
            trajlist = [[node] for node in nodes]

        else:
            nodedict = dict([(node.split(":")[0], node.split(":")[1].split(";")) for node in nodes])

            trajlist = [traj+[item] for traj in trajlist for item in nodedict[traj[-1]]]

        count += 1
    infile.close()
    outfile = open(trajfolder+"/"+chromID+"_Full_trajectory_set.csv","w")
    outfile.writelines(",".join(trajlist[i]) + "\n" for i in range(len(trajlist)))
    outfile.close()
    return(trajlist)
def optimal_haplotype(filename):
    infile = open(filename)
    totgtypes = [line.strip().split(",") for line in infile]
    infile.close()
    max_haps = []
    totcounts = Counter()

    for snp in range(len(totgtypes[0])):
        haps = [[totgtypes[j][snp].split("_")[0]+"_"+totgtypes[j][snp].split("_")[1], totgtypes[j][snp].split("_")[2]+"_"+totgtypes[j][snp].split("_")[3]] for j in range(len(totgtypes))]

        haps = [hap for gen in haps for hap in gen]

        counts = Counter(haps).most_common()
        totcounts = totcounts + Counter(haps)
        max_hap = [item[0] for item in counts if item[1] == counts[0][1]]
        max_haps.append(max_hap)

    hapdict = {}
    for trajidx in range(len(totgtypes)):
        haps = [[totgtypes[trajidx][j].split("_")[0]+"_"+totgtypes[trajidx][j].split("_")[1], totgtypes[trajidx][j].split("_")[2]+"_"+totgtypes[trajidx][j].split("_")[3]] for j in range(len(totgtypes[trajidx]))]

        haps = list(set([hap for gen in haps for hap in gen]))
        for key in haps:
            if key not in hapdict.keys():
                hapdict[key] = 1.0/len(totgtypes)
            else:
                temp = hapdict[key]
                hapdict[key] = temp + 1.0/len(totgtypes)

    countmax = []
    for i,item in enumerate(max_haps):
        if i == 0:
            curr = item
            count = 1
        else:
            if item == curr:
                count += 1
            else:
                countmax.append((curr,count))
                curr = item
                count = 1
    countmax.append((curr,count))
    countmax.reverse()
    return [max_haps,totcounts,countmax,hapdict]

def cross_chromosome_optimum(trajfolder,chromosomeIDs):

    totcounts = Counter()
    files = []
    totdict = {}
    outfile = open("Optimal_haplotype_at_each_locus.txt",'w')
    for chrom in chromosomeIDs:
        filename = trajfolder+"/"+chrom+"_Full_trajectory_set.csv"
        files.append(filename)
        [max_haps,chromcounts,countmax,hapdict] = optimal_haplotype(filename)
        outfile.write(chrom+":\n")
        outfile.write("\t".join([",".join(item) for item in max_haps])+"\n")
        totcounts = totcounts + chromcounts
        for key in hapdict.keys():
            if key not in totdict.keys():
                totdict[key] = hapdict[key]/len(chromosomeIDs)
            else:
                temp = totdict[key]
                totdict[key] = temp + hapdict[key]/len(chromosomeIDs)

    outfile.close()

    outfile = open(f'{trajfolder}/Consensus_trajectories_each_chromosome.csv','w')
    for filename in files:
        infile = open(filename)
        chromID = os.path.basename(filename).split("_")[0]
        totgtypes = [line.strip().split(",") for line in infile]
        infile.close()
        scorearray = np.zeros(len(totgtypes),dtype = np.int16)
        for trajidx in range(len(totgtypes)):
            haps = [[totgtypes[trajidx][j].split("_")[0]+"_"+totgtypes[trajidx][j].split("_")[1], totgtypes[trajidx][j].split("_")[2]+"_"+totgtypes[trajidx][j].split("_")[3]] for j in range(len(totgtypes[trajidx]))]
            haps = [hap for gen in haps for hap in gen]
            counts = Counter(haps)
            countdict = dict(counts)
            score = sum([totdict[key] for key in countdict.keys()])
            scorearray[trajidx] = score

        maxval = np.max(scorearray)
        maxargs = np.argwhere(scorearray == maxval)
        outfile.writelines(chromID+","+",".join(totgtypes[item[0]])+"\n" for item in maxargs)
    outfile.close()

def plot_traj(plotfolder, chromID, trajfile):
    infile = open(trajfile)
    count = 0
    traj = Trajectory()
    genotypes = []
    nodelist = []
    for line in infile:
        nodes = line.strip().split("\t")[1:]
        print(nodes)
        if count==0:
            membership = {}
            for node in nodes:
                membership[node] = 'Tail'
            traj.addLayer(nodes,count,membership)
            nodelist.append(traj.head)
            genotypes.append([item for item in nodes])
        else:
            nextnodes = list(set([item for node in nodes for item in node.split(":")[1].split(";")]))
            genotypes.append([item for item in nextnodes])
            prevnodes = [node.split(":")[0] for node in nodes]
            membership = {}
            for nextnode in nextnodes:
                membership[nextnode] = [node.split(":")[0] for node in nodes if nextnode in node.split(":")[1].split(";")]

            traj.addLayer(nextnodes,count,membership)
            nodelist.append(traj.head)

        print(count)
        count += 1
    infile.close()
    map_gtype = []
    x_gtype = []
    nlyrnodes = []
    nodelist.reverse()
    genotypes.reverse()
    scalefactor = 2
    for lyr in range(len(genotypes)):
        nnodes = len(genotypes[lyr])
        print(nnodes)
        tag_gtype = [item+":"+str(lyr) for item in genotypes[lyr]]
        tag_gtype.sort()
        genotypes[lyr] = tag_gtype
        nlyrnodes.append(nnodes)
        lower = int(round(-(scalefactor-1)*nnodes/2.0 + 0.5))
        map_gtype.append(range(lower, lower+scalefactor*nnodes, scalefactor))
        x_gtype.append([scalefactor*lyr+1 for i in range(nnodes)])
    x_gtype = [item for row in x_gtype for item in row]
    map_gtype = [item for row in map_gtype for item in row]
    genotypes = [item for row in genotypes for item in row]

    xmapdict = dict(zip(genotypes,x_gtype))
    ymapdict = dict(zip(genotypes,map_gtype))


    for cindex, currlist in enumerate(nodelist[:-1]):

        for node in currlist:
            currname = node.getName()+":"+str(cindex)
            links = node.getLinks()
            nextnames = [item.getName()+":"+str(cindex+1) for item in links[0]]

            xvalcurr = xmapdict[currname]

            yvalcurr = ymapdict[currname]

            xvalnext = [xmapdict[nnames] for nnames in nextnames]

            yvalnext = [ymapdict[nnames] for nnames in nextnames]

            for i in range(len(nextnames)):
                plt.arrow(xvalcurr, yvalcurr, xvalnext[i]-xvalcurr,yvalnext[i]-yvalcurr, width = 0.01,length_includes_head = True, color = 'b')

    maxnodes = max(nlyrnodes)

    [max_haps,totcounts,countmax,hapdict] = optimal_haplotype(trajfolder+"/"+chromID+"_Full_trajectory_set.csv")
    xcoord = [min(x_gtype)]
    currval = 0
    for item in countmax:
        nextval = item[1]+currval
        xcoord.append(scalefactor*nextval)
        currval = nextval

    plt.plot(xcoord,[scalefactor*round(0.5*maxnodes)+scalefactor-1 ]*len(xcoord),'g|', markersize = 10, mew = 1)
    for i in range(len(xcoord)-1):
        plt.arrow(xcoord[i], scalefactor*round(0.5*maxnodes)+scalefactor-1 , (xcoord[i+1]-xcoord[i]),0, width = 0.01,length_includes_head = True, color = 'g')
        plt.annotate(countmax[i][0], xy=(round(0.5*(xcoord[i+1]+xcoord[i])),scalefactor*round(0.5*maxnodes)+scalefactor-1 ), xytext=(-20,10), fontsize=8,textcoords='offset points', rotation=45,va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='lightblue', alpha=0.5))
    plt.plot(x_gtype,map_gtype,'ro')
    plt.axis([-2, scalefactor*count+scalefactor+1, -2-scalefactor*round(0.5*maxnodes), scalefactor*round(0.5*maxnodes) + 2*scalefactor + 1])
    topgenotypes = [item.split(":")[0].split("_")[0]+"_"+item.split(":")[0].split("_")[1] for item in genotypes]
    botgenotypes = [item.split(":")[0].split("_")[2]+"_"+item.split(":")[0].split("_")[3] for item in genotypes]
    for toplabel, x, y in zip(topgenotypes, x_gtype, map_gtype):
        plt.annotate(toplabel, xy=(x+scalefactor,y), xytext=(-20,5), fontsize=4,textcoords='offset points',rotation=45, va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5))
    for botlabel, x, y in zip(botgenotypes, x_gtype, map_gtype):
        plt.annotate(botlabel, xy=(x+scalefactor,y), xytext=(-20,-5), fontsize=4,textcoords='offset points', rotation=-45, va='top',
        bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5))
    plt.gca().axes.yaxis.set_ticklabels([])
    plt.gca().axes.xaxis.set_ticklabels([])
    plt.savefig(plotfolder+"/"+chromID+"_plot.png", dpi=300)
    #plt.show()
def consensus_plotting(chromosomeIDs, trajfolder, plotfolder, consensusfile):
    consensusdict = {}
    snplength = {}
    with open(consensusfile,'r') as infile:
        for line in infile:
            terms = line.strip().split(",")
            if terms[0] not in consensusdict.keys():
                consensusdict[terms[0]] = [terms[1:]]
                snplength[terms[0]] = len(terms[1:])
            else:
                consensusdict[terms[0]].append(terms[1:])

    for key,val in consensusdict.items():
        outfile = open(f"{trajfolder}/{key}_Consensus_Best_trajectories.tsv",'w')
        chromarr = np.array(val)
        startgen = "\t".join(list(set(chromarr[:,0])))
        snplength_chrom = snplength[key]
        outfile.write(f"SNP_{snplength_chrom}\t{startgen}\n")
        for SNP_ID in range(snplength_chrom-1):
            listSNP = [tuple(item) for item in chromarr[:,SNP_ID:SNP_ID+2]]
            outfile.write(f"SNP_{snplength_chrom-SNP_ID-1}")
            trajdict = {}
            for traj in listSNP:
                if traj[0] not in trajdict.keys():
                    trajdict[traj[0]] = [traj[1]]
                else:
                    trajdict[traj[0]].append(traj[1])
            for key, val in trajdict.items():
                valstring = ";".join(set(val))
                outfile.write(f"\t{key}:{valstring}")
            outfile.write("\n")
        outfile.close()

    for chrom in chromosomeIDs:
        trajfile = f"{trajfolder}/{chrom}_Consensus_Best_trajectories.tsv"
        generate_trajset(trajfolder, chrom, trajfile)
        plot_traj(trajfolder, chrom, trajfile)
    return
if __name__=="__main__":
    parser = argparse.ArgumentParser(description='Unravel the trajectories, find the identities of the best-fit reference haplotypes, and plot the results')
    parser.add_argument('-C','--chromosomeIDs', required=False, nargs = '+', help='List of Chromosome IDs to be considered: format chr followed by number, listed singly, separated by spaces (for eg. -C chr1 chr2 chr19)' )
    parser.add_argument('-t','--trajectorypattern', required=False, help='Pattern for trajectory files, where the chromosome ID would fit into the curly braces', default="{}_Best_trajectories.tsv")
    parser.add_argument('-T','--trajectoryfolder', required=False, help='Folder for storing the trajectory files', default="Trajectories")
    parser.add_argument('-P','--plotfolder', required=False, help='Folder for storing the plots', default="HMM_Plots")

    args = parser.parse_args()

    chromosomeIDs = args.chromosomeIDs
    trajectorypattern = args.trajectorypattern
    if(chromosomeIDs == None):
        chromosomeIDs = ["chr"+str(item) for item in range(1,23)]

    trajfolder = args.trajectoryfolder
    plotfolder = "HMM_plots"
    if not os.path.exists(trajfolder):
        os.mkdir(trajfolder)
    if not os.path.exists(plotfolder):
        os.mkdir(plotfolder)

    if trajfolder[-1] == "/":
        trajfolder = trajfolder[:-1]

    tottraj = []
    for chrom in chromosomeIDs:
        trajfile = trajfolder+"/"+trajectorypattern.format(chrom)
        ####This function generates a comprehensive set of independent best-fit trajectories, by enumerating all the trajectories written to the "Best_trajectories" file. For each chromosome
        #### these files are stored with a "Full_trajectory_set.csv" suffix
        generate_trajset(trajfolder, chrom, trajfile)

        ####This function builds a framework of linked nodes for each best-fit trajectory, and then plots all the trajectories in a png file saved to the plotfolder location.
        #### Additionally, the function calls the "optimal_haplotype" function, which outputs several metrics: a list of most frequent haplotypes at each locus; a Python Counter object of
        #### the haplotypes present int he trajectories; a running list of haplotype counts including the entire SNP interval over which that haplotype is the most frequent; and a dictionary
        #### containing the fraction of the number of trajectories over which a particular haplotype is found.
        plot_traj(plotfolder, chrom, trajfile)

    #### This function looks at consensus trajectories (as described in the paper) that are found by weighting trajectories for each chromosome according to the presence of certain
    #### high frequency haplotypes. These consensus trajectories are written to a file called "Consensus_trajectories_each_chromosome.csv", which, as the name suggests, has the best consensus
    #### trajectories for each chromosome.
    cross_chromosome_optimum(trajfolder,chromosomeIDs)

    consensusfile = f"{trajfolder}/Consensus_trajectories_each_chromosome.csv"
    #### This function runs the same trajectory plotting scheme, but on the consensus trajectories calculated above.
    consensus_plotting(chromosomeIDs, trajfolder, plotfolder, consensusfile)
