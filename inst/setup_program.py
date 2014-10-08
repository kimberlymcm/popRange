import numpy as np
import sys
import math
import random as rn
import config
import popRange_main


##########################################
##Initial set up of the world           ##
##########################################

def setup():
    #Making the lattice
    # lattice = np.empty( (makePopulations.max_X, makePopulations.max_Y), dtype=object)
    lattice = np.empty( (config.max_X, config.max_Y), dtype=object)
    for i in range(0,len(config.populations)):
        x = config.populations[i][0]
        y = config.populations[i][1]
        lattice[x, y] = popRange_main.site()
    
    if config.inFile_YEP != None:
        setupFromInputFile( lattice, config.inFile_YEP )
        #The parameters set in the input are "allele freqs", "genArose", and "s"
        if config.s != None:
            lattice[x,y].snps.sel = startingSel(lattice[x,y])
        elif config.sDiff != None:
            #print "hihi"
            lattice = spatialSel(lattice)
    
    #theta = 2 * config.gSize * config.mutRate * config.Ne #for 50 million bp, Ne=100, mutRate=2*10^-8, theta is 200
    #In this case, I am getting the snm for all the populations at once, and then dividing it
    elif config.SNP_model == 1: #1 is the stanford neutral model
        
        lattice = startingNe( lattice )
        SNP_mat = setSNM( lattice )
        
        prevNe = 0
        for i in range(0,len(config.populations)):
            x = config.populations[i][0]
            y = config.populations[i][1]
            if lattice[x,y].ne > 0: #This is set in 'starting Ne'
                lattice[x,y].snps.freq = np.sum(SNP_mat[prevNe:(lattice[x,y].ne+prevNe)], axis=0)
                prevNe = lattice[x,y].ne + prevNe
                config.nSNPs = len(lattice[x,y].snps.freq)
                config.snpNum = config.nSNPs - 1
                if config.sDiff == None:
                    lattice[x,y].snps.sel = startingSel(lattice[x,y])
                elif i == 0: #You only do this once, since it goes through all the populations!
                    lattice = spatialSel(lattice)
                qq = np.arange(0, config.nSNPs).tolist()
                lattice[x,y].snps.pos = list(qq)
                lattice[x,y].snps.arose = [0] * int(config.nSNPs)
    
    
    elif config.SNP_model == 0: #SNP model must be 0 and it is fixed
        config.snpNum = int(config.nSNPs) - 1 #config.nSNPs is already set, in this case, by the user
        lattice = startingNe(lattice)
        for i in range(0, len(config.populations)):
            x = config.populations[i][0]
            y = config.populations[i][1]
            if lattice[x,y].ne > 0:
                lattice[x,y].snps.freq = [ round(float(config.SNPs_starting_freq)*lattice[x,y].ne) ] * int(config.nSNPs)
                if config.sDiff == None:
                    lattice[x,y].snps.sel = startingSel(lattice[x,y])
                elif i == 0:
                    lattice = spatialSel(lattice)
                lattice[x,y].snps.arose = [0] * int(config.nSNPs)
                qq = np.arange(0, int(config.nSNPs) ).tolist()
                lattice[x,y].snps.pos = list(qq)
    else:
        sys.exit("ERROR: SNP_model must be either 0 (snm) or 1 (fixed)")
    
    createAlleleHistory(lattice)
    return lattice



def createAlleleHistory(lattice):
    global allele_history
    x = config.populations[0][0] #arbitrarily picked the first population to get the details from
    y = config.populations[0][1]
    
    a = len(config.populations)*2 + 4 #this is the # of columns (rows for now) I will need
    
    temp_list = np.zeros( (a,len(lattice[x,y].snps.freq)) ) #a rows, len(lattice.snps.freq) columns
    temp_list[1] = list(lattice[x,y].snps.pos) #SNP nums for everything.  Hopefully a deep copy
    temp_list[2] = list(lattice[x,y].snps.arose)
    
    #print len(lattice[x,y].snps.pos)
    #print len(lattice[x,y].snps.arose)
    #print len(lattice[x,y].snps.sel)
    #print lattice[x,y].snps.sel
    temp_list[3] = list(lattice[x,y].snps.sel) #SNP sel for everything.  Hopefully a deep copy.
    k = 4
    for i in range(0, len(config.populations)):
        x = config.populations[i][0]
        y = config.populations[i][1]
        if len(lattice[x,y].snps.freq) > 0:
            temp_list[k] = list(lattice[x,y].snps.freq)
            temp_list[k+1] = np.tile( lattice[x,y].ne, len(lattice[x,y].snps.freq) )
        k = k + 2
    
    allele_history = temp_list.transpose()
    #print allele_history


def startingNe (lattice):
    for i in range(0, len(config.populations)):
        x = config.populations[i][0]
        y = config.populations[i][1]
        if isinstance(config.Ne, int):
            if config.popMat[x][y] != -1:
                lattice[x,y].ne = config.Ne
                lattice[x,y].exists = True
        else:
            if config.Ne[x][y] > 0: #If it is given a pop size, then it exists
                lattice[x,y].exists = True
                lattice[ x, y ].ne = config.Ne[x][y]
    return lattice


def setSNM (lattice):
    total = 0
    for i in range(0,len(config.populations)):
        x = config.populations[i][0]
        y = config.populations[i][1]
        if lattice[x,y].exists == True:
            total = total + lattice[x,y].ne
    theta = 2 * config.gSize * config.mutRate * total
    #  print 'theta ' + str(theta) #E(s) = theta(1+(1/2)+(1/3)+(1/4)+...+(1/(n-1)))
    numSingle = theta
    w = np.arange(1.0,float(total),1)
    ww = 1/ w
    www = numSingle * ww
    x = np.rint(www)
    A = []
    for i in range(0,len(x)):
        AA = A
        B = [(i+1)/float( total)] * x[i]
        A = np.hstack((AA,B))
    AAAA = A * total #Now this is the model with all the SNPs and freqs
    #Therefore, I want to divide this into the amount of subsets by the amount of pops
    SNP_mat = np.zeros( (len(AAAA), total) ) #1 row per snp, 1 column per person
    for i in range( len(AAAA) ):
        SNP_mat[i][0:AAAA[i]] = 1
    map(np.random.shuffle, SNP_mat) #This shuffles the values in each row separately
    SNP_mat = np.array(SNP_mat).transpose()
    if (len(SNP_mat) == 0):
        sys.exit("Error: This simulation can't start with 0 SNPs. Perhaps try a higher mutation rate.")
    return SNP_mat


##This functions sets the selection coefficients if they are specifically different for each one.

def startingSel(pop): #This I am using if the SNPs are a following a gamma distribution
    if config.s != None: #If there is no selection, config.s == 0.
        pop.snps.sel = [config.s] * config.nSNPs
    elif config.s_mat != None:
        for line in config.s_mat:
            pop.snps.sel[(int(line[0])-1):int(line[1])] = [float(line[2])] * (int(line[1]) - int(line[0]) + 1)  #This is the number of snps. 1 based. includes both #s.
    else:
        pop.snps.sel = np.random.gamma( shape = config.g_a, scale = config.g_b, size=len(pop.snps.freq))
        aa = np.where(pop.snps.sel > 1) #restrict s to 0-1
        pop.snps.sel[aa] = 1
        bb = np.where(pop.snps.sel < 0)
        pop.snps.sel[bb] = 0
    return pop.snps.sel


######Starting selection if there is spatially variable selection
##BUT I THINK I ONLY WANT THIS TO BE AS LONG AS THE CURRENT SNPS
def spatialSel(lattice):
    cc = config.nSNPs
    for i in range(1,len(config.sDiff)):
        #sSNP, fSNP = config.sDiff[i][0:2] #The first and last SNP with the s, not using
        for k in range(0,len(config.populations)):
            x = config.populations[k][0]
            y = config.populations[k][1]
            popStr = str(x) + str(y)
            # print popStr + " " + str(config.sDiff[0][2])
            a = np.where(np.array(config.sDiff[0]) == popStr) #This is the index you want to use for the s
            cc_max = min(int(config.sDiff[i][1]), cc)
            #print cc_max
            #print int(config.sDiff[i][1])
            total = cc_max - (int(config.sDiff[i][0])-1)
            if int(config.sDiff[i][1]) > (cc + 1):
                print "User has indicated more SNPs (" + str(config.sDiff[i][1]) + ") in the sDiff file than in the nSNPs parameter (" + str(cc) + "). Please retry."
                sys.exit()
            #print total
            if len(a[0]) > 0:
                lattice[x,y].snps.sel[(int(config.sDiff[i][0])-1):cc_max] = np.tile( float(config.sDiff[i][a[0]]), total)
    return lattice


def setupFromInputFile( lattice, infile ):
    parsedFile = config.readFile( infile, remove=False )
    k = 1 #Now I'm reading in the pop info
    for i in range(0, len(config.populations)): #This is filling in all the info
        x = config.populations[i][0]
        y = config.populations[i][1]
        if parsedFile[k][0] == "None":
            pass
        else:
            lattice[x,y].ne = int(parsedFile[k][0]) #This should be the Ne
            lattice[x,y].snps.pos = map(int, parsedFile[k+1])
            lattice[x,y].snps.sel = map(float, parsedFile[k+2])
            lattice[x,y].snps.freq = map(int, parsedFile[k+3])
            lattice[x,y].exists = True
            lattice[x,y].snps.arose = map(int, parsedFile[k+4]) #idk if this should be like this!
        k = k+5
