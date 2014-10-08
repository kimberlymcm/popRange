#!/usr/bin/python

import numpy as np
import sys
import math
import random as rn
import config_3 #This imports my config_3 files with all my config_3 variables
import writeOutput_3 #This inputs the functions to write output in different ways
import os
import setup_program_3
import time


#########################################
##Defining initial variables of objects##
#########################################
#A "site" is a population.  Each population has a population size (self.ne), either exists or doesn't (self.exists),
#Has an allele tragectory for the selection alleles (self.alleleTrag), an allele tragectory for the neural alleles that I am keeping track of (self.alleleTragNeut),
#and has all of the currently segregating snps somewhere in the world (and their proportions in the population) in (self.snps).
#The (self.snps) class include the frequency of each snp in the population (self.freq),
#its selection coefficient (self.sel), and its position (self.pos)(what exactly does pos mean??)
class site:
    def __init__(self):
        self.ne = None
        self.exists = False
        self.snps = site.snps()
    class snps:
        def __init__(self):
            self.freq = []
            self.sel = []
            self.pos = []
            self.arose = []



################################################
## This is the main code the runs the program ##
################################################

def main():

    #start = time.time()

    config_3.readValues(sys.argv[1:])
    lattice = setup_program_3.setup()
    endLattice = sim(lattice) #Runs the simulations
    writeToFile(endLattice, config_3.nGens) #Write the final results to file

################################
## Now the generations happen ##
################################

def sim(lattice):
    c = 0
    for i in range(0, config_3.nGens):
        lattice = cat(lattice, i) #Sees if a catastrophe happens
        if len(config_3.populations) > 1 and config_3.mig == True: #Migration if more than 1 population
            lattice = mig(lattice, i)
        for j in range(0,len(config_3.populations)):
            x = config_3.populations[j][0]
            y = config_3.populations[j][1]

            if lattice[x,y].exists == True and lattice[x,y].ne > 0:
                if config_3.mutRate != None and config_3.mutRate != 0:
                    lattice = mut(lattice, x, y, i) #Does the mutations for each population that currently exists
                lattice = drift(lattice, x, y, i) #Does the drift for each population that currently exists
        if config_3.SNP_model == 1:
            lattice = loseZeroes(lattice) #Gets rid of the SNPs that are 0 in all populations, because these are unnecessary and take up extra time
        if config_3.recordTrag != 0 and ( i % config_3.recordTrag == 0):
            updateAlleleHistory(lattice, i)
    return lattice #Returns the lattice at the end of all of the generations

#################
## Catastrophe ##
#################

def cat(lattice, gen):
    #This is the probability of catastrophe in any grid point
    #(including grid points that don't currently have any individuals) during this generation
    for i in range(0,len(config_3.populations)):
        x = config_3.populations[i][0]
        y = config_3.populations[i][1]
        cat = 0
        if isinstance(config_3.catProb, (int,float,complex)):
            if config_3.catProb > 0:
                cat = np.random.binomial(1,config_3.catProb)
        elif config_3.catProb[x][y] > 0:
            cat = np.random.binomial(1,config_3.catProb[x][y])
        if cat == 1:
            print('Catastrophe!')
            if lattice[x,y].exists == True:
                lattice = extinctPop(lattice,x,y,gen)
            if isWorldStillAlive(lattice,gen) == False:
                worldEnded()
    return lattice

################
## Counting to ensure at least 1 population is alive
################

def isWorldStillAlive(lattice, gen):
    i = 0
    for j in range(0,len(config_3.populations)):
        x = config_3.populations[j][0]
        y = config_3.populations[j][1]
        if lattice[x,y].exists == True:
            return True
    return False

########################
## Extinct Population ##
########################

#This method is ran in 2 scenarios: when a catastrophe occurs in the 
#populations or when there becomes fewer than 4 haploid chromosomes

def extinctPop(lattice, i, j, gen):
    lattice[i,j].exists = False #since the population no longer exists
    lattice[i,j].ne = 0 #There are no individuals in it now
    lattice[i,j].snps.freq = [] #And no snps now
    lattice[i,j].snps.sel = []
    lattice[i,j].snps.pos = []
    lattice[i,j].snps.arose = []
    return lattice

###############
## Migration ##
###############

def mig(lattice, gen):
    if isinstance(config_3.migProb, (int,float,complex)):
        origins = rn.sample(config_3.populations, len(config_3.populations)) #This is just shuffling them
        for f in range(0, len(config_3.populations)): #For each pop
            x = origins[f][0]
            y = origins[f][1]
            if lattice[x,y].exists == True and lattice[x,y].ne > 0:
                prob = config_3.migProb
                if config_3.diploid == True:
                    a = np.random.binomial( int(lattice[x,y].ne/2.0), prob) #Num migrates away
                else:
                    a = np.random.binomial(lattice[x,y].ne, prob)
                
                #Now I must choose where the migrated to.
                for i in range(0,a): #Need to determine where each migrant goes, separately
                    success = False
                    while success == False:
                        dest = np.random.random_integers(0, (len(origins) - 1)) #Pick dest pop (slow if the pop grid is sparse)
                        g = origins[dest][0] #g + h are the coordinates of the destination grid point
                        h = origins[dest][1]
                        if g != x or h != y:
                            if abs(g-x) <= 1 and abs(h-y) <= 1 and config_3.popMat[g][h] != -1: #If it isn't migrating back to itself, but is migrating to a neighbor
                                success = True
                                
                                if lattice[g,h].exists == True and lattice[g,h].ne > 0: #If dest has peeps, we change the allele freqs
                                #Above changed to add "exists" part 8 Oct 2014
                                    if config_3.diploid == False:
                                        lattice = changeAlleleFreqs(lattice, x, y, g, h, 1, gen)
                                    else:
                                        lattice = changeAlleleFreqs(lattice, x, y, g, h, 2, gen)
                                elif lattice[g,h].ne == 0 or lattice[g,h].exists == False or lattice[g,h].ne == None: #If it doesn't, we must found a population with this individual
                                    if config_3.diploid == False:
                                        lattice = foundPop(lattice, x, y, g, h, 1, gen)
                                    else:
                                        lattice = foundPop(lattice, x, y, g, h, 2, gen)
                                else:
                                    print("Issue with migration")
                                    sys.exit()
                c = time.time()
    
    else: #If migProb is instead a matrix
        shuffled_migProb = rn.sample(config_3.migProb, len(config_3.migProb)) #shuffles the order

        for i in range(0,len(shuffled_migProb)): #I believe this should give the number of rows
            sx, sy, fx, fy = list(map(int, shuffled_migProb[i][0:4]))
            if lattice[sx,sy].exists == True and lattice[sx,sy].ne > 0 and config_3.popMat[fx][fy] != '-1':
                prob = shuffled_migProb[i][4]
                if config_3.diploid == True:
                    a = np.random.binomial( int(lattice[sx,sy].ne/2.0), prob) #Num migrated away
                else:
                    a = np.random.binomial( lattice[sx,sy].ne, prob )
                
                if a > 0:
                    if lattice[fx,fy].ne > 0:
                        if config_3.diploid == False:
                            lattice = changeAlleleFreqs(lattice, sx, sy, fx, fy, 1, gen)
                        else:
                            lattice = changeAlleleFreqs(lattice, sx, sy, fx, fy, 2, gen)
                    elif lattice[fx,fy].ne == 0:
                        if config_3.diploid == False:
                            lattice = foundPop(lattice, sx, sy, fx, fy, 1, gen)
                        else:
                            lattice = foundPop(lattice, sx, sy, fx, fy, 2, gen)
    
    return lattice

############################################
## Changing allele freqs due to migration ##
############################################

def changeAlleleFreqs(lattice, j, k, g, h, nmig, gen): # migrating from [j,k] to [g,h]
    gh_ne_old = lattice[g,h].ne #This is the Ne before the migrants come in
    lattice[g,h].ne = lattice[g,h].ne + nmig #This is the Ne with the new migrant.  In the way the code is currently set up, nmig always == 1.
    jk_ne_old = lattice[j,k].ne #This is the Ne before the migrants leave
    lattice[j,k].ne = lattice[j,k].ne - nmig #This is the Ne after the migrants leave
    a = np.zeros((1,config_3.nSNPs), dtype='int') #A vector of 0s as long as the current nSNPs, I think?
    b = a[0] #The first 0
    
    for i in range(0, nmig): #Nmig always == 1, the way the code is currently set up.
        s = np.random.random_integers(1,(jk_ne_old-i), config_3.nSNPs) #This returns nSNPs int results, all in the range of 1-jk_ne_old.
        C = (s <= (lattice[j,k].snps.freq - b)) #If the s < (#1s - #1stGone), then it is 1.  If it is greater, then it is a 0.  Why am I subtracting b?..Isn't b just a zero??
        C = C.astype('int') #Turns the "Trues" to 1 and the "Falses" to 0.  I think.
        D = (0 == (lattice[j,k].snps.freq - b)) ##This returns a vector that is "True" for all the SNPs that are at freq 0 in the source pop.
        D = D.astype('int') #This turns all the "Trues" into 1s. And all the "Falses" to 0. I think.
        C = C - D #Subtract 1 from all the SNPs that are at 0 frequency.
        C = C.clip(min=0) #If the above command made anything negative, then it is reset to 0.
        #This is a check to get rid of the boundary condition where C will be '1' if lattice[j,k].snps.freq == 0, but noo
        b = b + C #This is the vector of all of the SNPs in the new migrant
        del C, D
    lattice[g,h].snps.freq = lattice[g,h].snps.freq + b #This is added to the receiving pop
    lattice[j,k].snps.freq = lattice[j,k].snps.freq - b #And subtracted from the source pop
    if lattice[j,k].ne == 0: #If there are now no individuals in the populations, it is extinct.
        extinctPop(lattice,j,k,gen)
    del b, a
    return lattice

#################
## World ended ##
#################

#Method is called if there are no longer any individuals in the whole world.

def worldEnded():
    print('End of the world.  Everyone died.')
    sys.exit()


##############################
## Founding new populations ##
##############################

def foundPop(lattice, j, k, g, h, nmig, gen): #n diploids migrating from [j,k] to [g,h]..currently nmig always == 1
    lattice[g,h].exists = True
    lattice[g,h].ne = nmig
    lattice[j,k].ne = lattice[j,k].ne - nmig
    lattice[g,h].snps.pos = list(lattice[j,k].snps.pos) #makes it a COPY

    if config_3.sDiff == None:
        lattice[g,h].snps.sel = list(lattice[j,k].snps.sel) #makes it a COPY also
    else:
        lattice[g,h].snps.sel = getSpatialSel(lattice, g, h)

    lattice[g,h].snps.arose = list(lattice[j,k].snps.arose)
    a = np.zeros((1,config_3.nSNPs), dtype="int") #Vector of 0s for each nSNPs
    lattice[g,h].snps.freq = a[0] #Vector of 0s for each nSNPs
    b = a[0]
  
    for i in range(0, nmig):
        s = np.random.random_integers(0,(lattice[j,k].ne-i),config_3.nSNPs)
        C = (s <= (lattice[j,k].snps.freq - b))#if the rn < (# 1s-#1stgone), then its a 1. if rn # > 1s, then its a 0.
        C = C.astype('int')
        D = (0 == (lattice[j,k].snps.freq - b)) #1 for all the SNPs at freq 0 in the population.  0 for all the SNPs that don't.
        D = D.astype('int')
        C = C - D
        C = C.clip(min=0)
        #This is a check to get rid of the boundary condition where C will be '1' if lattice[j,k].snps.freq == 0, but noo
        b = b + C
        del C,D
    lattice[g,h].snps.freq = lattice[g,h].snps.freq + b #Added to the receiving pop
    lattice[j,k].snps.freq = lattice[j,k].snps.freq - b #Subtracted from the source pop
    if lattice[j,k].ne == 0 :
        lattice = extinctPop(lattice,j,k)
    del b
    return lattice

##############
## Mutation ##
##############

def getSpatialSel(lattice, g, h):
    popStr = str(g) + str(h)
    a = np.where(config_3.sDiff[0] == popStr)

    snpPositions = np.array(lattice[g,h].snps.pos)
    lattice[g,h].snps.sel = np.tile([0],len(lattice[g,h].snps.pos))
    
    snpStart = 0
    for i in range(1,len(config_3.sDiff)):
        sSNP, fSNP = config_3.sDiff[i][0:2] #This is the beginning and end SNP
        b = np.where(np.logical_and(snpPositions >= sSNP, snpPositions <= fSNP))
        if len(b[0]) > 0:
            snpEnd = snpStart + len(b)
            lattice[g,h].snps.sel[snpStart:snpEnd] = float(config_3.sDiff[i][a])
            snpStart = snpEnd
    return lattice[g,h].snps.sel

def getSpatialSelTwo(num, snpNum, x, y):
    popStr = str(x) + str(y)
    a = np.where(config_3.sDiff[0] == popStr)
    newSel = np.tile([0],num)
    snpPositions = np.arange(snpNum, snpNum + num) #Make sure no off by 1 error!
    snpStart = 0
    for i in range(1,len(config_3.sDiff)):
        sSNP, fSNP = config_3.sDiff[i][0:2] #This is the beginning and end SNP
        b = np.where(np.logical_and(snpPositions >= sSNP, snpPositions <= fSNP))
        if len(b[0]) > 0:
            snpEnd = snpStart + len(b[0])
            newSel[snpStart:snpEnd] = float(config_3.sDiff[i][a])
            snpStart = snpEnd
    snpPositions = []
    return newSel



def mut(lattice, x, y, gen): #Each allele has a 1.0*10^-9 probability of mutation each year
    #Assuming infinitely many sites model
    l = config_3.gSize * lattice[x,y].ne * config_3.mutRate #This is lamdba, the mean
    mutations = np.random.poisson(lam=1.0) #This returns the number of mutations in this pop in this gen
    newSNPNums = np.arange(config_3.snpNum+1, config_3.snpNum+mutations+1)
    config_3.snpNum = config_3.snpNum + mutations
    config_3.nSNPs = config_3.nSNPs + mutations #current # of SNPs after all these mutations in this population
    a = np.zeros((1,mutations),dtype='int')
    b = a[0]
    if config_3.sDiff == None:
        sel_b = sel(mutations)
    else:
        sel_b = getSpatialSelTwo(mutations, config_3.snpNum, x, y)

    dddd = np.ones(mutations,dtype='int')
    lattice[x,y].snps.arose = np.hstack((lattice[x,y].snps.arose, np.repeat([gen],mutations)))
    lattice[x,y].snps.freq = np.hstack((lattice[x,y].snps.freq, dddd)) #This adds a new SNP to the frequency
    lattice[x,y].snps.pos = np.hstack((lattice[x,y].snps.pos,newSNPNums))
    lattice[x,y].snps.sel = np.hstack((lattice[x,y].snps.sel,sel_b))

    for i in range(0, len(config_3.populations)):
        u = config_3.populations[i][0]
        o = config_3.populations[i][1]
        if (int(u) != int(x)) or (int(o) != int(y)):
            if lattice[u,o].exists == True:
                lattice[u,o].snps.freq = np.hstack((lattice[u,o].snps.freq,b))#Adds 0 to all pop besides the one it arose in
                lattice[u,o].snps.pos = np.hstack((lattice[u,o].snps.pos,newSNPNums))#This adds snpNums+1 - snpNums+mutations
                if config_3.sDiff != None:
                    sel_b = getSpatialSelTwo(mutations, config_3.snpNum, u, o)
                lattice[u,o].snps.sel = np.hstack((lattice[u,o].snps.sel,sel_b))#This is a zero for each mutation
                lattice[u,o].snps.arose = np.hstack((lattice[u,o].snps.arose, np.repeat([gen],mutations)))
    n = time.time()
    
    return lattice


#####################################
#Selection
#####################################
def sel(num):
    #s = np.random.beta(0.5,0.5) - 0.5 #The expected value of a beta(0.5,0.5) appears to be about 0.5
	#For now, a selected allele will occur at different points.  I will set the s and position manually.
	#For the selection phase.  In the drift phase...if the
    if config_3.s != None: #If there is no selection, config_3.s == 0.
        s = [config_3.s] * num
    elif config_3.s_mat != None:
        #I have to find the snpNums..which is the running total.  0 based, while input is 1 based.
        found = False
        s = None
        for line in s_mat:
            if (snpNum >= line[0] and snpNum <= line[1]) or ((snpNum + num) >= line[0] and (snpNum + num) <= line[1]):
                a = max(snpNum, int(s_mat[0]))
                b = min(int(s_mat[1]), snpNum + num)
                added_s = [ float(line[2]) ] * (b - a + 1)
                s = np.hstack( (s,added_s) )
        if len(s) < num:
            s = np.hstack( (s, [0] * (num-len(s)) ) ) #if s runs out, but there are still more
    
    else:
        s = np.random.gamma( shape = config_3.g_a, scale = config_3.g_b, size=num)
        aa = np.where(s > 1) #restrict s to 0-1
        s[aa] = 1
        bb = np.where(s < 0)
        s[bb] = 0
    return s

#Selection changing in one population
def selChange(pop, selected_Positions, oldPopNe):
    pop.snps.sel = np.array(pop.snps.sel)
    pop.snps.freq = np.array(pop.snps.freq)
    s = pop.snps.sel[ selected_Positions ]
    p = pop.snps.freq [ selected_Positions ] / float(oldPopNe)
    q = 1 - p
    h_array = [ config_3.h ] * len(p)
    h_array = np.array(h_array)
    w = 1 - (2*p*q*h_array*s) - ((q**2) * s)
    w_zeros = np.where(w == 0)
    w_good = np.where(w > 0)
    w_good_1 = w_good[0]
    w_good_2 = w_good_1.astype(int)
    A = np.zeros(len(w))
    A[w_good_2] = (p[w_good_2]**2 + (p[w_good_2]*q[w_good_2]*(1-h_array[w_good_2]*s[w_good_2]))) / w[w_good_2] #A is the new introduced selected allele
    nextGen = np.random.binomial(pop.ne, A)
    #iii =  np.where(nextGen == min(nextGen))
    p = nextGen
    #hh = np.where(p < 0)
    #p[hh] = 0
    h = np.where(p > pop.ne) #Why is this here? What am I doing with this? 7 Oct 2014
    if len(h[0] > 0):
        print("ERROR: PROBLEM WITH SELECTION")
        sys.exit()
    return p

#############
## Growth  ##
#############
def growth(j,k,ne):
    
    if isinstance(config_3.rMean, (int, float, complex)):
        rMean_toUse = config_3.rMean
    else:
        rMean_toUse = config_3.rMean[j][k]
    if isinstance(config_3.rVar, (int, float, complex)):
        rVar_toUse = config_3.rVar
    else:
        rVar_toUse = config_3.rVar[j][k]
    if isinstance(config_3.A, (int, float, complex)):
        A_toUse = config_3.A
    else:
        A_toUse = config_3.A[j][k] 
    if isinstance(config_3.K, (int, float, complex)):
        K_toUse = config_3.K
    else:
        K_toUse = config_3.K[j][k]
    
    r = rn.normalvariate(rMean_toUse, math.sqrt(rVar_toUse))
    newNe = r * ne * (1 - ne / float(K_toUse))* ((ne - A_toUse) / float(K_toUse))
    
    newNe = ne + newNe
    newNe = int(round(newNe))

    if newNe < 0:
        newNe = 0
    
    return int(newNe)

###########
## Drift ##
###########
def drift(lattice, j, k, gen):
    # newPs1 = np.zeros((1,config_3.nSNPs), dtype=int)
    # newPs = newPs1[0] #array of 0s
    newPs = np.zeros((1,config_3.nSNPs), dtype=int)[0]
    oldPopNe = lattice[j,k].ne #Ne before drift and growth
    
    if config_3.rMean != 0 or config_3.rVar !=0 or config_3.A != 0:
        lattice[j,k].ne = growth(j, k, lattice[j,k].ne)
    
    if config_3.diploid == True and lattice[j,k].ne < 4:
        lattice = extinctPop(lattice, j, k, gen)
    elif config_3.diploid == False and lattice[j,k].ne < 2:
        lattice = extinctPop(lattice, j, k, gen)
    else:
        newPs = np.random.binomial(lattice[j,k].ne, ( np.array(lattice[j,k].snps.freq) / float(oldPopNe) ) )
        selected_Positions = [i for i,v in enumerate(lattice[j,k].snps.sel) if v > 0 or v < 0]
        if len(selected_Positions) > 0:
            newPs[ selected_Positions ] = selChange(lattice[j,k], selected_Positions, oldPopNe) #This just updates those selected positions, without changing the other new ones, I think :P
        lattice[j,k].snps.freq = newPs
    newPs = []
    selected_Positions = []
    return lattice

#############################################################
##Lose all SNPs that are at frequency 0 in all populations ##
#############################################################

def loseZeroes(lattice):
    if len(lattice[0,0].snps.freq) > 1:
        A = [] #freq
        pos1 = [] #pos
        sel1 = []
        arose1 = []
        for i in range(0, len(config_3.populations)):
            x = config_3.populations[i][0]
            y = config_3.populations[i][1]
            if lattice[x,y].exists == True:
                A.append(lattice[x,y].snps.freq)
                pos1.append(lattice[x,y].snps.pos)
                sel1.append(lattice[x,y].snps.sel)
                arose1.append(lattice[x,y].snps.arose)
        B = np.sum(A, axis=0) #This returns a vector with the col sums
        z1 = np.zeros((1,( config_3.nSNPs )), dtype='int')
        z = z1[0]
        C = B > z #This is the 0/1 for the results...C is TRUE if the col sum is greater than 0, and FALSE if it is 0.  (It should never be below 0)
        AA = np.asarray(A)
        D = AA[:,C] #This is now the A matrix without the zeros
        pos2 = np.asarray(pos1)
        sel2 = np.asarray(sel1)
        arose2 = np.asarray(arose1)
        pos11 = pos2[:,C]
        sel11 = sel2[:,C]
        arose11 = arose2[:,C]
        k = 0
        config_3.nSNPs = len(D[0,]) #nSNPs is the running counter of the number of snps currently active
        for i in range(0, len(config_3.populations)):
            x = config_3.populations[i][0]
            y = config_3.populations[i][1]
            if lattice[x,y].exists == True:
                lattice[x,y].snps.freq = list(D[k,])
                lattice[x,y].snps.pos = list(pos11[k,])
                lattice[x,y].snps.sel = list(sel11[k,])
                lattice[x,y].snps.arose = list(arose11[k,])
                k = k + 1
    return lattice

#A = np.array([ [0.1,0.1,0,0.5],[0.1,0.2,0,0],[0,0.9,0,1],[0.2,0.1,0,0.3] ])
# B = np.sum(A, axis=0)
# z = [0] * 4
# C = B > z
# D = A[:,C] #IT WORKS


##################################
## Updating Allele History
## columns = (1.Gen) (2.siteNum) (3.genArose) (4.fitness) (5.# derived copies in pop(0,0)) (6. pop(0,0) ne) ....
##Somehow I also need to output which is which
##################################
def updateAlleleHistory(lattice, gen):
    x = config_3.populations[0][0] #if there are no snps here, this won't work for keeping track of numSNPs
    y = config_3.populations[0][1]
    
    a = len(config_3.populations)*2 + 4 #this is the # of columns (rows for now) I will need
    
    temp_list = np.zeros( (a,len(lattice[x,y].snps.freq)) ) #a rows, len(lattice.snps.freq) columns
    temp_list[0] = np.tile( gen, len(lattice[x,y].snps.freq) ) #generation
    temp_list[1] = list(lattice[x,y].snps.pos) #SNP nums for everything.  Hopefully a deep copy
    temp_list[2] = list(lattice[x,y].snps.arose)
    #temp_list[2] is the genArose, which is 0 for all of these, change this to do something...make another structure in "site"
    temp_list[3] = list(lattice[x,y].snps.sel) #SNP sel for everything.  Hopefully a deep copy.
    k = 4
    for i in range(0, len(config_3.populations)):
        x = config_3.populations[i][0]
        y = config_3.populations[i][1]
        if len(lattice[x,y].snps.freq) > 0:
            temp_list[k] = list(lattice[x,y].snps.freq)
            temp_list[k+1] = np.tile( lattice[x,y].ne, len(lattice[x,y].snps.freq) )
        k = k + 2
    global allele_history
    config_3.allele_history.append(temp_list.transpose())


############################
##Checking for edge cases ##
############################
##Make sure there are SNPs
##Make sure there are people
##If no people, end the world
##If no SNPs, give an error
##Fix the output for the allele trajectories


########################################
## Writing final allele freqs to file ##
########################################

def writeToFile(lattice, gen):
    
    f = open(str(config_3.outFile) + ".results.gen" + str(gen), 'wb') #added the b, 8 Oct 2014
    lineToWrite = str(config_3.max_X) + ", " + str(config_3.max_Y) + ", " + str(config_3.nSNPs) + "\n"
    f.write(bytes(lineToWrite, 'UTF-8'))
    #f.write(str(config_3.max_X) + ", " + str(config_3.max_Y) + ", " + str(config_3.nSNPs) + "\n")
    for i in range(0, len(config_3.populations)):
        x = config_3.populations[i][0]
        y = config_3.populations[i][1]
        lineToWrite2 = str(lattice[x,y].ne) + ", " + str(x) + ", " + str(y) + "\n"
        f.write(bytes(lineToWrite2, 'UTF-8'))
        f.write(bytes(', '.join(map(str, lattice[x,y].snps.pos)) + "\n", 'UTF-8'))
        f.write(bytes(', '.join(map(str, lattice[x,y].snps.sel)) + "\n", 'UTF-8'))
        f.write(bytes(', '.join(map(str, lattice[x,y].snps.freq)) + "\n", 'UTF-8'))
        f.write(bytes(', '.join(map(str, lattice[x,y].snps.arose)) + "\n", 'UTF-8'))
    f.close()
    
    if config_3.recordTrag != 0:
        aaa = config_3.allele_history[0]
        for j in range(0,len(config_3.allele_history)): #FIX THIS. Allele_history isn't a vector!
            aaa1 = np.vstack((aaa,config_3.allele_history[j]))
            aaa = np.around(aaa1, decimals=3) #6 October 2014 added
        alleleFileName = str(config_3.outFile) + ".alleleHistory.gen" + str(gen)
        headerToWrite = 'gen siteNum genAlleleArose s '
        for i in range(0, len(config_3.populations)):
            headerToWrite = headerToWrite + 'pop' + str(config_3.populations[i][0]) + str(config_3.populations[i][1]) + 'Freq ' + 'pop' + str(config_3.populations[i][0]) + str(config_3.populations[i][1]) + 'Ne '
        headerToWrite = headerToWrite + "\n"
        with open(alleleFileName, 'wb') as f:
            f.write(bytes(headerToWrite, 'UTF-8'))
            #np.savetxt(f, aaa, fmt='%i', delimiter='\t')
            np.savetxt(f, aaa, fmt='%i', delimiter='\t')
        f.close()
        
        ##Writing to different types of files
    if len(lattice[config_3.populations[0][1],config_3.populations[0][1]].snps.freq) == 0:
        print("NO MORE SNPS")
        sys.exit()
    else:
        for i in range(0, len(config_3.populations)):
            x = config_3.populations[i][0]
            y = config_3.populations[i][1]
            if lattice[x,y].exists == True:
                indivMat = writeOutput_3.makeIndiv(lattice, x, y, config_3.diploid)
                if config_3.GENEPOP:
                    writeOutput_3.GENEPOP(lattice, x, y, indivMat, gen)
                    writeOutput_3.makeCoorFile(lattice, config_3.populations , x, y, gen)
                if config_3.GENELAND:
                    writeOutput_3.GENELAND(lattice, x, y, indivMat, gen)
                if config_3.PLINK:
                    writeOutput_3.PLINK(lattice, x, y, indivMat, gen)


#This is the standard boilerplate that calls the main() function.
if __name__ == '__main__':
    main()