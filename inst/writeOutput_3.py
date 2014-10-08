########################
## 28 October 2013
## This code contains functions to create various output files from the data
########################
#import global_file
import config_3
import numpy as np
import os
import math
import setup_program_3

#######################
## Make individuals
## The goal of this is to convert the population level data into individual level data
## This is just done by sampling w/o replacement for the genotypes at each locus.
## Since they are all unlinked loci, I think this is a valid thing to do.
##
## If haploid, it returns an array where each column is a site, and each row in an individual
## If diploid, it returns an array where there are 2 columns for each site (1 for each of the two alleles), and
######################

def makeIndiv(lattice, i, j, diploid):
    a = lattice[i,j].snps.freq #currently this is an array with the # of individuals with that allele.
    indivSNPs=0
    if diploid == False:
        for k in range(0,len(a)):
            # b = np.hstack( ( np.ones(math.ceil(a[k]*lattice[i,j].ne)), np.zeros(lattice[i,j].ne - math.ceil(lattice[i,j].ne*a[k]) ) ) )
            b = np.hstack( ( np.ones(a[k]), np.zeros(lattice[i,j].ne - a[k]) ) )
            #  print b
            np.random.shuffle(b)#This just randomly shuffles b in place
            if k == 0:
                indivSNPs = b
            else:
                indivSNPs = np.column_stack( [ indivSNPs , b ] ) #This is creating the array of everything
    else: #Make diploids.
        if type(a) == int: #Then there will just be 1 iteration
            iterations = a
        else:
            iterations = len(a)
        for k in range(0, iterations):
            dipVector = []
            if type(a) == int:
                currSNP = a
            else:
                currSNP = a[k]
            b = np.hstack( (np.ones(currSNP), np.zeros(lattice[i,j].ne - currSNP)) )
            np.random.shuffle(b)
            for l in range(0, int(len(b)/2.0)):
                if l == 0:
                    dipVector = (b[l*2], b[l*2+1])
                else:
                    dipVector = np.vstack( (dipVector, (b[l*2], b[l*2+1])) )
            if k == 0:
                indivSNPs = dipVector
            else:
                indivSNPs = np.column_stack( [ indivSNPs , dipVector ] ) #This is creating the array of everything
    return indivSNPs #This is the individual data, but haploid.

##GENEPOP
##Ex
#Microsat on Chiracus radioactivus, a pest species
#Loc1, Loc2, Loc3, Y-linked, Loc4
#POP
#AA8, 0405 0711 0304 0000      0505
#AA9, 0405 0609 0208 0000      0505
#A10, 0205 0609 0101 0000      0305
#Pop
#AF, 0000 0000 0000 0000      0505
#AF, 0205 0307 0102 0000      0505
#AF, 0202 0609 0202 0000      0505
#pop
#C45, 0505 0606 0202 0000      0505
#C45, 0505 0909 0202 0000      0505
#C45, 0505 0306 0202 0000      0505
#This can be done one at a time and appended to

def GENEPOP(lattice, i, j, indivMat, gen):
    a = lattice[i,j].snps.freq #currently this is an array with the # of individuals with that allele.
    inputFile = str(config_3.outFile) + ".GENEPOP.gen" + str(gen)
    f = open(inputFile, 'a')
    if is_non_zero_file(inputFile) == False:
        f.write("GENEPOP output file \n")
        if type(a) == int:
            f.write("Loc1 \n")
        else:
            locNames_1 = ["Loc"] * len(a)
            locNames_2 = np.arange(1,len(a)+1)
            locNames = [ ( locNames_1[0] + str(locNames_2[0]) ) ]
            for k in range(0,len(a)):
                locNames.append( str( locNames_1[k] + str(locNames_2[k]) ) )
            #arr = ', '.join(map(str, locNames))
            arr = "\n".join(map(str, locNames)) #Kimberly added 22 September 2014
            f.write(arr)
            f.write("\n")
    f.write("POP\n")
    if config_3.diploid == False:
        numRows = lattice[i,j].ne
        for k in range(0,numRows):
            f.write("POP" + str(i) + str(j) +"_" + str(k) + ", ")
            if isinstance(indivMat[0],(int,float,complex)):
                if indivMat[k][0] == 0:
                    f.write(" 01")
                else:
                    f.write(" 02")
            else:
                for l in range(0, len(indivMat[0])): #for each snp in the kth individual
                    if indivMat[k][l] == 0:
                        f.write(" 01")
                    else:
                        f.write(" 02")
            f.write("\n")
    else: #If its diploid
        numRows = int(lattice[i,j].ne / 2) #If they are diploid, there will be half the # of rows

        for k in range(0,numRows):
            f.write("POP" + str(i) + str(j) +"_" + str(k) + ", ")
            for l in range(0, int(len(indivMat[0])/2) ): #for every two alleles of the snps in the kth individual
                if indivMat[k][l*2] == 0:
                    a = "01"
                if indivMat[k][l*2] == 1:
                    a = "02"
                if indivMat[k][l*2+1] == 0:
                    b = "01"
                if indivMat[k][l*2+1] == 1:
                    b = "02"
                f.write(a + b)
                f.write(" ")
            f.write("\n")
    f.close()

def GENELAND(lattice, i, j, indivMat, gen):
    a = lattice[i,j].snps.freq
    inputFile = str(config_3.outFile) + ".GENELAND.gen" + str(gen)
    f = open(inputFile, 'a')
    if config_3.diploid == False:
        numRows = lattice[i,j].ne
        for k in range(0,numRows):
            for l in range(0, len(indivMat[0])): #for each snp in the kth individual
                if indivMat[k][l] == 0:
                    f.write("01 ")
                else:
                    f.write("02 ")
            f.write("\n")
    
    else: #If its diploid
        numRows = int(lattice[i,j].ne / 2) #If they are diploid, there will be half the # of rows
        for k in range(0,numRows):
            for l in range(0, ( int(len(indivMat[0])/2)) ): #for every two alleles of the snps in the kth individual
                if indivMat[k][l*2] == 0:
                    a = "01"
                if indivMat[k][l*2] == 1:
                    a = "02"
                if indivMat[k][l*2+1] == 0:
                    b = "01"
                if indivMat[k][l*2+1] == 1:
                    b = "02"
                f.write(a + "/" + b + " ")
            f.write("\n")
    
    f.close()

def makeCoorFile(lattice, nPops, i, j, gen):
    inputFile = str(config_3.outFile) + ".GENEPOP.PopCoor.gen" + str(gen)
    f = open(inputFile, 'a')
    if config_3.diploid == False:
        for k in range(0,lattice[i,j].ne):
            f.write(str(i) + " " + str(j) + "\n")
    else:
        for k in range(0,(lattice[i,j].ne/2)):
            f.write(str(i) + " " + str(j) + "\n")
    f.close()

#FamID IndvID PatID MatID Sex Phenotype Alleles, (genotypes can be 1,2,3,4)
#FAM001  1  0 0  1  2  A A  G G  A C
#FAM001  2  0 0  1  2  A A  A G  0 0

def PLINK(lattice, i, j, indivMat, gen):
    a = lattice[i,j].snps.freq
    
    ##Make PED file
    inputFile = str(config_3.outFile) + ".ADMIXTURE.PED.gen" + str(gen)
    f = open(inputFile, 'a')
    if config_3.diploid == False:
        numRows = lattice[i,j].ne
        for k in range(0,numRows):
            f.write("POP" + str(i) + str(j) + " " + str(k) + " 0 0 -9 -9") #POP00 k 0 0 -9 -9
            for l in range(0, len(indivMat[0])): #for each snp in the kth individual
                if indivMat[k][l] == 0:
                    f.write(" 1")
                else:
                    f.write(" 2")
            f.write("\n")
    else: #If its diploid
        numRows = int(lattice[i,j].ne / 2) #If they are diploid, there will be half the # of rows
        for k in range(0,numRows):
            f.write("POP" + str(i) + str(j) + " " + str(k) + " 0 0 -9 -9") #POP00 k 0 0 -9 -9
            for l in range(0, ( int(len(indivMat[0])/2)) ): #for every two alleles of the snps in the kth individual
            #Kimberly added "int" b/c python3 said it was a float
                if indivMat[k][l*2] == 0:
                    a = "1"
                if indivMat[k][l*2] == 1:
                    a = "2"
                if indivMat[k][l*2+1] == 0:
                    b = "1"
                if indivMat[k][l*2+1] == 1:
                    b = "2"
                f.write(" " + a + " " + b)
            f.write("\n")
    f.close()
    
    ##Make Map file
    #chr snpID    genDist bpPos
    # 1  rs123456  0  1234555
    # 1  rs234567  0  1237793
    # 1  rs233556  0  1337456
    a = lattice[i,j].snps.freq
    if i == 0 and j == 0:
        inputFile = str(config_3.outFile) + ".ADMIXTURE.MAP.gen" + str(gen)
        f = open(inputFile, 'a')
        if type(a) == int:
            f.write("0 Loc0 0 1\n")
        else:
            locNames_1 = ["Loc"] * len(a)
            locNames_2 = np.arange(1,len(a)+1)
            locNames = [ ( locNames_1[0] + str(locNames_2[0]) ) ]
            for k in range(1,len(a)):
                locNames.append( str( locNames_1[k] + str(locNames_2[k]) ) )
            for k in range(0, len(a)):
                f.write("0 " + str(locNames[k]) + " 0 " + str(k) + "\n")
        f.close()

def is_non_zero_file(fpath):
    return True if os.path.isfile(fpath) and os.path.getsize(fpath) > 0 else False

