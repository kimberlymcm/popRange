\name{popRangeSim}
\alias{popRangeSim}
\title{popRangeSim}
\description{
This method calls a python program that simulates stochastic demographic events based on the parameters input by the user.
}
\usage{
popRangeSim(world, popSize, rMean = 0, rVar = 0, A = 0, K = 100, catProb = 0, 
diploid = TRUE, nGens = 100, migration = 0, SNP_model = 1, h = 0.5, s = 0, 
gamma_shape = NULL, gamma_scale = NULL, gSize = NULL, mutRate = NULL, 
nSNPs = NULL, SNPs_starting_freq = NULL, GENEPOP = FALSE, GENELAND = FALSE, 
PLINK = FALSE, outfile = "", infile = NULL, recordTrag = 0, sDiff = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{world}{
This matrix grid defines the grid of populations.  Each point in the matrix defines a potential population.  A 1 indicates a population may exist at that point. -1 indicates no population can exist at that grid point.  This allows the user to define a population grid that is not rectangular.  This can be a data frame, a matrix, or an int. 
}
  \item{popSize}{
Starting population size.  If numeric, every existing population has the same size.  If matrix, each population size is defined separately.  Input must be an integer >= 0 or a matrix of integers >= 0.
}
  \item{rMean}{
Mean r (from exponential growth equation).  Float between c(0,1).  If numeric, all populations have the same r.  If matrix, each population has its own r.
}
  \item{rVar}{
Variance in r.  Float between c(0,1). If numeric, all populations have the same variance.  If matrix, each population has its own r.  Set to 0 if you want r to be exactly rMean each generation.
}
  \item{A}{
Allee Effect. If numeric, all populations have same Allee effect.  If matrix, each population has its own Allee effect.  Set to 0 if you do not want to incorporate the Allee effect
}
  \item{K}{
Carry capacity. If int, all populations have same carry capacity.  If matrix, each population has its own carrying capacity.
}
  \item{catProb}{
%%     ~~Describe \code{catProb} here~~
Probability of population extinction each generation (range:0-1).  If a numeric, all populations have the same catProb. If matrix, each population has its own catProb.

}
  \item{diploid}{
Ploidy of individuals.  Two options: haploid and diploid.  If TRUE, diploid.  If FALSE, haploid.
}
  \item{nGens}{
Int, number of generations for the simulations to run.
}
  \item{migration}{
Probability of each individual migrating each generation.  If float, every population has same migration rate. Each migration randomly picks on adjacent viable grid point to migrate to (diagonals included). Alternatively, a matrix can be provided.  The first 2 columns are the coordinate of the initial population.  The next two columns are the coordinates of the final population.  The fifth column is the migration rate.  Population pairs not included are assumed to have no migration. 
}
  \item{SNP_model}{
There are two options. The first option (0) indicates that the user will provide starting SNP frequencies.  The second option (1) determines starting SNP frequencies according to the standard neutral model with user provided genome size and mutation rate.  
}
  \item{h}{
Dominance parameter.  (Only relevant for diploids.)
}
  \item{s}{
If numeric, all SNPs have the same selection coefficient.  If matrix, first 2 columns define starting and ending SNP.  The 3rd column is the selection coefficient for those SNPs.
}
  \item{gamma_shape}{
Alpha parameters to the gamma distribution.  NOTE: Provide EITHER an s input OR gamma_a, gamma_b values.  
}
  \item{gamma_scale}{
Beta parameters to the gamma distribution.  NOTE: Provide EITHER an s input OR gamma_a, gamma_b values.  
}
  \item{gSize}{
Number of base pairs in the genome. (Required if SNP_Model = 1)
}
  \item{mutRate}{
Mutation rate per generation per site. (Required if SNP_Model = 1)
}
  \item{nSNPs}{
Number of SNPs, if SNP_Model = 0.  (Required if SNP_Model = 0)
}
  \item{SNPs_starting_freq}{
The initial frequency of the SNPs (range:0-1).  If int, all SNPs are given the same starting frequencies.  If matrix, columns 1 and 2 indicate the first and last SNPs.  Column 3 indicates the allele frequency of these SNPs.  (Required if SNP_Model = 0)
}
  \item{GENEPOP}{
If TRUE, write results to GENEPOP formatted file.
}
  \item{GENELAND}{
If TRUE, write results to GENELAND formatted file.
}
  \item{PLINK}{
If TRUE, write results to an PLINK formatted file
}
  \item{outfile}{
A string of the outfile name.
}
  \item{infile}{
If you have already ran a simulation, you can use the results file as input.  This allows for temporally variable parameters and to output the state of the population at various time points.  
}
  \item{recordTrag}{
If n, this records the trajectory of all of the alleles n generations and outputs them to a file.  NOTE: If you have many SNPs or populations, this option may significantly increase simulation time.  Also if you are using SNP_model=1, alleles that have a frequency of 0 in all population will be dropped from simulation and thus will stop being recorded in this file.  This decreases simulation time.
}
  \item{sDiff}{
Defines the selection coefficients for each SNP in each population.  This allows for spatially variable selection coefficients.  The first line of the matrix consists of the strings sSNP and  fSNP.  These are the first and last SNPs with a specific set of selection coefficients.  The coordinates of each population are subsequently listed. 
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
No value is returned.  Output files are written.
}
\references{
%% ~put references to the literature/web site here ~
Please see http://stanford.edu/~kfm/html/research.html 
for extended manual and examples.
}
\author{
%%  ~~who you are~~
\tabular{ll}{
AUTHOR: \tab Kimberly F. McManus \cr
}
}

\examples{ 
## Example: In this scenario, we are simulating a 3x3 grid of populations for 50 
## generations. Each population starts with 100 diploid individuals. Each individual 
## has a 0.01 probability of migrating away from their populations. There are 100 
## SNPs that all have a starting frequency of 0.5, and the program outputs the standard 
## "results" file, as well as a PLINK file. 

## Commands

mat = matrix(1,nrow=3,ncol=3) 
popRangeSim(world = mat, popSize = 100, diploid = TRUE, nGens = 50,
              mig = 0.01, SNP_model = 0, nSNPs = 100, SNPs_starting_freq = 0.5, 
              outfile= "outFile1", PLINK=TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
