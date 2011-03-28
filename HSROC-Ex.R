pkgname <- "HSROC"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('HSROC')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("HSROC")
### * HSROC

flush(stderr()); flush(stdout())

### Name: HSROC
### Title: A function for joint meta-analysis of sensitivity and
###   specificity of a diagnostic test.
### Aliases: HSROC
### Keywords: models

### ** Examples


#===============================================================
#TO SET UP THE REFERENCE STANDARD
#===============================================================


#There were four different reference standards for the In.house dataset.  
#The first reference standard was used in study 1 and 2, 
#the second was used in studies 3 and 4, the third in study 5 and the 
#fourth in studies 6 to 12.
REFSTD = list(4, 1:2, 3:4, 5, 6:12) 

#===============================================================
#TO SET UP DATA AND INITIAL VALUES
#===============================================================

data(In.house)
M = length(In.house[,1])


#Initial values for the within-study parameters
init.alpha = rep(2.5, M) ;	init.theta = rep(1, M) ;
init.s1 =  rep(0.5, M) ;	init.c1 = rep(0.5, M) ;
init.pi = rep(0.5, M)

#Initial values for the between-study parameters
init.THETA = 1 ;	init.sd.theta = 0.5 ;
init.LAMBDA = 2.5 ;	init.sd.alpha = 0.5 ;
init.beta = 0 ;

#Initial values for the reference standard sensitivities and specificities
init.s2 = rep(0.5, REFSTD[[1]]) ;	init.c2 = rep(0.5, REFSTD[[1]])

#The ordering of the initial values is important!
init1 = cbind(init.alpha, init.theta, init.s1, init.c1, init.pi)
init2 = c(init.THETA, init.sd.theta, init.LAMBDA, init.sd.alpha, init.beta)
init3 = rbind(init.s2, init.c2)

init = list(init1, init2, init3)
#===============================================================
#TO PROVIDE PRIOR INFORMATION
#===============================================================

S2.a = c(0.2, 0.2, 0.6, 0.7) ; 	S2.b = c(0.6, 0.7, 0.8, 0.9)
C2.a = rep(0.9, 4) ;	C2.b = rep(1, 4)

#===============================================================
#TO RUN GIBBS SAMPLER
#===============================================================

estimates = HSROC(data=In.house, init=init, iter.num=10,  
   prior.SEref=c(S2.a,S2.b), prior.SPref=c (C2.a,C2.b), sub_rs=REFSTD) 





cleanEx()
nameEx("HSROCSummary")
### * HSROCSummary

flush(stderr()); flush(stdout())

### Name: HSROCSummary
### Title: Summary statistics for HSROC models.
### Aliases: HSROCSummary
### Keywords: models

### ** Examples


#REAL-LIFE EXAMPLES
#
#Example 1
#To get descriptive statistics and graphical summaries for the MRI data 
#(Scheidler et al. 1997) after dropping the first 5,000 iterations.

data(MRI)	#load the data
## Not run: 
##D HSROCSummary(data = MRI, burn_in=5000, print_plot=TRUE )
## End(Not run)


#Example 2
#To get descriptive statistics and graphical summaries for the In.house 
#data (Pai et al. 2004) coming from 2 different chains.  
#We provide the path to each chain's directory, i.e. the directory where 
#all files created during the Gibbs sampler process are stored for 
#each chain.  Let's assume there are two fictional directoies 
#chain_path = list("C:/path_to_chain_1", "C:/path_to_chain_2").
#Let's assume we drop the first 5,000 iterations and we use a thinning 
#interval of 10.

data(In.house)	#load the data
## Not run: 
##D HSROCSummary(data = In.house, burn_in=5000, Thin=10, 
##D 		chain=chain_path, print_plot=TRUE )
## End(Not run)

## Don't show: 

x <- rnorm(1000)
y <- as.mcmc(x)	
z <- HPDinterval(y)

## End Don't show




cleanEx()
nameEx("In.house")
### * In.house

flush(stderr()); flush(stdout())

### Name: In.house
### Title: IN-HOUSE NUCLEIC ACID AMPLIFICATION TESTS (INH) FOR TB PLEURITIS
### Aliases: In.house
### Keywords: datasets

### ** Examples

data(In.house)
In.house




cleanEx()
nameEx("MRI")
### * MRI

flush(stderr()); flush(stdout())

### Name: MRI
### Title: MAGNETIC RESONANCE IMAGING TEST (MRI) for evaluation of lymph
###   node metastasis in women with invasive cervical cancer
### Aliases: MRI
### Keywords: datasets

### ** Examples

data(MRI)
MRI




cleanEx()
nameEx("beta.parameter")
### * beta.parameter

flush(stderr()); flush(stdout())

### Name: beta.parameter
### Title: A function that returns the shape parameters of the beta
###   distribution
### Aliases: beta.parameter
### Keywords: methods

### ** Examples

  

## Not run: beta.parameter(-1, 0.5) #Returns error!
## Not run: beta.parameter(0, 0) #Not allowed.  Returns error!
## Not run: beta.parameter(0.75, 0.25) #Returns error!

beta.parameter(0, 1)
beta.parameter(0.5, 1) 
beta.parameter(0.1, 0.7)            




cleanEx()
nameEx("simdata")
### * simdata

flush(stderr()); flush(stdout())

### Name: simdata
### Title: Simulate a dataset based on a HSROC model
### Aliases: simdata
### Keywords: datagen

### ** Examples


#EXAMPLE 1
#We want to simulate data for 10 studies based on an HSROC model assuming 
#each study uses the same imperfect reference standard.
 
N = 10
LAMBDA = 2
sd_alpha = 0.75
THETA = 1.5
sd_theta = 0.5
beta = 0
pi = runif(10,0,1)

REFSTD = list(1, 1:10)  #Only 1 reference standard ...
s2 = c(0.5)	 #Sensitivity of the reference test
c2 = c(0.85) 	 #Specificity of the reference test	


sim.data = simdata(N=N, n = c(50,50,60,60,70,70,80,80,90,90), 
   sub_rs = REFSTD, prev=pi, se_ref=s2, sp_ref=c2, T=THETA, 
   L=LAMBDA, sd_t=sd_theta, sd_a=sd_alpha, b=beta)



#EXAMPLE 2
#We want to simulate data for 15 studies based on an HSROC model such that 
#the first 5 studies share a common reference standard and the remaining 
#10 studies also share a common reference standard.
 
N = 15
LAMBDA = 3.6
sd_alpha = 1.15
THETA = 2.3
sd_theta = 0.75
beta = 0.15
pi = runif(15,0.1,0.5)

REFSTD = list(2, 1:5, 6:15)  #Two different reference standards ...
s2 = c(0.40, 0.6)	 #Sensitivity of the reference tests
c2 = c(0.75,0.95) 	 #Specificity of the reference tests	

#Thus, for the first 5 studies, S2 = 0.40 and C2 = 0.75 while for the last 
#10 studies s2 = 0.6 and c2 = 0.95


sim.data = simdata(N=N, n=seq(30,120,1), n.random=TRUE, sub_rs = REFSTD, 
   prev=pi,  se_ref=s2, sp_ref=c2, T=THETA, L=LAMBDA, sd_t=sd_theta, 
   sd_a=sd_alpha, b=beta)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
