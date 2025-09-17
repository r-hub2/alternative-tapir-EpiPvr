## ----setup, include=FALSE-----------------------------------------------------
library(EpiPvr)
knitr::opts_chunk$set(echo = TRUE)

## ----sec1---------------------------------------------------------------------
# EXAMPLE OF PT VIRUS # EXAMPLE OF PT VIRUS # EXAMPLE OF PT VIRUS # EXAMPLE OF PT VIRUS 
virusType="PT"
# set assay structure
nReps=30 # number of reps
numWF=20 # number of insects in cohorts

## ----sec2---------------------------------------------------------------------

# virus/insect rates per hr
alrate=0.1 # acquisition rate
berate=1 # inoculation rate
gamrate=0.5 # latency progression rate  (NA if SPT virus)
murate=0.01 # virus clearance rate


## ----sec3---------------------------------------------------------------------
AAP_lens=c(1,1.5,1.75,2,2.25,2.5,3,4); # variable duration vectors, hours AAP sub-assay
LAP_lens=c(0.5,1,2,3,4,5,6,7,8);                                        # LAP sub-assay
IAP_lens=c(5,10,15,20,25,30,40,50,60)/60;                               # IAP sub-assay

AAP_Reps=rep(nReps,length(AAP_lens));  # vectors for number of replicates
LAP_Reps=rep(nReps,length(LAP_lens));
IAP_Reps=rep(nReps,length(IAP_lens));

AAP_Infs_zeros=rep(0,length(AAP_lens));  # vectors for num infected test plants-initially 0
LAP_Infs_zeros=rep(0,length(LAP_lens));
IAP_Infs_zeros=rep(0,length(IAP_lens));

AAPinput=rbind(AAP_lens, AAP_Reps, AAP_Infs_zeros)  # group structural vectors for input
LAPinput=rbind(LAP_lens, LAP_Reps, LAP_Infs_zeros)
IAPinput=rbind(IAP_lens, IAP_Reps, IAP_Infs_zeros)



## ----sec4---------------------------------------------------------------------

# default durations of acquisition, latent and inoculation periods
T_A=2
T_L=0.5
T_I=1

# Place '-1' in the varied assay component
AAPfixedComponent=c(-1,T_L,T_I)
LAPfixedComponent=c(T_A,-1,T_I)
IAPfixedComponent=c(T_A,T_L,-1)
ddur_mat=rbind(AAPfixedComponent,LAPfixedComponent,IAPfixedComponent)



## ----sec5---------------------------------------------------------------------
# simulate and hence populate the number of infected test plants in the AAP subassay
assay1=AP_assay_simulator(AAPinput,
                          AAPfixedComponent,
                          numWF, c(alrate,berate,gamrate,murate),isVerbose=0,'PT')  
# simulate and hence populate the number of infected test plants in the LAP subassay
assay2=AP_assay_simulator(LAPinput,
                          LAPfixedComponent,
                          numWF, c(alrate,berate,gamrate,murate),isVerbose=0,'PT') 
# simulate and hence populate the number of infected test plants in the IAP subassay
assay3=AP_assay_simulator(IAPinput,
                          IAPfixedComponent,
                          numWF, c(alrate,berate,gamrate,murate),isVerbose=0,'PT') 



## ----sec6---------------------------------------------------------------------
ap_data_sim=list(d_AAP=assay1,      # AAP sub-assay structure 
                 d_LAP=assay2,      # LAP sub-assay structure  (OMIT IF SPT VIRUS) 
                 d_IAP=assay3,      # IAP sub-assay structure 
                 d_durations=ddur_mat,  # fixed durations
                 d_vectorspp=numWF,     # insects in a cohort
                 d_virusType=virusType) # virus code - element of {'SPT','PT'}



## ----sec7---------------------------------------------------------------------
# Adding virus parameters as attributes nb. -1 indicates unknown virus parameters
attr(ap_data_sim, "alpha") <- alrate # acquisition rate
attr(ap_data_sim, "beta") <- berate  # inocualation rate
if (!is.na(gamrate)) {
  attr(ap_data_sim, "gamma") <- gamrate  # latent progression rate(OMIT IF SPT VIRUS) 
}
attr(ap_data_sim, "mu") <- murate  # insect recovery rate rate

print(ap_data_sim)

# example save command
#save(ap_data_sim,file="myvirusdataset.rda"))


## ----sec8, results='asis', echo=FALSE-----------------------------------------

output <- capture.output(print(ap_data_sim))  # Capture console output
writeLines(c("\\begin{verbatim}", output, "\\end{verbatim}"), "console_output.tex")


## ----sec9, echo=FALSE---------------------------------------------------------
#pdf("output_plot1.pdf")
plot(ap_data_sim$d_AAP[1,],ap_data_sim$d_AAP[3,]/ap_data_sim$d_AAP[2,], xlab = "AAP duration, mins", ylab = "Prop. test plants PT infected",ylim=c(0,1))
#dev.off()  

if (!is.null(ap_data_sim$d_LAP)) {
  #pdf("output_plot2.pdf")
  plot(ap_data_sim$d_LAP[1,],ap_data_sim$d_LAP[3,]/ap_data_sim$d_LAP[2,], xlab = "LAP duration, mins", ylab = "Prop. test plants PT infected",ylim=c(0,1))
  #dev.off()  
}

#pdf("output_plot3.pdf")
plot(ap_data_sim$d_IAP[1,],ap_data_sim$d_IAP[3,]/ap_data_sim$d_IAP[2,], xlab = "IAP duration, mins", ylab = "Prop. test plants PT infected",ylim=c(0,1))
#dev.off()  

## ----sec10--------------------------------------------------------------------

# EXAMPLE OF SPT VIRUS # EXAMPLE OF SPT VIRUS # EXAMPLE OF SPT VIRUS  

virusType="SPT"
nReps=30 # number of reps
numWF=20 # number of insects in cohorts

# virus rates per hr
alrate=0.1 # acquisition rate
berate=1 # inoculation rate
gamrate=NA # latency progression rate
murate=1 # virus clearance rate

AAP_lens=c(2,3,3.5,4,4.5,5,6,8); # hours
IAP_lens=c(5,10,15,20,25,30,40,50,60)/60; # hours

AAP_Reps=rep(nReps,length(AAP_lens)); 
IAP_Reps=rep(nReps,length(IAP_lens));

AAP_Infs_zeros=rep(0,length(AAP_lens)); 
IAP_Infs_zeros=rep(0,length(IAP_lens));

AAPinput=rbind(AAP_lens, AAP_Reps, AAP_Infs_zeros)
IAPinput=rbind(IAP_lens, IAP_Reps, IAP_Infs_zeros)

# default durations of acquisition and inoculation periods 
T_A=4
T_I=6
# Place '-1' in active assay component
AAPfixedComponent=c(-1,0,T_I) # 0 hours for LAP fixed duration since SPT virus
IAPfixedComponent=c(T_A,0,-1)



assay1=AP_assay_simulator(AAPinput,
                          AAPfixedComponent,
                          numWF, c(alrate,berate,gamrate,murate),isVerbose=0,virusType)  
assay3=AP_assay_simulator(IAPinput,
                          IAPfixedComponent,
                          numWF, c(alrate,berate,gamrate,murate),isVerbose=0,virusType) 

# dataset format for SPT does not include LAP
ddur_mat=rbind(AAPfixedComponent[c(1,3)],IAPfixedComponent[c(1,3)])

ap_data_sim=list(d_AAP=assay1,
                 d_IAP=assay3,
                 d_durations=ddur_mat,
                 d_vectorspp=numWF,
                 d_virusType=virusType)


## ----sec11--------------------------------------------------------------------
# Adding virus parameters as attributes nb. -1 indicates unknown virus parameters
# note that dataset format does not include gamma atttibute for SPT virus
attr(ap_data_sim, "alpha") <- alrate 
attr(ap_data_sim, "beta") <- berate 
attr(ap_data_sim, "mu") <- murate 


print(ap_data_sim)
#save(ap_data_sim,file=paste0("data/",virusType,"_DELETE_SIM_AccessExp_new.rda"))


## ----sec12, results='asis', echo=FALSE----------------------------------------

output <- capture.output(print(ap_data_sim))  # Capture console output
writeLines(c("\\begin{verbatim}", output, "\\end{verbatim}"), "console_output2.tex")


## ----sec13, echo=FALSE--------------------------------------------------------
#pdf("output_plot1.pdf")
plot(ap_data_sim$d_AAP[1,],ap_data_sim$d_AAP[3,]/ap_data_sim$d_AAP[2,], xlab = "AAP duration, mins", ylab = "Prop. test plants PT infected",ylim=c(0,1))
#dev.off()  

if (!is.null(ap_data_sim$d_LAP)) {
  #pdf("output_plot2.pdf")
  plot(ap_data_sim$d_LAP[1,],ap_data_sim$d_LAP[3,]/ap_data_sim$d_LAP[2,], xlab = "LAP duration, mins", ylab = "Prop. test plants PT infected",ylim=c(0,1))
  #dev.off()  
}

#pdf("output_plot3.pdf")
plot(ap_data_sim$d_IAP[1,],ap_data_sim$d_IAP[3,]/ap_data_sim$d_IAP[2,], xlab = "IAP duration, mins", ylab = "Prop. test plants PT infected",ylim=c(0,1))
#dev.off()  

