#########################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: rf_2pops.R                                                    ##
## Contains: rf_2pops                                                  ##
##                                                                     ##
## Written by Rodrigo R Amadeu and Marianella Quezada                  ##
##                                                                     ##
## First version: 20/12/2017                                           ##
## Last update: 23/04/2018                                             ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#########################################################################

## Function to estimate multipoint recombination fraction for composite populations

##' Multipoint recombination fraction based on two mapping populations
##'
##' The function `rf_2pops` estimates the recombination fraction based on two mapping populations.
##' It estimates the recombination fractions based on a multipoint approach implemented using the
##' methodology of Hidden Markov Models (HMM) with the Expectation Maximization (EM) algorithm as
##' explained in the supplementary material of Quezada et al. (2017).
##'
##'
##' @aliases rf_2pops
##' @param marker_names an vector of strings with the marker order.
##' @param data_P1 an onemap object class with the information of the first population.
##' @param data_P2 an onemap object class with the information of the second population.
##' @param rftwopoints_P1 an rf_2pts object with the information of rf from the first population.
##' @param rftwopoints_P2 an rf_2pts object with the information of rf from the second population.
##' @param LOD minimum LOD Score to declare linkage (defaults to \code{3}).
##' @param max.rf maximum recombination fraction to declare linkage (defaults
##' to \code{0.50}).
##' @param log10.mintol log base 10 of the minimum tolerance to converge the EM (defaults to \code{-6})
##' @param max_it maximum number of iteration to converge the EM (defaults to \code{50})
##'
##' @return A list with two objects of class \code{rf_2pts}, one for each population and the
##' log-likelihood of the composite  map
##'
##' @note The thresholds used for \code{LOD} and \code{max.rf} will be used in
##' subsequent analyses, but can be overriden.
##' @author Rodrigo R Amadeu \email{rramadeu@@gmail.com} and Marianella Quezada \email{marianellaquezada@@gmail.com}
##'
##' @references Quezada et al. 2018
##'
##' @keywords utilities
##'
##' @examples
##' notrun{
##' data(onemap2pop)
##' #Extracting the marker names:
##' order_LG1POP1 <- colnames(POP1.geno$geno)[LG1_POP1_final$seq.num]
##' #Computing the rf and likelihood considering information of POP1 and POP2
##' LG1_POP1order <- rf_2pops(markers_names=order_LG1POP1,
##'                           data_P1=POP1_geno,
##'                           data_P2=POP2_geno,
##'                           rftwopoints_P1 = twopts_POP1,
##'                           rftwopoints_P2 = twopts_POP2,
##'                           LOD=3,
##'                           max.rf=0.5,
##'                           log10.mintol = -6,
##'                           max_it = 60)
##' }
##'

rf_2pops <- function(markers_names,
                     data_P1,
                     data_P2,
                     rftwopoints_P1,
                     rftwopoints_P2,
                     LOD=3,
                     max.rf=0.5,
                     log10.mintol = -6,
                     max_it = 50){

    ## Checking inputed data
    if(!any(c("onemap", "outcross") %in% class(data_P1)))
        stop(deparse(substitute(data_P1))," is not an object of class 'onemap'.")

    if(!any(c("onemap", "outcross") %in% class(data_P2)))
        stop(deparse(substitute(data_P2))," is not an object of class 'onemap'.")

    if(!is(rftwopoints_P1, "rf_2pts"))
        stop(deparse(substitute(rftwopoints_P1))," is not an object of class 'rf_2pts'")

    if(!is(rftwopoints_P2, "rf_2pts"))
        stop(deparse(substitute(rftwopoints_P2))," is not an object of class 'rf_2pts'")

    ## Getting position inside each geno data
    seqnum_P1 <- match(markers_names,colnames(data_P1$geno))
    seqnum_P2 <- match(markers_names,colnames(data_P2$geno))

    if(!prod(!is.na(seqnum_P1))){
        print(paste(markers_names[which(!is.na(seqnum_P1))],sep=" "))
        stop("Above markers not in population 1. Double-check it.")
    }

    if(!prod(!is.na(seqnum_P2))){
        print(paste(markers_names[which(!is.na(seqnum_P2))],sep=" "))
        stop("Above markers not in population 2. Double-check it.")
    }

    ## Counting total mark information (not NA) between consecutive markers
    ## to be used as weighted on M step, eg.:
    ## first element is the total shared markers information between marker 1 and marker 2
    ## second element between marker 2 and marker 3, and so on
    ## For weighted mean
    ## n_P1 <- data_P1$geno[,seqnum_P1]!=0
    ## n_P1 <- diag(crossprod(n_P1)[-1,])

    ## n_P2 <- data_P2$geno[,seqnum_P2]!=0
    ## n_P2 <- diag(crossprod(n_P2)[-1,])

    ## For a scalar mean
    n_P1 <- nrow(data_P1$geno)
    n_P2 <- nrow(data_P2$geno)

    ## Transforming in weights
    n_P1P2 <- n_P1+n_P2 ## Total
    n_P1 <- n_P1/n_P1P2 ## Weights for P1
    n_P2 <- n_P2/n_P1P2 ## Weights for P2

    ## make_seq and map on both populations
    makeseq_P1 <- make_seq(rftwopoints_P1,seqnum_P1)
    makeseq_P1$twopt <- deparse(substitute(rftwopoints_P1))
    mapseq_P1 <- map(makeseq_P1)
    makeseq_P2 <- make_seq(rftwopoints_P2,seqnum_P2)
    makeseq_P2$twopt <- deparse(substitute(rftwopoints_P2))
    mapseq_P2<- map(makeseq_P2)

    ## getting initial rf for separated pop
    rfinit_P1 <- onemap:::get_vec_rf_out(mapseq_P1, acum=FALSE)
    rfinit_P2 <- onemap:::get_vec_rf_out(mapseq_P2, acum=FALSE)

    ## EM: Make the first E step outside de loop, than M and E progressivaly loop by loop
    ## tol limit control the loop by loop step

    ## First E step - within populations ##
    Emap_P1 <- onemap:::est_map_hmm_out(geno=t(get(mapseq_P1$data.name, pos=1)$geno[,mapseq_P1$seq.num]),
                               type=get(mapseq_P1$data.name, pos=1)$segr.type.num[mapseq_P1$seq.num],
                               phase=mapseq_P1$seq.phases,
                               rf.vec=rfinit_P1,
                               verbose=FALSE,
                               tol=0.1)

    Emap_P2<- onemap:::est_map_hmm_out(geno=t(get(mapseq_P2$data.name, pos=1)$geno[,mapseq_P2$seq.num]),
                               type=get(mapseq_P2$data.name, pos=1)$segr.type.num[mapseq_P2$seq.num],
                               phase=mapseq_P2$seq.phases,
                               rf.vec=rfinit_P2,
                               verbose=FALSE,
                               tol=0.1)


    log10.mintol <- log10.mintol+1
    seq.tol <- c(rep(seq(1,.2,-.1),abs(log10.mintol))*
                 rep(10^(seq(-1,log10.mintol,-1)),each=9),10^(log10.mintol-1))


    tol <- seq.tol[1]
    j <- 1
    loglike.old <- Emap_P1$loglike + Emap_P2$loglike

## EM Alghoritm
    for(k in 1:max_it){
        ## M STEP - between populations ##

        ## Simple Mean
        ## M_rf <- rowMeans(cbind(Emap_P2$rf, Emap_P1$rf), na.rm=TRUE)

        ## Weighted Mean
        M_rf <- Emap_P1$rf*n_P1 + Emap_P2$rf*n_P2

        ## E STEP ##
        Emap_P1 <- onemap:::est_map_hmm_out(geno=t(get(mapseq_P1$data.name, pos=1)$geno[,seqnum_P1]),
                                           type=get(mapseq_P1$data.name, pos=1)$segr.type.num[seqnum_P1],
                                           phase=mapseq_P1$seq.phases,
                                           rf.vec=M_rf,
                                           verbose=FALSE,
                                           tol=tol)

        Emap_P2 <- onemap:::est_map_hmm_out(geno=t(get(mapseq_P2$data.name, pos=1)$geno[,seqnum_P2]),
                                           type=get(mapseq_P2$data.name, pos=1)$segr.type.num[seqnum_P2],
                                           phase=mapseq_P2$seq.phases,
                                           rf.vec=M_rf,
                                           verbose=FALSE,
                                           tol=tol)


        loglike <- Emap_P1$loglike + Emap_P2$loglike
        if(round(loglike,1) <= round(loglike.old,1)){
            j <- j+1
            tol = seq.tol[j]
        }
        if(j == length(seq.tol)){
            message("reached min tolerance value")
            message(paste("interaction:",k,"; loglike:", loglike))
            break
        }
        loglike.old <- loglike
        if( k%%10 == 0 )
            message(paste("interaction:",k,"; loglike:", loglike, "; tol:", tol))
    }

    M_rf <- Emap_P1$rf*n_P1 + Emap_P2$rf*n_P2

    ## Map with joining rf
    list(
        P1=structure(list(seq.num=seqnum_P1,
                          seq.phases=mapseq_P1$seq.phases,
                          seq.rf=Emap_P1$rf,
                          seq.like=Emap_P1$loglike,
                          data.name=mapseq_P1$data.name,
                          twopt=mapseq_P1$twopt),
                     class = "sequence"),

        P2=structure(list(seq.num=seqnum_P2,
                          seq.phases=mapseq_P2$seq.phases,
                          seq.rf=Emap_P2$rf,
                          seq.like=Emap_P2$loglike,
                          data.name=mapseq_P2$data.name,
                          twopt=mapseq_P2$twopt),
                     class = "sequence"),

        P1P2_seq.like = loglike)

      #  P1P2=structure(list(seq.num=seqnum_P1,
      #                      seq.phases=mapseq_P1$seq.phases,
      #                      seq.rf=M_rf,
      #                      seq.like=loglike,
      #                      data.name=mapseq_P1$data.name,
      #                      twopt=mapseq_P1$twopt),
      #                 class = "sequence")
#    )
}

