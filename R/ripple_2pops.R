#########################################################################
##                                                                     ##
## Package: onemap                                                     ##
##                                                                     ##
## File: ripple_2pops.R                                                ##
## Contains: ripple_2pops                                              ##
##                                                                     ##
## Written by Rodrigo R Amadeu and Marianella Quezada                  ##
##                                                                     ##
## First version: 20/12/2017                                           ##
## Last update: 23/04/2018                                             ##
## License: GNU General Public License version 2 (June, 1991) or later ##
##                                                                     ##
#########################################################################

## This function searches for alternative orders, by comparing all possible
## orders of subsets of markers based on a composite map (rf_2pop)

##' Compares and displays plausible alternative orders for a given linkage
##' group based on a composite map
##'
##' For a given sequence of ordered markers, computes the multipoint likelihood
##' of alternative orders, by shuffling subsets (windows) of markers within the
##' sequence. For each position of the window, all possible \eqn{(ws)!}{(ws)!}
##' orders are compared.
##'
##' Large values for the window size make computations very slow, specially if
##' there are many partially informative markers.
##'
##' @aliases ripple_2pops
##' @param marker_names an vector of strings with the marker order.
##' @param data_P1 an onemap object class with the information of the first population.
##' @param data_P2 an onemap object class with the information of the second population.
##' @param twopts_POP1 an rf_2pts object with the information of rf from the first population.
##' @param twopts_POP2 an rf_2pts object with the information of rf from the second population.
##' @param LOD minimum LOD Score to declare linkage (defaults to \code{3}).
##' @param max.rf maximum recombination fraction to declare linkage (defaults
##' to \code{0.50}).
##' @param log10.mintol log base 10 of the minimum tolerance to converge the EM (defaults to \code{-6})
##' @param max_it maximum number of iteration to converge the EM (defaults to \code{50})
##' @param window an integer specifying the length of the window size (defaults to 4).
##'
##' @return A list with two objects. The first is a matrix with the tested ordered markers,
##' the second is a vector with the respectivally log-likelihood
##'
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
##' #RIPPLE algorithm:
##' LG1_POP1order <- ripple_2pops(markers_names = order_LG1POP1,
##'                               data_P1 = POP1_geno,
##'                               data_P2 = POP2_geno,
##'                               rftwopoints_P1 = twopts_POP1,
##'                               rftwopoints_P2 = twopts_POP2,
##'                               LOD = 3,
##'                               max.rf = 0.5,
##'                               log10.mintol = -6,
##'                               max_it = 50,
##'                               window = 4)
##' }
##'
##'
ripple_2pops <- function(markers_names,
                         data_P1,
                         data_P2,
                         twopts_POP1,
                         twopts_POP2,
                         LOD=3,
                         max.rf=0.5,
                         log10.mintol = -6,
                         max_it = 50,
                         window = 3){
  n <- length(markers_names)
  ripple_list <- vector("list",2)
  for(i in window:n){
    ind <- is.na(match(1:n,(i-window+1):(i)))
    markers_names_perm <- permn(markers_names[(i-window+1):(i)])
    m<-length(markers_names_perm)
    for(j in 1:m){
      print(paste("interaction",i,"/",n,"subset",j,"/",m))
      markers_names_temp <- rep(NA,n)
      markers_names_temp[ind] <- markers_names[ind]
      markers_names_temp[!ind] <- markers_names_perm[[j]]

      like_temp <- rf_2pops(markers_names=markers_names_temp,
                            data_P1=data_P1,
                            data_P2=data_P2,
                            rftwopoints_P1=twopts_POP1,
                            rftwopoints_P2=twopts_POP2,
                            LOD=LOD,
                            max.rf=max.rf,
                            log10.mintol = log10.mintol,
                            max_it = max_it)$P1P2$seq.like
      ripple_list[[2]] <- c(ripple_list[[2]], like_temp)
      ripple_list[[1]] <- rbind(ripple_list[[1]],markers_names_temp)
    }
  }
  return(ripple_list)
}
