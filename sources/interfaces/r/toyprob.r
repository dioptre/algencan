#    =================================================================
#    File: toyprob.r
#    =================================================================
#
#    =================================================================
#    Module: Subroutines that define the problem
#    =================================================================
#
#    Last update of any of the component of this module:
#
#    September 11, 2009.
#
#    Users are encouraged to download periodically updated versions of
#    this code at the TANGO home page:
#
#    www.ime.usp.br/~egbirgin/tango/
#
#    ******************************************************************
#    ******************************************************************

     inip <- function(n,x,l,u,m,lambda,equatn,linear,coded,checkder) {

#    Number of variables

     n <- 2

#    Initial point

     x[1] <- 0.0
     x[2] <- 0.0

#    Lower and upper bounds

     l[1] <- - 10.0
     u[1] <-   10.0

     l[2] <- - 1.0e20
     u[2] <-   1.0e20

#    Number of constraints (equalities plus inequalities)

     m <- 2

#    Lagrange multipliers approximation.

     lambda <- array(0.0, c(m))

#    For each constraint i, set equatn(i) = 1. if it is an equality
#    constraint of the form c_i(x) = 0, and set equatn(i) = 0 if
#    it is an inequality constraint of the form c_i(x) <= 0.

     equatn[1] <- 0
     equatn[2] <- 0

#    For each constraint i, set linear(i) = 1 if it is a linear
#    constraint, otherwise set linear(i) = 0

     linear[1] <- 0
     linear[2] <- 1

#    In this R interface evalf, evalg, evalh, evalc, evaljac and
#    evalhc are present. evalfc, evalgjac, evalhl and evalhlp are
#    present only for instructive purposes but they are not being
#    used at all.

     coded[1]  <- 1 # evalf
     coded[2]  <- 1 # evalg
     coded[3]  <- 1 # evalh
     coded[4]  <- 1 # evalc
     coded[5]  <- 1 # evaljac
     coded[6]  <- 1 # evalhc
     coded[7]  <- 0 # evalfc
     coded[8]  <- 0 # evalgjac
     coded[9]  <- 0 # evalhl
     coded[10] <- 0 # evalhlp

     checkder  <- 1

     }

#    ******************************************************************
#    ******************************************************************

     evalf <- function(n,x,f,flag) {

     flag <- 0

     f <- x[2]

     }

#    ******************************************************************
#    ******************************************************************

     evalg <- function(n,x,g,flag) {

     flag <- 0

     g[1] <- 0
     g[2] <- 1.0

     }

#    ******************************************************************
#    ******************************************************************

     evalh <- function(n,x,hlin,hcol,hval,hnnz,flag) {

     flag <- 0

     hnnz <- 0

     }

#    ******************************************************************
#    ******************************************************************

     evalc <- function(n,x,ind,cind,flag) {

     flag <- 0

     if (ind == 1)
         cind <- x[1] * x[1] + 1.0 - x[2]

     else if (ind == 2)
         cind <- 2.0 - x[1] - x[2]

     else
         flag <- -1

     }

#    ******************************************************************
#    ******************************************************************

     evaljac <- function(n,x,ind,jcvar,jcval,jcnnz,flag) {

     flag <- 0

     if ( ind == 1 ) {
         jcnnz <- 2

         jcvar[1] <- 1
         jcval[1] <- 2.0 * x[1]

         jcvar[2] <- 2
         jcval[2] <- - 1.0
     }

     else if ( ind == 2 ) {
         jcnnz <- 2

         jcvar[1] <- 1.0
         jcval[1] <- - 1.0

         jcvar[2] <- 2.0
         jcval[2] <- - 1.0
     }

     else
         flag <- -1

     }

#    ******************************************************************
#    ******************************************************************

     evalhc <- function(n,x,ind,hclin,hccol,hcval,hcnnz,flag) {

     flag <- 0

     if ( ind == 1 ) {
         hcnnz <- 1

         hclin[1] <- 1
         hccol[1] <- 1
         hcval[1] <- 2.0
     }

     else if ( ind == 2 )
         hcnnz <- 0

     else
         flag <- -1

     }

#    *****************************************************************
#    *****************************************************************

     evalfc <- function(n,x,f,m,constr,flag) {

     flag <- -1

     }

#    *****************************************************************
#    *****************************************************************

     evalgjac <- function(n,x,g,m,jcfun,jcvar,jcval,jcnnz,flag) {

     flag <- -1

     }

#    *****************************************************************
#    *****************************************************************

     evalhl <- function(n,x,m,lambda,sf,sc,hllin,hlcol,hlval,hlnnz,
     flag) {

     flag <- -1

     }


#    ******************************************************************
#    ******************************************************************
     evalhlp <- function(n,x,m,lambda,sf,sc,p,hp,gothl,flag) {

     flag <- -1

     }

#    ******************************************************************
#    ******************************************************************

     endp <- function(n,x,l,u,m,lambda,equatn,linear) {

     }

