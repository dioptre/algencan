#   =================================================================
#   File: algencanma.tcl
#   =================================================================
#
#   =================================================================
#   Module: Main program
#   =================================================================
#
#   Last update of any of the component of this module: 
# 
#   April 28th, 2010.
#
#   Users are encouraged to download periodically updated versions of 
#   this code at the TANGO home page:
# 
#   www.ime.usp.br/~egbirgin/tango/
#
#   *****************************************************************
#   *****************************************************************
#
#   ALGENCAN solves problems of the form
#   ------------------------------------
#
#   min f(x)
#
#   subject to
#
#           c_j(x)  = 0, j in E,
#           c_j(x) <= 0, j in I,
#           l <= x <= u,
#
#   where E is the set of indices of the equality constraints, I is
#   the set of indices of the inequality constraints, and there are
#   n variables and m constraints.
#
#   ALGENCAN is an Augmented Lagrangian method that uses GENCAN to
#   solve the bound-constrained problems.
#
#   ALGENCAN is part of the TANGO Project.
#
#   Visit the TANGO home page in the web:
#
#   www.ime.usp.br/~egbirgin/tango/
#
#   *****************************************************************
#
#   TANGO LICENSE:
#   --------------
#
#   See file 'license.txt'.
#
#   ******************************************************************
#   ******************************************************************


proc param {pepsfeas pepsopt piprint pncomp} {

    upvar $pepsfeas epsfeas
    upvar $pepsopt   epsopt
    upvar $piprint   iprint
    upvar $pncomp     ncomp

    set epsfeas 1.0e-08	
    set epsopt  1.0e-08
    set iprint       10
    set ncomp         6

}

# LOAD THE PROBLEM DEFINITION FILE
source "toyprob.tcl"

# LOAD THE SOLVER WRAPPER
load ./tclwrapper[info sharedlibextension]

# SET SOME SOLVER PARAMETERS
param epsfeas epsopt iprint ncomp

# SET UP PROBLEM DATA
inip n x l u m lambda equatn linear coded checkder

# CALL OPTMIZATION SOLVER
algencan "evalf" "evalg" "evalh" "evalc" "evaljac" "evalhc" "evalfc" \
"evalgjac" "evalhl" "evalhlp" $epsfeas $epsopt $iprint $ncomp $n x l \
u $m lambda equatn linear coded $checkder f cnorm snorm nlpsupn inform

# WRITE ADDITTIONAL INFORMATION CODED BY THE USER
endp n x l u m lambda equatn linear
