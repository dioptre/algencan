#   =================================================================
#   File: algencanma.r
#   =================================================================
#
#   =================================================================
#   Module: Main program
#   =================================================================
#
#   Last update of any of the component of this module:
#
#   February 14, 2006.
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
#   TANGO is free software; you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by
#   the Free Software Foundation. Non-free versions of TANGO are
#   available under terms different from those of the General Public
#   License. Professors J. M. Martínez (martinez@ime.unicamp.br,
#   martinezimecc@gmail.com) or E. G. Birgin (egbirgin@ime.usp.br,
#   egbirgin@gmail.com) should be contacted for more information
#   related to such a license, future developments and/or technical
#   support.
#
#   Every published work that uses ALGENCAN should cite:
#
#   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt,
#   "On Augmented Lagrangian methods with general lower-level
#   constraints", to appear in SIAM Journal on Optimization.
#
#   and
#
#   R. Andreani, E. G. Birgin, J. M. Martínez and M. L. Schuverdt,
#   "Augmented Lagrangian methods under the Constant Positive Linear
#   Dependence constraint qualification", Mathematical
#   Programming, 111, pp. 5-32, 2008.
#
#   Every published work that uses GENCAN should cite:
#
#   E. G. Birgin and J. M. Martínez, "Large-scale active-set
#   box-constrained optimization method with spectral projected
#   gradients", Computational Optimization and Applications 23, pp.
#   101-125, 2002.
#
#   (See other related works at the TANGO home page.)
#
#   ******************************************************************
#   ******************************************************************

    param <- function(epsfeas,epsopt,iprint,ncomp) {

#   ******************************************************************
#   FROM HERE ON YOU MUST MODIFY THE SUBROUTINE TO SET SOME ALGENCAN
#   ARGUMENTS RELATED TO DERIVATIVES, STOPPING CRITERIA AND OUTPUT:
#   ******************************************************************

    epsfeas <- 1.0e-08
    epsopt  <- 1.0e-08
    iprint  <- 10
    ncomp   <- 6


#   ******************************************************************
#   STOP HERE YOUR MODIFICATIONS OF SUBROUTINE PARAM.
#   ******************************************************************
}

#   Load the problem definition file

    source("toyprob.r")

#   Load the solver wrapper

    dyn.load("algencan.so")

#   Call the solver

   .Call("ralgencan",body(evalf),body(evalg),body(evalh),body(evalc),
    body(evaljac),body(evalhc),body(evalfc),body(evalgjac),body(evalhl),
    body(evalhlp),body(inip),body(endp),body(param),sys.frame(0))
