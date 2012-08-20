#   =================================================================
#   File: algencan.py
#   =================================================================
#
#   =================================================================
#   Module: Python module
#   =================================================================
#
#   Last update of any of the components of this module:
#
#   March 25, 2008.
#
#   Users are encouraged to download periodically updated versions of
#   this code at the TANGO home page:
#
#   www.ime.usp.br/~egbirgin/tango/
#
#   ******************************************************************
#   ******************************************************************

from numpy import *

#   Load the solver wrapper

from pywrapper import solver

def solvers(evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc,evalgjac,evalhl,
            evalhlp,inip,endp):
    """Call the solver."""

    solver(evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc,evalgjac,evalhl,
           evalhlp,inip,endp,param)

#   Parameters of the solver:
#   =========================

param = {
    'epsfeas': 1.0e-08,
    'epsopt' : 1.0e-08,

    'iprint': 10,
    'ncomp' : 6
}
