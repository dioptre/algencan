#!/usr/bin/env python
# -*- coding: iso-8859-15 -*-

#   =================================================================
#   File: runalgencan
#   =================================================================
#
#   =================================================================
#   Module: Main program
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
#   *****************************************************************
#   *****************************************************************

#   Load the problem definition file

from toyprob import *

#   Import the Python module

import algencan

#   Set some optional params

algencan.param['iprint'] = 10

#   Call the solver

algencan.solvers(evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc,evalgjac,evalhl,
                 evalhlp,inip,endp)
