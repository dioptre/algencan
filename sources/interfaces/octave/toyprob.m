#   =================================================================
#   File: toyprob.m
#   =================================================================
#   Coded by: Ricardo Andrade
#
#   =================================================================
#   Module: Subroutines that define the problem
#   =================================================================
#
#   Last update of any of the component of this module:
#
#   September 11, 2009.
#
#   Users are encouraged to download periodically updated versions of
#   this code at the TANGO home page:
#
#   www.ime.usp.br/~egbirgin/tango/
#
#   ******************************************************************
#   ******************************************************************

warning("off","Octave:undefined-return-values")

function [n,x,l,u,m,lambda,equatn,linear,coded,checkder]=inip()

#   Number of variables

    n = 2;

#   Initial point

    x = [0.0,0.0];

#   Lower and upper bounds

    l = [-10.0,-1.0e20];

    u = [ 10.0, 1.0e20];

#   Number of constraints (equalities plus inequalities)

    m = 2;

#   Lagrange multipliers approximation.

    lambda = zeros(1,m);

#   For each constraint i, set equatn(i) = 1. if it is an equality
#   constraint of the form c_i(x) = 0, and set equatn(i) = 0 if
#   it is an inequality constraint of the form c_i(x) <= 0.

    equatn = [0,0];

#   For each constraint i, set linear(i) = 1 if it is a linear
#   constraint, otherwise set linear(i) = 0

    linear = [0,1];

#   In this octave interface evalf, evalg, evalh, evalc, evaljac
#   and evalhc are present. evalfc, evalgjac, evalhl and evalhlp
#   are present only for instructive purposes but they are not
#   being used at all.

    coded(1)  = 1; # evalf
    coded(2)  = 1; # evalg
    coded(3)  = 1; # evalh
    coded(4)  = 1; # evalc
    coded(5)  = 1; # evaljac
    coded(6)  = 1; # evalhc
    coded(7)  = 0; # evalfc
    coded(8)  = 0; # evalgjac
    coded(9)  = 0; # evalhl
    coded(10) = 0; # evalhlp

    checkder  = 1;

endfunction

#   ******************************************************************
#   ******************************************************************

function [f,flag]=evalf(n,x)

    flag = 0;

    f = x(n);

endfunction

#   ******************************************************************
#   ******************************************************************

function [g,flag]=evalg(n,x)

    flag = 0;

    g = zeros(1,n);

    g(n) = 1.0;

endfunction

#   ******************************************************************
#   ******************************************************************

function [hlin,hcol,hval,hnnz,flag]=evalh(n,x)

    flag = 0;

    hnnz = 0;

endfunction

#   ******************************************************************
#   ******************************************************************

function [c,flag]=evalc(n,x,ind)

    flag = 0;

    if (ind == 1)
        c = x(1) * x(1) + 1.0 - x(n);

    elseif (ind == 2)
        c = 2.0 - x(1) - x(n);

    else
        flag = -1;
    endif

endfunction

#   ******************************************************************
#   ******************************************************************

function [jcvar,jcval,jcnnz,flag]=evaljac(n,x,ind)

    flag = 0;

    if (ind == 1)
        jcnnz = 2;

        jcvar(1) = 1;
        jcval(1) = 2.0 * x(1);

        jcvar(2) = n;
        jcval(2) = -1.0;

    elseif (ind == 2)
        jcnnz = 2;

        jcvar(1) =  1;
        jcval(1) = -1.0;

        jcvar(2) =  n;
        jcval(2) = -1.0;

    else
        flag = -1;
    endif

endfunction

#   ******************************************************************
#   ******************************************************************

function [hclin,hccol,hcval,hcnnz,flag]=evalhc(n,x,ind)

    flag = 0;

    if (ind == 1)
        hcnnz = 1;

        hclin = [];
        hccol = [];
        hcval = [];

        hclin(1) = 1;
        hccol(1) = 1;
        hcval(1) = 2.0;

    elseif (ind == 2)
        hcnnz = 0;

    else
        flag = -1;
    endif

endfunction

#   ******************************************************************
#   ******************************************************************

function [f,constr,flag]=evalfc(n,x,m)

    flag = -1;

endfunction

#   ******************************************************************
#   ******************************************************************

function [g,jcfun,jcvar,jcval,jcnnz,flag]=evalgjac(n,x,m)

    flag = -1;

endfunction

#   ******************************************************************
#   ******************************************************************

function [hllin,hlcol,hlval,hlnnz,flag]=evalhl(n,x,m,lambda,sf,sc)

    flag = -1;

endfunction


#   ******************************************************************
#   ******************************************************************

function [hp,gothl,flag]=evalhlp(n,x,m,lambda,sf,sc,p,gothl)

    flag = -1;

endfunction

#   ******************************************************************
#   ******************************************************************

function endp(n,x,l,u,m,lambda,equatn,linear)

endfunction
