#   =================================================================
#   File: toyprob.py
#   =================================================================
#
#   =================================================================
#   Module: Subroutines that define the problem
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

def inip():
    """This subroutine must set some problem data.

    For achieving this objective YOU MUST MODIFY it according to your
    problem. See below where your modifications must be inserted.

    Parameters of the subroutine:

    On Entry:

    This subroutine has no input parameters.

    On Return

    n        integer,
             number of variables,

    x        double precision x(n),
             initial point,

    l        double precision l(n),
             lower bounds on x,

    u        double precision u(n),
             upper bounds on x,

    m        integer,
             number of constraints (excluding the bounds),

    lambda   double precision lambda(m),
             initial estimation of the Lagrange multipliers,

    equatn   logical equatn(m)
             for each constraint j, set equatn(j) = .true. if it is an
             equality constraint of the form c_j(x) = 0, and set
             equatn(j) = .false. if it is an inequality constraint of
             the form c_j(x) <= 0,

    linear   logical linear(m)
             for each constraint j, set linear(j) = .true. if it is a
             linear constraint, and set linear(j) = .false. if it is a
             nonlinear constraint.
    """

#   Number of variables

    n = 2

#   Number of constraints (equalities plus inequalities)

    m = 2

#   Initial point

    x = zeros(n)

#   Lower and upper bounds

    l = array([-10.0,-1.0e20])

    u = array([ 10.0, 1.0e20])

#   Lagrange multipliers approximation. Most users prefer to use the
#   null initial Lagrange multipliers estimates. However, if the
#   problem that you are solving is "slightly different" from a
#   previously solved problem of which you know the correct Lagrange
#   multipliers, we encourage you to set these multipliers as initial
#   estimates. Of course, in this case you are also encouraged to use
#   the solution of the previous problem as initial estimate of the
#   solution. Similarly, most users prefer to use rho = 10 as initial
#   penalty parameters. But in the case mentioned above (good
#   estimates of solution and Lagrange multipliers) larger values of
#   the penalty parameters (say, rho = 1000) may be more useful. More
#   warm-start procedures are being elaborated.

    lambda_ = zeros(m)

#   For each constraint i, set equatn[i] = 1. if it is an equality
#   constraint of the form c_i(x) = 0, and set equatn[i] = 0 if
#   it is an inequality constraint of the form c_i(x) <= 0.

    equatn = [False, False]

#   For each constraint i, set linear[i] = 1 if it is a linear
#   constraint, otherwise set linear[i] = 0.

    linear = [False, True]

#   In this Python interface evalf, evalg, evalh, evalc, evaljac and
#   evalhc are present. evalfc, evalgjac, evalhl and evalhlp are not.

    coded = [True,  # evalf
             True,  # evalg
             True,  # evalh
             True,  # evalc
             True,  # evaljac
             True,  # evalhc
             False, # evalfc
             False, # evalgjac
             False, # evalhl
             False] # evalhlp

#   Set checkder = 1 if you code some derivatives and you would
#   like them to be tested by finite differences. It is highly
#   recommended.

    checkder = True

    return n,x,l,u,m,lambda_,equatn,linear,coded,checkder

#   ******************************************************************
#   ******************************************************************

def evalf(x):
    """This subroutine must compute the objective function.

    For achieving this objective YOU MUST MODIFY it according to your
    problem. See below where your modifications must be inserted.

    Parameters of the subroutine:

    On Entry:

    x        double precision x(n),
             current point,

    On Return

    f        double precision,
             objective function value at x,

    flag     integer,
             You must set it to any number different of 0 (zero) if
             some error ocurred during the evaluation of the objective
             function. (For example, trying to compute the square root
             of a negative number, dividing by zero or a very small
             number, etc.) If everything was o.k. you must set it
             equal to zero.
    """

    flag = 0

    n = len(x)

    f = x[n-1]

    return f,flag

#   ******************************************************************
#   ******************************************************************

def evalg(x):
    """This subroutine must compute the gradient vector of the objective
    function.

    For achieving these objective YOU MUST MODIFY it in the way specified
    below. However, if you decide to use numerical derivatives (we dont
    encourage this option at all!) you dont need to modify evalg.

    Parameters of the subroutine:

    On Entry:

    x        double precision x(n),
             current point,

    On Return

    g        double precision g(n),
             gradient vector of the objective function evaluated at x,

    flag     integer,
             You must set it to any number different of 0 (zero) if
             some error ocurred during the evaluation of any component
             of the gradient vector. (For example, trying to compute
             the square root of a negative number, dividing by zero or
             a very small number, etc.) If everything was o.k. you
             must set it equal to zero.
    """

    flag = 0

    n = len(x)

    g = zeros(n)

    g[n-1] = 1.0

    return g,flag

#   ******************************************************************
#   ******************************************************************

def evalh(x):
    """This subroutine might compute the Hessian matrix of the objective
    function.

    For achieving this objective YOU MAY MODIFY it according to your
    problem. To modify this subroutine IS NOT MANDATORY. See below
    where your modifications must be inserted.

    Parameters of the subroutine:

    On Entry:

    x        double precision x(n),
             current point,

    On Return

    hnnz     integer,
             number of perhaps-non-null elements of the computed
             Hessian,

    hlin     integer hlin(nnzh),
             see below,

    hcol     integer hcol(nnzh),
             see below,

    hval     double precision hval(nnzh),
             the non-null value of the (hlin(k),hcol(k)) position
             of the Hessian matrix of the objective function must
             be saved at hval(k). Just the lower triangular part of
             Hessian matrix must be computed,

    flag     integer,
             You must set it to any number different of 0 (zero) if
             some error ocurred during the evaluation of the Hessian
             matrix of the objective funtion. (For example, trying
             to compute the square root of a negative number,
             dividing by zero or a very small number, etc.) If
             everything was o.k. you must set it equal to zero.
    """

    flag = 0

    hnnz = 0

    hlin = zeros(hnnz, int)
    hcol = zeros(hnnz, int)
    hval = zeros(hnnz, float)

    return hlin,hcol,hval,hnnz,flag

#   ******************************************************************
#   ******************************************************************

def evalc(x,ind):
    """This subroutine must compute the ind-th constraint.

    For achieving this objective YOU MUST MOFIFY it according to your
    problem. See below the places where your modifications must be
    inserted.

    Parameters of the subroutine:

    On Entry:

    x        double precision x(n),
             current point,

    ind      integer,
             index of the constraint to be computed,

    On Return

    c        double precision,
             i-th constraint evaluated at x,

    flag     integer
             You must set it to any number different of 0 (zero) if
             some error ocurred during the evaluation of the
             constraint. (For example, trying to compute the square
             root of a negative number, dividing by zero or a very
             small number, etc.) If everything was o.k. you must set
             it equal to zero.
    """

    flag = 0

    n = len(x)

    if ind == 1:
        c = x[0] * x[0] + 1.0 - x[n-1]

    elif ind == 2:
        c = 2.0 - x[0] - x[n-1]

    else:
        flag = -1

    return c,flag

#   ******************************************************************
#   ******************************************************************

def evaljac(x,ind):
    """This subroutine must compute the gradient of the ind-th constraint.

    For achieving these objective YOU MUST MODIFY it in the way specified
    below.

    Parameters of the subroutine:

    On Entry:

    x        double precision x(n),
             current point,

    ind      integer,
             index of the constraint whose gradient will be computed,

    On Return

    jcnnz    integer,
             number of perhaps-non-null elements of the computed
             gradient,

    jcvar    integer jcvar(jcnnz),
             see below,

    jcval    double precision jcval(jcnnz),
             the non-null value of the partial derivative of the i-th
             constraint with respect to the jcvar(k)-th variable must
             be saved at jcval(k).

    flag     integer
             You must set it to any number different of 0 (zero) if
             some error ocurred during the evaluation of the
             constraint. (For example, trying to compute the square
             root of a negative number, dividing by zero or a very
             small number, etc.) If everything was o.k. you must set
             it equal to zero.
    """

    flag = 0

    n = len(x)

    if ind == 1:
        jcnnz = 2

        jcvar = zeros(jcnnz, int)
        jcval = zeros(jcnnz, float)

        jcvar[0] = 0
        jcval[0] = 2.0 * x[0]

        jcvar[1] = n - 1
        jcval[1] = - 1.0

    elif ind == 2:
        jcnnz = 2

        jcvar = zeros(jcnnz, int)
        jcval = zeros(jcnnz, float)

        jcvar[0] =  0
        jcval[0] = -1.0

        jcvar[1] = n - 1
        jcval[1] = - 1.0

    else:
        flag = -1

    return jcvar,jcval,jcnnz,flag

#   ******************************************************************
#   ******************************************************************

def evalhc(x,ind):
    """This subroutine might compute the Hessian matrix of the ind-th
    constraint.

    For achieving this objective YOU MAY MODIFY it according to your
    problem. To modify this subroutine IS NOT MANDATORY. See below
    where your modifications must be inserted.

    Parameters of the subroutine:

    On Entry:

    x        double precision x(n),
             current point,

    ind      integer,
             index of the constraint whose Hessian will be computed,

    On Return

    hcnnz    integer,
             number of perhaps-non-null elements of the computed
             Hessian,

    hclin    integer hclin(hcnnz),
             see below,

    hccol    integer hccol(hcnnz),
             see below,

    hcval    double precision hcval(hcnnz),
             the non-null value of the (hclin(k),hccol(k)) position
             of the Hessian matrix of the ind-th constraint must
             be saved at hcval(k). Just the lower triangular part of
             Hessian matrix must be computed,

    flag     integer,
             You must set it to any number different of 0 (zero) if
             some error ocurred during the evaluation of the Hessian
             matrix of the ind-th constraint. (For example, trying
             to compute the square root of a negative number,
             dividing by zero or a very small number, etc.) If
             everything was o.k. you must set it equal to zero.
    """

    flag = 0

    n = len(x)

    if ind == 1:
        hcnnz = 1

        hclin = zeros(hcnnz, int)
        hccol = zeros(hcnnz, int)
        hcval = zeros(hcnnz, float)

        hclin[0] = 0
        hccol[0] = 0
        hcval[0] = 2.0

    elif ind == 2:
        hcnnz = 0

        hclin = zeros(hcnnz, int)
        hccol = zeros(hcnnz, int)
        hcval = zeros(hcnnz, float)

    else:
        flag = -1

    return hclin,hccol,hcval,hcnnz,flag

#   ******************************************************************
#   ******************************************************************

def evalfc(x,m):

    flag = -1

    f = 0.0

    c = zeros(m)

    return f,c,flag

#   ******************************************************************
#   ******************************************************************

def evalgjac(x,m):

    flag = -1

    n = len(x)

    g = zeros(n)

    jcnnz = 0

    jcfun = zeros(jcnnz, int)
    jcvar = zeros(jcnnz, int)
    jcval = zeros(jcnnz, float)

    return g,jcfun,jcvar,jcval,jcnnz,flag

#   ******************************************************************
#   ******************************************************************

def evalhl(x,m,lambda_,scalef,scalec):

    flag = -1

    hlnnz = 0

    hllin = zeros(hlnnz, int)
    hlcol = zeros(hlnnz, int)
    hlval = zeros(hlnnz, float)

    return hllin,hlcol,hlval,hlnnz,flag

#   ******************************************************************
#   ******************************************************************

def evalhlp(x,m,lambda_,scalef,scalec,p,goth):

    flag = -1

    n = len(x)

    hp = zeros(n)

    return hp,goth,flag

#   ******************************************************************
#   ******************************************************************

def endp(x,l,u,m,lambda_,equatn,linear):
    """This subroutine can be used to do some extra job.

    This subroutine can be used to do some extra job after the solver
    has found the solution, like some extra statistics, or to save the
    solution in some special format or to draw some graphical
    representation of the solution. If the information given by the
    solver is enough for you then leave the body of this subroutine
    empty.

    Parameters of the subroutine:

    The parameters of this subroutine are the same parameters of
    subroutine inip. But in this subroutine there are not output
    parameter. All the parameters are input parameters.
    """

    pass
