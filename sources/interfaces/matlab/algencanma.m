evalf    = 'evalf';
evalg    = 'evalg';
evalh    = 'evalh';
evalc    = 'evalco';
evaljac  = 'evaljac';
evalhc   = 'evalhc';
evalfc   = 'evalfco';
evalgjac = 'evalgjac';
evalhl   = 'evalhl';
evalhlp  = 'evalhlp';

%% SET UP PROBLEM DATA
[n x l u m lambda equatn linear coded checkder] = inip;

%% SET SOME SOLVER ARGUMENTS
[epsfeas,epsopt,iprint,ncomp] = param;

%% CALL OPTIMIZATION SOLVER
[x,lambda,f,cnorm,snorm,nlpsupn,inform] =                   ...
    algencan(evalf,evalg,evalh,evalc,evaljac,evalhc,evalfc, ...
    evalgjac,evalhl,evalhlp,epsfeas,epsopt,iprint,ncomp,n,  ...
    x,l,u,m,lambda,equatn,linear,coded,checkder);

%% WRITE ADDITIONAL OUTPUT INFORMATION CODED BY THE USER
endp(n,x,l,u,m,lambda,equatn,linear);

