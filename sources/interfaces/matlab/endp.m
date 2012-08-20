function endp(n,x,l,u,m,lambda,equatn,linear)

%     [] = endp(n,x,l,u,m,lambda,equatn,linear)
%
%
%     This function can be used to do some extra job after the solver
%     has found the solution,like some extra statistics, or to save the
%     solution in some special format or to draw some graphical
%     representation of the solution. If the information given by the
%     solver is enough for you then leave the body of this subroutine
%     empty.
%
%     Parameters of the subroutine:
%
%     On Entry:
%
%     n        integer,
%              number of variables,
%
%     x        double precision x(n),
%              initial point,
%
%     l        double precision l(n),
%              lower bounds on x,
%
%     u        double precision u(n),
%              upper bounds on x,
%
%     m        integer,
%              number of constraints (excluding the bounds),
%
%     lambda   double precision lambda(m),
%              initial estimation of the Lagrange multipliers,
%
%     equatn   logical equatn(m)
%              for each constraint j, set equatn(j) = .true. if it is an
%              equality constraint of the form c_j(x) = 0, and set
%              equatn(j) = .false. if it is an inequality constraint of
%              the form c_j(x) <= 0,
%
%     linear   logical linear(m)
%              for each constraint j, set linear(j) = .true. if it is a
%              linear constraint, and set linear(j) = .false. if it is a
%              nonlinear constraint.


end
