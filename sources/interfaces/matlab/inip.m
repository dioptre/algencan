function [n,x,l,u,m,lambda,equatn,linear,coded,checkder] = inip

% [n,x,l,u,m,lambda,equatn,linear,coded,checkder] = inip

% Number of variables

  n = 2;

% Initial point

  x = zeros(n-1,1);

  x(n) = 0;

% Lower and upper bounds

      for i = 1:n - 1
          l(i) = - 10;
          u(i) =   10;
      end

      l(n) = - 1.0d+20;
      u(n) =   1.0d+20;

%     Number of constraints (equalities plus inequalities)

      m = 2;

%     Lagrange multipliers approximation. Most users prefer to use the
%     null initial Lagrange multipliers estimates. However, if the
%     problem that you are solving is "slightly different" from a
%     previously solved problem of which you know the correct Lagrange
%     multipliers, we encourage you to set these multipliers as initial
%     estimates. Of course, in this case you are also encouraged to use
%     the solution of the previous problem as initial estimate of the
%     solution. Similarly, most users prefer to use rho = 10 as initial
%     penalty parameters. But in the case mentioned above (good
%     estimates of solution and Lagrange multipliers) larger values of
%     the penalty parameters (say, rho = 1000) may be more useful. More
%     warm-start procedures are being elaborated.

      lambda = zeros(m,1);

%     For each constraint i, set equatn(i) = .true. if it is an equality
%     constraint of the form c_i(x) = 0, and set equatn(i) = .false. if
%     it is an inequality constraint of the form c_i(x) <= 0.

      equatn = zeros(m,1);

%     For each constraint i, set linear(i) = .true. if it is a linear
%     constraint, otherwise set linear(i) = .false.

      linear(1) = 0;
      linear(2) = 1;

%     Indicate which subroutines did you code.

      coded( 1) = 1;  % evalf
      coded( 2) = 1;  % evalg
      coded( 3) = 1;  % evalh
      coded( 4) = 1;  % evalc
      coded( 5) = 1;  % evaljac
      coded( 6) = 1;  % evalhc
      coded( 7) = 0;  % evalfc
      coded( 8) = 0;  % evalgjac
      coded( 9) = 0;  % evalhl
      coded(10) = 0;  % evalhlp

%     Set checkder = TRUE if you code some derivatives and you would
%     like them to be tested by finite differences. It is highly
%     recommended.

      checkder = 1;

end
