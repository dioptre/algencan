function [jac,flag] = evaljac(n,x,ind)

% [jac,flag] = evaljac(n,x,ind)

  flag = 0;

  if (ind == 1)
    jac = spalloc(n,1,2);
    jac(1) = 2 * x(1);
    jac(n) = - 1;

  elseif (ind == 2)
    jac = spalloc(n,1,2);
    jac(1) = - 1;
    jac(n) = - 1;

  else
    flag = - 1;
  end

end
