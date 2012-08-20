function [c,flag] = evalco(n,x,ind)

% [c,flag] = evalco(n,x,ind)

  flag = 0;

  if (ind == 1)
    c = x(1)^2 + 1 - x(n);

  elseif (ind == 2)
    c = 2 - x(1) - x(n);

  else
    flag = - 1;
  end

end
