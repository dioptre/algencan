function [h,flag] = evalhc(n,x,ind)

% [h,flag] = evalhc(n,x,ind)

  flag = 0;

  if (ind == 1)
    h = spalloc(n,n,1);
    h(1,1) = 2;

  elseif (ind == 2)
    h = spalloc(n,n,0);

  else
    flag = - 1;
  end

end
