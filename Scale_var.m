function [y] = Scale_var(param,u)

    y = u;
    for j=1:size(u,7)
      y(:,:,:,:,:,:,j) =y(:,:,:,:,:,:,j)./ param.ukscale(j);
      if j==1
        y(:,:,:,:,:,:,j) = y(:,:,:,:,:,:,j).*max(param.ukscale)./param.ratio;
      else
       y(:,:,:,:,:,:,j) = y(:,:,:,:,:,:,j).*max(param.ukscale);
      end
    end