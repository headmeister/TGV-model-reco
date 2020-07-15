function [param] = set_scale(param,u)

    for j = 1:size(u,7)
      param.ukscale(j) = norm(reshape(u(:,:,:,:,:,:,j),1,[]));
   
    end
    
    
    
  
end

