function [sens_map] = mask(velikost,bc,param)

A=zeros(velikost);

[X,Y]=meshgrid(1:velikost,1:velikost);

A=(X-velikost/2+10).^2+(Y-velikost/2+10).^2<=(0.75*velikost/2)^2;

A=imgaussfilt(double(A),4);

A=A+1e-12;
A=A./max(max(A));

 sens_map=bc.*repmat(A,1,1,param.nCh);
% 
% 
% sens_map_scale = max(abs(sens_map(:)));
% sens_map = sens_map/sens_map_scale;
% sens_map_conj = conj(sens_map);
% 
% sens_correct_term = 1./sum(sens_map_conj.*sens_map,4);
% 
% sens_correct_term = sqrt(sens_correct_term);
% sens_map = bsxfun(@times,sens_correct_term,sens_map);
% 
% sens_map = single(sens_map);


end

