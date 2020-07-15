function [y] = L122E(x)


x(isnan(x))=0;
y=abs(x).^2;
y(:,:,:,:,:,:,:,3)=2*y(:,:,:,:,:,:,:,3);
a=sum(y,8);

a=sum(a,7);

y=sqrt(a);

end
