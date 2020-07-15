function [y] = L122V(x)

x(isnan(x))=0;
y=abs(x).^2;

y=sum(y,8);
y=sum(y,7);
y=sqrt(y);

end
