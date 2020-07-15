function [param] = setup_dcf(param)
 w = abs(linspace(-param.nFE./2,param.nFE./2,param.nFE))';
w = w.*pi./4./param.nSpokes;
w   = repmat(w, [1, param.nSpokes]);
w=sqrt(w);

param.w=w;
end

