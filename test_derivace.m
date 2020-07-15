clear all
Model.param.Model_scale=ones(1,2);
sin_FA=sin(pi/4);
cos_FA=cos(pi/4);
syms M T;

f=(Model.param.Model_scale(1).*M.*sin_FA.*(1-Model.param.Model_scale(2).*T))./(1-cos_FA.*Model.param.Model_scale(2).*T);

d_M=diff(f,M);
d_T=diff(f,T);

out_syms=cat(2,(Model.param.Model_scale(1).*sin_FA.*(1-Model.param.Model_scale(2).*T))./(1-cos_FA.*Model.param.Model_scale(2).*T),Model.param.Model_scale(2).*Model.param.Model_scale(1).*(M.*(-sin_FA+sin_FA.*cos_FA))./(-2.*Model.param.Model_scale(2).*T.*cos_FA+(Model.param.Model_scale(2).*T).^2.*cos_FA.^2+1));

d_M_r=(Model.param.Model_scale(1).*sin_FA.*(1-Model.param.Model_scale(2).*T))./(1-cos_FA.*Model.param.Model_scale(2).*T);
d_T_r=Model.param.Model_scale(2)*Model.param.Model_scale(1)*(M.*(-sin_FA+sin_FA*cos_FA))./(1-Model.param.Model_scale(2).*T*cos_FA).^2;
testvalM=-100:100;

d_M_f=double(subs(d_M,T,testvalM));
d_M_fr=double(subs(d_M_r,T,testvalM));

M0=1;

E1=Model.param.Model_scale(2)*T;

grad_T1 = M*Model.param.Model_scale(2).*Model.param.Model_scale(1)*(-E1 + 1)*sin_FA*cos_FA./(-E1*cos_FA + 1).^2 -...
    M*Model.param.Model_scale(2).*Model.param.Model_scale(1)*sin_FA./(-E1*cos_FA + 1);


vys_T=subs(d_T,M,M0);
vys_Tr=subs(d_T_r,M,M0);
vys_graz=subs(grad_T1,M,M0);
vys_T=double(subs(vys_T,T,testvalM));
vys_Tr=double(subs(vys_Tr,T,testvalM));

vys_graz=double(subs(vys_graz,T,testvalM));