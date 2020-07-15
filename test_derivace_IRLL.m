clear all
FA=5/180*pi;
tau=6;
tr=6;,
TR=tr;
td=6;
scale_2=2;
scale_1=1;
scale=100;

cos_phi=cos(FA);
sin_phi=sin(FA);
syms Efit n

Etau = Efit.^(tau/scale);
Etr = Efit.^(tr/scale);
Etd = Efit.^(td/scale);
M0 = 1;
M0_sc =scale_1;
F = (1 - Etau)./(1-Etau*cos_phi);
Q = (-Etr.*Etd.*F*(-(Etau.*cos_phi).^(800 - 1) + 1)*cos_phi + Etr.*Etd - 2*Etd + 1)./(Etr.*Etd.*(Etau.*cos_phi).^(800 - 1).*cos_phi + 1);
Q_F = Q-F;


mod=M0.*M0_sc.*(((Etau*cos_phi).^(n - 1)).*Q_F + F)*sin_phi;
testval=0.5:0.005:1;
out=zeros(size(testval));

out_d=zeros(size(testval));
d_M=diff(mod,Efit);
out_dg=zeros(size(testval));

tmp1 = ((-TR.*Etr.*Etd.*(-Etau + 1).*(-(Etau.*cos_phi).^(800 - 1) + 1)...
    .*cos_phi./(Efit.*scale.*(-Etau.*cos_phi + 1)) + TR.*Etr.*Etd./(Efit.*scale) - tau.*Etr.*Etau.*Etd.*(-Etau + 1).*...
    (-(Etau.*cos_phi).^(800 - 1) + 1).*cos_phi.^2./(Efit.*scale.*(-Etau.*cos_phi + 1).^2) + ...
    tau.*Etr.*Etau.*Etd.*(-(Etau.*cos_phi).^(800 - 1) + 1).*cos_phi./(Efit.*scale.*(-Etau.*cos_phi + 1)) + ...
    tau.*Etr.*Etd.*(Etau.*cos_phi).^(800 - 1).*(800 - 1).*(-Etau + 1).*cos_phi./(Efit.*scale.*(-Etau.*cos_phi + 1)) - ...
    td.*Etr.*Etd.*(-Etau + 1).*(-(Etau.*cos_phi).^(800 - 1) + 1).*cos_phi./(Efit.*scale.*(-Etau.*cos_phi + 1)) +...
    td.*Etr.*Etd./(Efit.*scale) - 2.*td.*Etd./(Efit.*scale))./(Etr.*Etd.*(Etau.*cos_phi).^(800 - 1).*cos_phi + 1) +...
    (-TR.*Etr.*Etd.*(Etau.*cos_phi).^(800 - 1).*cos_phi./(Efit.*scale) - tau.*Etr.*Etd.*(Etau.*cos_phi).^(800 - 1).*(800 - 1)...
    .*cos_phi./(Efit.*scale) - td.*Etr.*Etd.*(Etau.*cos_phi).^(800 - 1).*cos_phi./(Efit.*scale)).*(-Etr.*Etd.*(-Etau + 1)...
    .*(-(Etau.*cos_phi).^(800 - 1) + 1).*cos_phi./(-Etau.*cos_phi + 1) + Etr.*Etd - 2.*Etd + 1)./...
    (Etr.*Etd.*(Etau.*cos_phi).^(800 - 1).*cos_phi + 1).^2 - tau.*Etau.*(-Etau + 1).*cos_phi./(Efit.*scale.*(-Etau.*cos_phi + 1).^2) ...
    + tau.*Etau./(Efit.*scale.*(-Etau.*cos_phi + 1)));

tmp2 = tau.*Etau.*(-Etau + 1).*cos_phi./(Efit.*scale.*(-Etau.*cos_phi + 1).^2) - ...
    tau.*Etau./(Efit.*scale.*(-Etau.*cos_phi + 1));

tmp3 =  (-(-Etau + 1)./(-Etau.*cos_phi + 1) ...
    + (-Etr.*Etd.*(-Etau + 1).*(-(Etau.*cos_phi).^(800 - 1) + 1).*cos_phi./(-Etau.*cos_phi + 1) + Etr.*Etd - 2.*Etd + 1)...
    ./(Etr.*Etd.*(Etau.*cos_phi).^(800 - 1).*cos_phi + 1))./(Efit.*scale);

d_Mgraz=M0.*M0_sc.*((Etau.*cos_phi).^(n - 1).*tmp1 + tmp2 + tau.*(Etau.*cos_phi).^(n - 1).*(n - 1).*tmp3).*sin_phi;

for(i=1:length(testval))
    pom= (subs(mod,Efit,testval(i)));
    out(i)=double(subs(pom,n,i));
    pom=(subs(d_M,Efit,testval(i)));
    out_d(i)=double(subs(pom,n,i));
    pom=(subs(d_Mgraz,Efit,testval(i)));
    out_dg(i)=double(subs(pom,n,i));
    
    
end

disp(['RMS=' num2str(sqrt(sum((out_d-out_dg).^2,'all')))]);

