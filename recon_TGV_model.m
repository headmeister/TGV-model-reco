%
clear all
close all
addpath(genpath(pwd));

% % load and reformat data
% o(1)= load('data1');
% o(2)= load('data2');
% o(3)= load('data3');
% o(4)= load('data5');
% o(5)= load('data8');
% o(6)= load('data11');
% o(7)= load('data19');
% 
% 
% for(i=1:length(o))
%     
%     rawdata(:,:,:,i)= permute(o(i).data.signal,[2,3,1]);
%     k(:,:,i)=o(i).data.k;
% end
% vel=size(rawdata);
% vel=[vel(1:2),1,vel(3:end)];
% rawdata=reshape(rawdata,vel);
% 
% k=reshape(k,size(k,1).*size(k,2),size(k,3));
% 
% kx=real(k);
% ky=imag(k);
% 
% kx=kx./max(max(abs(kx))).*0.5;
% ky=ky./max(max(abs(ky))).*0.5;
% 
% kx=reshape(kx,[1,size(kx)]);
% ky=reshape(ky,[1,size(ky)]);
% k=[kx;ky];
% 
% clear kx ky
% 
% 
% 
% n=128;
% k=k(:,1:n.*size(rawdata,1),:);
% 
% rawdata=rawdata(:,1:n,:,:,:,:);
% nSlice=1;
% 
% param = setup_param(rawdata,nSlice);
% param=setup_dcf(param);
% 
% 
% param.k=k;
% clear k;

% setup IRLL
load('dataIRLL');
 nSlice=1;

k=data.k;
rawdata=permute(data.signal,[2,3,1]);
%k=reshape(k,size(k,1).*size(k,2),size(k,3));
% % load('data_moz.mat');
% %
% % load('k.mat');
% %
% % rawdata=permute(data,[1,4,2,3,5]);
% rawdata=reshape(rawdata,size(rawdata,1),size(rawdata,2),size(rawdata,3),size(rawdata,4),1,size(rawdata,5));
% rawdata=rawdata(:,:,:,24,:,:,:);
kx=real(k);
ky=imag(k);

kx=kx./max(max(abs(kx))).*0.5;
ky=ky./max(max(abs(ky))).*0.5;



k=kx + 1j.*ky;




n=800;
 kx=kx(:,1:n,:);
 ky=ky(:,1:n,:);


rawdata=rawdata(:,1:n,:,:,:,:);
param.NSpoke=5;
rawdata=reshape(rawdata,size(rawdata,1),param.NSpoke,1,1,size(rawdata,2)./param.NSpoke);

param = setup_param(rawdata,nSlice);

param.NSpoke=5;
param.type='I';
param=setup_dcf(param);

param.tau=data.tau;
param.tr=5;
param.td=data.td;
param.NProjPerInv=800;

param.FA=data.FA;


kx=reshape(kx,[1,size(kx,1).*param.NSpoke,size(rawdata,5)]);
ky=reshape(ky,[1,size(ky,1).*param.NSpoke,size(rawdata,5)]);
k=[kx;ky];



clear kx ky
param.k=k;
 clear k;
%% setup Model and scale data
p=size(rawdata);
p(1:2)=1;
rawdata=rawdata.*repmat(sqrt(param.w),p);
rawdata=rawdata./norm(col(rawdata)).*sqrt(2000)*2;

%
%  param.FA=[1,2,3,5,8,11,19.47]./180.*pi;
% param.TR=6.56;


images=prep_images(param,rawdata);
obr=sum(images,[4,6,7,8]);

bc=get_sens_map(obr,'2D');
bc=squeeze(bc);

 bc=mask(param.nFE,bc,param);
bc_p=reshape(bc,[size(bc,1),size(bc,2),1,1,size(bc,3)]);
images=sum(images.*conj(repmat(bc_p,1,1,1,1,1,param.nPar)),5);

param.coil=bc;
clear bc;

Model=Model_IRLL(param,images);



%% initiate recon variables


uk=(Model.guess);

u=(zeros(size(uk)));
u_new=u;
z1=(zeros([param.size_Model_data_grad 2]));
v=z1;
v_new=v;
z1_new=z1;


z2=(zeros([param.size_Model_data_grad 3]));
r=(zeros(param.size_kSpace));
r_new=r;
z2_new=z2;
IRGN=setup_IRGN();

operator=operator(param);
index=1;
for(i=1:IRGN.GNsteps)
    
    
    
    [Model,uk]=update_scale(Model,uk);
    
    operator=calc_rez(operator,Model,uk,rawdata);
    param=set_scale(param,uk);
    tau=IRGN.tau_init;
    theta=IRGN.theta_init;
    beta=IRGN.beta;
    
    %---------------------------------------
    %   Kh=@(x,y,z)(cat(9,adjoint_gradient(operator,Model,x)+Scale_var(param,param.div1(y)),-y+param.div2(z)));
    
    
    Kh1=@(x,y) adjoint_gradient(operator,Model,x)+Scale_var(param,param.div1(y));
    Kh1_old=Kh1(r,z1);
    
    Kh2=@(y,z) -y+param.div2(z);
    Kh2_old=Kh2(z1,z2);
    
    u=uk;
    gap_min=0;
    primal=0;
    Axold = forward_gradient(operator,Model,u);
    
    for(j=1:IRGN.PDsteps)
        
        
        Pg=@(x) (((tau./IRGN.gamma).*uk+x)./(1+tau./IRGN.gamma));
        Palfa=@(x,alfa) (x./max(1,repmat(L122V(x)./(IRGN.lambda.*alfa),1,1,1,1,1,1,size(x,7),size(x,8))));
        Palfa1=@(x,alfa)(x./max(1,repmat(L122E(x)./(IRGN.lambda.*alfa),1,1,1,1,1,1,size(x,7),size(x,8))));
        
        
        
        u_new=Pg(u-tau.*(Kh1_old));
        v_new=v-tau.*(Kh2_old);
        u_new=apply_constraints(Model,u_new);
        
        beta_new = beta.*(1+tau./IRGN.gamma);
        tau_new = tau.*sqrt(beta./beta_new.*(1+theta));
        beta = beta_new;
        if(j>1)
            IRGN.PDsteps=200;
        end
        
        gradx = param.grad(Scale_var(param,u_new));
        gradx_xold = gradx - param.grad(Scale_var(param,u));
        v_vold = v_new-v;
        symgrad_v = param.E(v_new);
        symgrad_v_vold = symgrad_v - param.E(v);
        
        Ax = forward_gradient(operator,Model,u_new);
        Ax_Axold = Ax-Axold;
        
        while(1)
            theta=tau_new./tau;
            sigma=tau_new.*beta;
            PL2=@(x) ((x-sigma.*operator.rez)./(1+sigma/10));
            
            
            
            
            
            z1_new=Palfa(gradx + theta.*gradx_xold- v_new - theta.*v_vold ,IRGN.alfa1);
            z2_new=Palfa1(z2+beta.*tau_new.*(symgrad_v + theta.*symgrad_v_vold ),IRGN.alfa0);
            
            
            tmp = Ax+theta.*Ax_Axold;
            
            r_new=PL2(r+beta.*tau_new.*(tmp));
            
            Kh1_new=Kh1(r_new,z1_new);
            Kh2_new=Kh2(z1_new,z2_new);
            
            
            
            Ls=norm([col(Kh1_new);col(Kh2_new)]-[col(Kh1_old);col(Kh2_old)]);
            
            Rs=(norm([r_new(:)-r(:);z1_new(:)-z1(:);z2_new(:)-z2(:)]));
            
            if(sqrt(beta).*tau_new.*Ls<=IRGN.delta.*Rs)
                break
            else
                tau_new=tau_new.*IRGN.mu;
            end
            
            
            
            
        end
        
        primal_new= real(1./2.*norm(col(Ax-operator.rez))^2+IRGN.alfa1.*IRGN.lambda.*sum(abs((gradx-v_new)),'all') +...
            IRGN.alfa0.*IRGN.lambda.*sum(abs(symgrad_v),'all') + 1./(2.*IRGN.gamma).*norm(col(u_new-uk))^2);
        
        dual = real(-1./2.*norm(col(-Kh1_new))^2 - sum(col(uk).*col(-Kh1_new),'all') + sum(Kh2_new,'all')...
            - 1./(2.*IRGN.gamma).*norm(col(r_new))^2 - sum(col(operator.rez).*col(r_new)));
        
        gap(index) = abs(primal_new - dual);
        
        index=index+1;
        u=u_new;
        v=v_new;
        z1=z1_new;
        Kh1_old=Kh1_new;
        Kh2_old=Kh2_new;
        z2=z2_new;
        r=r_new;
        tau=tau_new;
        Axold=Ax;
        
        subplot(411);
        %         plot(log10(primal_new));
        %         hold on
        %         plot(log10(dual));
        
        plot(log10(gap));
        P=abs(u);
        subplot(412);
        
        imshow(P(:,:,:,:,:,:,1).*Model.param.Model_scale(1),[]);
        
        subplot(413);
        imagesc(-Model.param.T1_sc./log(Model.param.Model_scale(2).*(P(:,:,:,:,:,:,2))));
        daspect([1 1 1]);
        colorbar
        
        pause(0.001);
        
        [primal,finish,gap_min]=check_PD_conv(primal,primal_new,dual,gap_min,j,IRGN);
        
        
        f_val(index)= (1/param.nFE^2)*((1/2*norm(col(rawdata - forward_forward(operator,Model,u))))^2+IRGN.lambda*sum(abs(param.grad(Scale_var(param,u)-v)),'all') +IRGN.lambda*(2)*sum(abs(param.E(v)),'all')+1/(2*IRGN.gamma)*norm(col(u-uk))^2);
        subplot(414);
        %         plot(log10(primal_new));
        %         hold on
        %         plot(log10(dual));
        
        plot(log10(f_val));
        
        if(finish==1)
            break
        end
        
    end
    uk=u;
    pom=Scale_var(param,u);
    scale=norm(reshape(param.grad((pom(:,:,:,:,:,:,1))),1,[]))./norm(reshape(param.grad((pom(:,:,:,:,:,:,2))),1,[]));
    
    if scale == 0 || isinf(scale)
        param.ratio = param.ratio;
    else
        param.ratio = scale.*param.ratio;
    end
    
    if(IRGN.gamma*IRGN.qgamma<=IRGN.gamma_max)
        IRGN.gamma=IRGN.gamma.*IRGN.qgamma;
    end
    
    if(IRGN.lambda*IRGN.qlambda>=IRGN.lambda_min)
        IRGN.lambda=IRGN.lambda.*IRGN.qlambda;
    end
    
end









