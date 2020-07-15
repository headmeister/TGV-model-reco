classdef Model_IRLL
    
    properties
        param
        guess
        DA
        rez
        constraints
    end
    
    methods
        
        function [Model] =Model_IRLL(param,images)
            
            
            Model.param=param;
            Model.param.Model_scale=ones(1,2);
            Model.param.T1_sc=100;
            
            test_T1 = reshape(linspace(10,5500,numel(images(:,:,1,1,1,1))),[size(images,1),size(images,2)]);
            test_M0 = ones(size(images(:,:,1,1,1,1)));
            test_T1 = 1/Model.param.Model_scale(2)*exp(-Model.param.T1_sc./(test_T1));
            
            %  u=cat(4,test_M0*ones(size(images(:,:,1,1,1,1))),test_T1);
            
            G_x = forward(Model,(cat(7,test_M0./Model.param.Model_scale(1),test_T1)));
            
            I=abs(images(:));
            Model.param.Model_scale(1) =( Model.param.Model_scale(1)*median(I(I>0))/median(abs((G_x(:)))));
            
            
            DG_x =  gradient(Model,(cat(7,test_M0./Model.param.Model_scale(1),test_T1)));
            
            
            Model.param.Model_scale(2) =  Model.param.Model_scale(2)*(norm(reshape(abs(DG_x(:,:,:,:,:,:,1)),1,[]))/norm(reshape(abs(DG_x(:,:,:,:,:,:,2)),1,[])));
            
           % Model.param.Model_scale(2) =  Model.param.Model_scale(2)/sqrt(Model.param.Model_scale(1));
            
            %   DG_x =  self.execute_gradient_3D(np.array([test_M0/self.uk_sc[0]*np.ones((Nislice,dimY,dimX),dtype=DTYPE),test_T1],dtype=DTYPE))
            
            
            Model.guess = cat(7,1/Model.param.Model_scale(1)*ones(Model.param.size_Model_data),1/Model.param.Model_scale(2)*exp(-100./(2500*ones(param.size_Model_data))));
            
            
            Model.constraints=[-300,300;exp(-Model.param.T1_sc/(10))/Model.param.Model_scale(2),exp(-Model.param.T1_sc/(5500))/Model.param.Model_scale(2)];
            Model.param.main_p=param.FA;
            
        end
        
        
        function [out] = forward(Model,x) % Forward model estimatiton S(M,T1,(T2))
            
            scale=100;
            
            cos_phi=cos(Model.param.FA);
            sin_phi=sin(Model.param.FA);
            
            Efit = x(:,:,:,:,:,:,2)*Model.param.Model_scale(2);
            
            Etau = Efit.^(Model.param.tau/scale);
            Etr = Efit.^(Model.param.tr/scale);
            Etd = Efit.^(Model.param.td/scale);
            M0 = x(:,:,:,:,:,:,1);
            M0_sc = Model.param.Model_scale(1);
            F = (1 - Etau)./(1-Etau*cos_phi);
            Q = (-Etr.*Etd.*F*(-(Etau.*cos_phi).^(Model.param.NProjPerInv - 1) + 1)*cos_phi + Etr.*Etd - 2*Etd + 1)./(Etr.*Etd.*(Etau.*cos_phi).^(Model.param.NProjPerInv - 1).*cos_phi + 1);
            Q_F = Q-F;
            
            out=zeros([Model.param.size_Data_cart,1,Model.param.NSpoke]);
            
            for i= 1:Model.param.NProjPerInv/Model.param.NSpoke
                for(j=1:Model.param.NSpoke)
                    
                    n = (i-1)*Model.param.NSpoke+j;
                    out(:,:,:,:,:,i,:,j) = M0.*M0_sc.*(((Etau*cos_phi).^(n - 1)).*Q_F + F)*sin_phi;
                end
            end
            
            out=mean(out,8);
            
            out(isnan(out))=1e-20;
            out(isinf(out))=1e-20;
        end
        
        
        
        function [out] = gradient(Model,x) % Model gradient estimatiton: [dS/dM,dS/dT1...]
            out=zeros([Model.param.size_Data_cart_grad,Model.param.NSpoke]);
            sin_phi=sin(Model.param.FA);
            M0_sc = (Model.param.Model_scale(1));
            cos_phi=(cos(Model.param.FA));
            scale=100;
            Efit = (x(:,:,:,:,:,:,2).*Model.param.Model_scale(2));
            Etau = Efit.^(Model.param.tau./scale);
            Etr = Efit.^(Model.param.tr./scale);
            Etd = Efit.^(Model.param.td./scale);
            
            M0 = (x(:,:,:,:,:,:,1));
            TR=(Model.param.tr);
            tau=(Model.param.tau);
            td=(Model.param.td);
            
            F = (1 - Etau)./(1-Etau.*cos_phi);
            Q = (-Etr.*Etd.*(-Etau + 1).*(-(Etau.*cos_phi).^(Model.param.NProjPerInv - 1) + 1).*cos_phi./(-Etau.*cos_phi + 1) + Etr.*Etd - 2.*Etd + 1)./(Etr.*Etd.*(Etau.*cos_phi).^(Model.param.NProjPerInv - 1).*cos_phi + 1);
            Q_F = Q-F;
            tmp1 = ((-TR.*Etr.*Etd.*(-Etau + 1).*(-(Etau.*cos_phi).^(Model.param.NProjPerInv - 1) + 1)...
                .*cos_phi./(x(:,:,:,:,:,:,2).*scale.*(-Etau.*cos_phi + 1)) + TR.*Etr.*Etd./(x(:,:,:,:,:,:,2).*scale) - tau.*Etr.*Etau.*Etd.*(-Etau + 1).*...
                (-(Etau.*cos_phi).^(Model.param.NProjPerInv - 1) + 1).*cos_phi.^2./(x(:,:,:,:,:,:,2).*scale.*(-Etau.*cos_phi + 1).^2) + ...
                tau.*Etr.*Etau.*Etd.*(-(Etau.*cos_phi).^(Model.param.NProjPerInv - 1) + 1).*cos_phi./(x(:,:,:,:,:,:,2).*scale.*(-Etau.*cos_phi + 1)) + ...
                tau.*Etr.*Etd.*(Etau.*cos_phi).^(Model.param.NProjPerInv - 1).*(Model.param.NProjPerInv - 1).*(-Etau + 1).*cos_phi./(x(:,:,:,:,:,:,2).*scale.*(-Etau.*cos_phi + 1)) - ...
                td.*Etr.*Etd.*(-Etau + 1).*(-(Etau.*cos_phi).^(Model.param.NProjPerInv - 1) + 1).*cos_phi./(x(:,:,:,:,:,:,2).*scale.*(-Etau.*cos_phi + 1)) +...
                td.*Etr.*Etd./(x(:,:,:,:,:,:,2).*scale) - 2.*td.*Etd./(x(:,:,:,:,:,:,2).*scale))./(Etr.*Etd.*(Etau.*cos_phi).^(Model.param.NProjPerInv - 1).*cos_phi + 1) +...
                (-TR.*Etr.*Etd.*(Etau.*cos_phi).^(Model.param.NProjPerInv - 1).*cos_phi./(x(:,:,:,:,:,:,2).*scale) - tau.*Etr.*Etd.*(Etau.*cos_phi).^(Model.param.NProjPerInv - 1).*(Model.param.NProjPerInv - 1)...
                .*cos_phi./(x(:,:,:,:,:,:,2).*scale) - td.*Etr.*Etd.*(Etau.*cos_phi).^(Model.param.NProjPerInv - 1).*cos_phi./(x(:,:,:,:,:,:,2).*scale)).*(-Etr.*Etd.*(-Etau + 1)...
                .*(-(Etau.*cos_phi).^(Model.param.NProjPerInv - 1) + 1).*cos_phi./(-Etau.*cos_phi + 1) + Etr.*Etd - 2.*Etd + 1)./...
                (Etr.*Etd.*(Etau.*cos_phi).^(Model.param.NProjPerInv - 1).*cos_phi + 1).^2 - tau.*Etau.*(-Etau + 1).*cos_phi./(x(:,:,:,:,:,:,2).*scale.*(-Etau.*cos_phi + 1).^2) ...
                + tau.*Etau./(x(:,:,:,:,:,:,2).*scale.*(-Etau.*cos_phi + 1)));
            
            tmp2 = tau.*Etau.*(-Etau + 1).*cos_phi./(x(:,:,:,:,:,:,2).*scale.*(-Etau.*cos_phi + 1).^2) - ...
                tau.*Etau./(x(:,:,:,:,:,:,2).*scale.*(-Etau.*cos_phi + 1));
            
            tmp3 =  (-(-Etau + 1)./(-Etau.*cos_phi + 1) ...
                + (-Etr.*Etd.*(-Etau + 1).*(-(Etau.*cos_phi).^(Model.param.NProjPerInv - 1) + 1).*cos_phi./(-Etau.*cos_phi + 1) + Etr.*Etd - 2.*Etd + 1)...
                ./(Etr.*Etd.*(Etau.*cos_phi).^(Model.param.NProjPerInv - 1).*cos_phi + 1))./(x(:,:,:,:,:,:,2).*scale);
            
            for i= 1:Model.param.NProjPerInv/Model.param.NSpoke
                for(j=1:Model.param.NSpoke)
                    n = (i-1).*Model.param.NSpoke+j;
                    
                   out(:,:,:,:,:,i,1,j) =( M0_sc.*((Etau.*cos_phi).^(n - 1).*Q_F + F).*sin_phi);
                    
                   out(:,:,:,:,:,i,2,j) =((M0.*M0_sc.*(Etau.*cos_phi).^(n - 1)).*tmp1 + tmp2 + (tau.*(Etau.*cos_phi).^(n - 1)).*(n - 1).*tmp3).*sin_phi;
                    
                    
                    
                    
                end
            end
            
            out=(mean(out,8));
            out(isnan(out))=1e-20;
            out(isinf(out))=1e-20;
        end
        
        function [Model] =update_DA(Model,x)
            
            Model.DA=gradient(Model,x);
            
        end
        
        function [Model] =update_constraints(Model,sc,index)
            
            
            Model.constraints(index,:)=Model.constraints(index,:)/sc;
            
        end
        
        
        function [Model,x] =update_scale(Model,x)
            
            Model=update_DA(Model,x);
            for(i=1:Model.param.nVar-1)
                sc=(norm(abs(reshape(Model.DA(:,:,:,:,:,:,1),1,[])))/norm(abs(reshape(Model.DA(:,:,:,:,:,:,i+1),1,[]))));
                
                x(:,:,:,:,:,:,i+1)=x(:,:,:,:,:,:,i+1)*Model.param.Model_scale(i+1);
                Model.param.Model_scale(i+1)=Model.param.Model_scale(i+1)*sc;
                x(:,:,:,:,:,:,i+1)=x(:,:,:,:,:,:,i+1)/Model.param.Model_scale(i+1);
                Model=update_constraints(Model,sc,i+1);
                p=x(:,:,:,:,:,:,i+1);
                p(p==Model.constraints(i+1,1))=1./Model.param.Model_scale(2).*exp(-Model.param.T1_sc./(1000));
                p(p==Model.constraints(i+1,2))=1./Model.param.Model_scale(2).*exp(-Model.param.T1_sc./(1000));
                x(:,:,:,:,:,:,i+1)=p;
            end
            
            Model=update_DA(Model,x);
            
            
            
            
        end
        
        function [x] =apply_constraints(Model,x)
            
            
            for(i=1:Model.param.nVar)
                
                pom=x(:,:,:,:,:,:,i);
                if(i==2 )
                    pom=real(pom);
                end
                
                pom(pom<Model.constraints(i,1))=Model.constraints(i,1);
                pom(pom>Model.constraints(i,2))=Model.constraints(i,2);
                x(:,:,:,:,:,:,i)=pom;
            end
            
        end
        
        
        
    end
    
end