classdef Model
    
    properties
        param
        guess
        DA
        rez
        constraints
    end
    
    methods
        
        function [Model] =Model(param,images)
            
            
            Model.param=param;
            Model.param.Model_scale=ones(1,2);
            
            
            test_T1 = reshape(linspace(10,6000,numel(images(:,:,1,1,1,1,1))),[size(images,1),size(images,2)]);
            test_M0 = sqrt((size(images,1).*pi./2)./param.nSpokes).*ones(size(images(:,:,1,1,1,1)));
            test_T1 = 1./Model.param.Model_scale(2).*exp(-param.TR./(test_T1));
            
            %  u=cat(4,test_M0.*ones(size(images(:,:,1,1,1,1))),test_T1);
            
            G_x = forward(Model,(cat(7,test_M0./Model.param.Model_scale(1),test_T1)));
            
            I=abs(images(:));
            
            Model.param.Model_scale(1) = Model.param.Model_scale(1).*median(I(I>0))./median(abs((G_x(:))));
             Model.param.Model_scale(1)= Model.param.Model_scale(1);
            
            
            DG_x =  gradient(Model,(cat(7,test_M0./Model.param.Model_scale(1),test_T1)));
            
            
            Model.param.Model_scale(2) =  Model.param.Model_scale(2).*(norm(reshape(abs(DG_x(:,:,:,:,:,:,1)),1,[]))./norm(reshape(abs(DG_x(:,:,:,:,:,:,2)),1,[])));
            
            Model.param.Model_scale(2) =  Model.param.Model_scale(2)./sqrt(Model.param.Model_scale(1));
         %   Model.param.Model_scale(2) = Model.param.Model_scale(2);
            %DG_x =  self.execute_gradient_3D(np.array([test_M0./self.uk_sc[0].*np.ones((Nislice,dimY,dimX),dtype=DTYPE),test_T1],dtype=DTYPE))
            
            
            Model.guess = cat(7,1./Model.param.Model_scale(1).*ones(Model.param.size_Model_data),1./Model.param.Model_scale(2).*exp(-param.TR./(500.*ones(param.size_Model_data))));
            
            
            Model.constraints=[-300./Model.param.Model_scale(1),300./Model.param.Model_scale(1);exp(-Model.param.TR./(5))./Model.param.Model_scale(2),exp(-Model.param.TR./(6000))./Model.param.Model_scale(2)];
            Model.param.main_p=param.FA;
            Model.param.T1_sc=param.TR;
            
        end
        
        
        function [out] = forward(Model,x) % Forward model estimatiton S(M,T1,(T2))
            
            sin_FA=sin(Model.param.FA);
            
            cos_FA=cos(Model.param.FA);
            out=(zeros(Model.param.size_Model_data));
            
            for(i=1:length(Model.param.FA))
                out(:,:,:,:,:,i)=(Model.param.Model_scale(1).*x(:,:,:,:,:,:,1).*sin_FA(i).*(1-Model.param.Model_scale(2).*x(:,:,:,:,:,:,2)))./(1-cos_FA(i).*Model.param.Model_scale(2).*x(:,:,:,:,:,:,2));
            end
            out(isnan(out))=1e-20;
            out(isinf(out))=1e-20;
        end
        
        
        
        function [out] = gradient(Model,x) % Model gradient estimatiton: [dS./dM,dS./dT1...]
            
            sin_FA=sin(Model.param.FA);
            
            cos_FA=cos(Model.param.FA);
            
            out=(zeros(Model.param.size_Model_data_grad));
            
            E1=Model.param.Model_scale(2).*x(:,:,:,:,:,:,2);


            
            for(i=1:length(Model.param.FA))
                out(:,:,:,:,:,i,:)=cat(7,(Model.param.Model_scale(1).*sin_FA(i).*(1-Model.param.Model_scale(2).*x(:,:,:,:,:,:,2)))./(1-cos_FA(i).*Model.param.Model_scale(2).*x(:,:,:,:,:,:,2)), (x(:,:,:,:,:,:,1).*Model.param.Model_scale(2).*Model.param.Model_scale(1).*(-E1 + 1).*sin_FA(i).*cos_FA(i))./(-E1.*cos_FA(i) + 1).^2 -...
    (x(:,:,:,:,:,:,1).*Model.param.Model_scale(2).*Model.param.Model_scale(1).*sin_FA(i))./(-E1.*cos_FA(i) + 1));
            end
            out(isnan(out))=1e-20;
            out(isinf(out))=1e-20;
        end
        
        function [Model] =update_DA(Model,x)
            
            Model.DA=gradient(Model,x);
            
        end
        
        function [Model] =update_constraints(Model,sc,index)
            
            
                Model.constraints(index,:)=Model.constraints(index,:)./sc;
            
        end
        
        
        function [Model,x] =update_scale(Model,x)
            
            Model=update_DA(Model,x);
            for(i=1:Model.param.nVar-1)
                sc=(norm(abs(reshape(Model.DA(:,:,:,:,:,:,1),1,[])))./norm(abs(reshape(Model.DA(:,:,:,:,:,:,i+1),1,[]))));
                
                x(:,:,:,:,:,:,i+1)=x(:,:,:,:,:,:,i+1).*Model.param.Model_scale(i+1);
                Model.param.Model_scale(i+1)=Model.param.Model_scale(i+1).*sc;
                x(:,:,:,:,:,:,i+1)=x(:,:,:,:,:,:,i+1)./Model.param.Model_scale(i+1);
                Model=update_constraints(Model,sc,i+1);
                
                p=x(:,:,:,:,:,:,i+1);
                p(p==Model.constraints(i+1,1))=1./Model.param.Model_scale(2).*exp(-Model.param.TR./(1000));
                p(p==Model.constraints(i+1,2))=1./Model.param.Model_scale(2).*exp(-Model.param.TR./(1000));
                x(:,:,:,:,:,:,i+1)=p;
            end
            
            Model=update_DA(Model,x);
            
            
            
        end
        
        function [x] =apply_constraints(Model,x)
            
            
                for(i=1:Model.param.nVar)
                
                pom=x(:,:,:,:,:,:,i);
                if(i==2)
                pom=real(pom);
                end
                
                pom(pom<Model.constraints(i,1))=Model.constraints(i,1);
                pom(pom>Model.constraints(i,2))=Model.constraints(i,2);
                x(:,:,:,:,:,:,i)=pom;
                end
            
        end
        
        
        
    end
    
end