classdef operator
    %OPERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    
    
    properties
        param
        rez
    end
    
    methods
        
        function [operator] =operator(param)
            
            operator.param=param;
            
            
            
        end
        
        
        function [out] = forward_forward(operator,model,x)
            
            data=repmat(forward(model,x),1,1,1,1,operator.param.nCh,1,1,1);
            out=zeros(operator.param.size_kSpace);
            for(i=1:size(data,6))
                
                FT = gpuNUFFT(operator.param.k(:,:,i),col(operator.param.w),operator.param.osf,operator.param.wg,operator.param.sw,[operator.param.nFE,operator.param.nFE,operator.param.nSlice],operator.param.coil,true);
                
                
                nft  = @(x) reshape(FT*x,[model.param.size_kSpace_gradient(1:3),model.param.size_kSpace_gradient(4)]);
                %inft = @(x) FT'*reshape(x,[nFE*nSpokes,nCh]);
                
                for(j=1:size(data,4))
                    
                    out(:,:,j,:,i)=nft((squeeze(data(:,:,:,j,:,i))));
                    
                end
                
                
                
                
            end
            
            out(isnan(out))=1e-20;
            out(isinf(out))=1e-20;
            
            
            
        end
        
        
        
        function [out] = forward_gradient(operator,model,x)
            
            
            data=repmat((model.DA).*repmat(x,[1,1,1,1,1,operator.param.nPar,1]),[1,1,1,1,operator.param.nCh,1,1,1,1]);
            
            out=(zeros(operator.param.size_kSpace_gradient));
            
            for(i=1:size(model.DA,6))
                
                FT = gpuNUFFT(operator.param.k(:,:,i),col(operator.param.w),operator.param.osf,operator.param.wg,operator.param.sw,[operator.param.nFE,operator.param.nFE,operator.param.nSlice],operator.param.coil,true);
                
                
                nft  = @(x) reshape(FT*x,operator.param.nFE,operator.param.nSpokes,operator.param.nCh);
                %inft = @(x) FT'*reshape(x,[nFE*nSpokes,nCh]);
                
                for(j=1:size(data,4))
                    
                    out(:,:,j,:,i,1)=nft((squeeze(sum(data(:,:,:,j,:,i,:),7))));
                    
                    
                end
                
                
                
            end
            
            out=sum(out,6);
            out(isnan(out))=1e-20;
            out(isinf(out))=1e-20;
            
            
        end
        
        function [out] = adjoint_gradient(operator,model,data)
            
            
            
            
            out=(zeros(operator.param.size_Model_data_grad));
            
            for(i=1:size(model.DA,6))
                
                FT = gpuNUFFT(operator.param.k(:,:,i),col(operator.param.w),operator.param.osf,operator.param.wg,operator.param.sw,[operator.param.nFE,operator.param.nFE,operator.param.nSlice],(operator.param.coil),true);
                
                
                % nft  = @(x) reshape(FT*x,operator.param.nFE,operator.param.nSpokes,operator.param.nCh);
                inft = @(x) FT'*reshape(x,[operator.param.nFE*operator.param.nSpokes,operator.param.nCh]);
                
                for(j=1:size(data,3))
                    
                    out(:,:,:,j,:,i,1)=inft((squeeze(data(:,:,j,:,i))));
                    
                    
                end
                
                
                
            end
            
            
            out(:,:,:,:,:,:,2)=out(:,:,:,:,:,:,1);
            
            out=sum(conj(model.DA).*out,6);
            out(isnan(out))=1e-20;
            out(isinf(out))=1e-20;
            
            
        end
        
        function [operator] = calc_rez(operator,Model,uk,rawdata)
            
            operator.rez=(rawdata-forward_forward(operator,Model,uk)+ forward_gradient(operator,Model,uk));
            
            
        end
        
        
        
    end
    
    
end

