function [images] = prep_images(param,rawdata)

    images=zeros(param.size_Data_cart);
    
    for(i=1:size(rawdata,5))
        
        FT = gpuNUFFT(param.k(:,:,i),col(param.w),param.osf,param.wg,param.sw,[param.nFE,param.nFE,param.nSlice],[],true);
        
        
        inft = @(x) FT'*reshape(x,[numel(rawdata)/param.nCh/size(rawdata,5),param.nCh]);
        
        
        images(:,:,:,:,:,i)=inft(rawdata(:,:,:,:,i));
        
    end


end


