function [param] = setup_param(rawdata,nSlice)

param.nFE=size(rawdata,1);
param.nCh=size(rawdata,4);
param.nPar=size(rawdata,5);
param.nt=size(rawdata,3);
param.nSlice=nSlice;
param.size_kSpace=size(rawdata);
param.nVar=2; % temp

param.nSpokes=size(rawdata,2);

param.size_kSpace_gradient=[size(rawdata),param.nVar];
param.size_Data_cart=[param.nFE,param.nFE,param.nSlice,param.nt,param.nCh,param.nPar];
param.size_Data_cart_grad=[param.size_Data_cart,param.nVar];
param.size_Model_data=param.size_Data_cart;
param.size_Model_data(5)=1;
param.size_Model_data(6)=1;
param.size_Model_data_grad=[param.size_Model_data,param.nVar];

param.osf= 2;
param.wg= 3;
param.sw= 8;

param.ratio=200;
m=param.nFE;
n=param.nFE;

param.dx     = @(x) [diff(x,1,2), x(:,1,:,:,:,:,:,:)-x(:,m,:,:,:,:,:,:)];       % discrete x-derivative
param.dy     = @(x) [diff(x,1,1); x(1,:,:,:,:,:,:,:)-x(n,:,:,:,:,:,:,:)];       % discrete y-derivative
param.dxt    = @(x) [x(:,m,:,:,:,:,:,:)-x(:,1,:,:,:,:,:,:),-diff(x,1,2)];      % transpose x-derivative
param.dyt    = @(x) [x(n,:,:,:,:,:,:,:)-x(1,:,:,:,:,:,:,:);-diff(x,1,1)];      % transpose y-derivative



param.grad=@(x) cat(8,param.dx(x),param.dy(x));


param.div1=@(x) (param.dxt(x(:,:,:,:,:,:,:,1))+param.dyt(x(:,:,:,:,:,:,:,2)));

param.E=@(x) cat(8,param.dx(x(:,:,:,:,:,:,:,1)),param.dy(x(:,:,:,:,:,:,:,2)),0.5*(param.dy(x(:,:,:,:,:,:,:,1))+param.dx(x(:,:,:,:,:,:,:,2))) );

param.div2=@(x) (cat(8,param.dxt(x(:,:,:,:,:,:,:,1))+param.dyt(x(:,:,:,:,:,:,:,3)),param.dyt(x(:,:,:,:,:,:,:,2))+param.dxt(x(:,:,:,:,:,:,:,3))));





end

