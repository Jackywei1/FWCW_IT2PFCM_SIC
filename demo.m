%% BEST RUN WITH MATLAB R2018b!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interval type-2 possibilistic fuzzy clustering noisy image segmentation algorithm 
% with adaptive spatial constraints and local feature weighting & clustering weighting 
% International Journal of Approximate Reasoning
% This code was solely written by Tongyi Wei.
%
% Basically, you can run this code SEVERAL times to acquire the most desired result.
% It is welcomed to change the following parameters as you like to see what gonna happen.
%
% CUDA is required in this version. 
% However there is no need to install CUDA seperately since MATLAB has done all the work.
%
% Inputs:
% m - membership factor
% error - Minimum Error
% max_iter - Maximum iterations
% phi - Variance control parameters of Eq.(16)
% tao - Weighting factor
% sigm - Weighted regular parameter
% density - Mixed noise density
% cluster_num - Number of Clustering
% ============== Parameters of non-local spatial information================
% l - Side length of local block
% S - Side length of non-local block
% g - Attenuation of exponential function in Eqs. (10)-(11)
% sigma - Gaussian standard deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Intialization
clc
clear
warning off
close all
%% Input image
f_uint8 = imread('118035.jpg');
[row,col,d] = size (f_uint8); N = row * col;
figure(1),imshow(f_uint8),title('The original image');
%% Add noise
density = 0.05; f = (double(f_uint8))/255;
f = imnoise(f,'gaussian',0,density); f = imnoise(f,'salt & pepper',density); f = imnoise(f,'speckle',density);
f = f*255;
figure(2),imshow(uint8(f)),title('The noisy image');
%% Computing non-local spatial information
l = 7; s = 15; g=10; sigma = 4; 
f = gpuArray(f);
non_local_infomation = non_local_information(f, l, s, g, sigma);
figure(3),imshow(uint8(non_local_infomation)),title('The non-local spatial information');
%% Recombined Pixels
all_pixels=gather(reshape(double(f), N ,d));
all_pixels_xi=gather(reshape(double(non_local_infomation), N ,d));
%% Calculate alpha and beta values
difference =exp(20*(mean(abs( mean(all_pixels)-mean(all_pixels_xi)))).^2 + eps);
alpha = 1 ./ difference; beta = difference;
%% Parameter settings
k = 3 ;
m1 = 2; m2 = 4; theta1 = 3; theta2 = 5;
Cf=0.6; Cp=0.4; 
K=1; 
q=2; p=0.4; 
t_max=100;
beta_memory=0.3;
gamma_y=1./var(all_pixels); gamma_y(gamma_y==inf)=1;
gamma_xi=1./var(all_pixels_xi); gamma_xi(gamma_xi==inf)=1;
%% Initialize η
[ETA,C] = Initialization_ETA (all_pixels,all_pixels_xi,gamma_y,gamma_xi,alpha,beta,mean(m1,m2),k,K);
%% Start clustering
[C1,C2,UT_upper,UT_low, E1,E2]=FWCW_IT2PFCM_SIC(k,all_pixels,all_pixels_xi,C,gamma_y,gamma_xi,alpha,beta,q,p,m1,m2, ...
    theta1,theta2,Cp,Cf,ETA,t_max,beta_memory);
%% Assign labels
UT =(UT_upper+UT_low)/2;
best_clustering=zeros(1,N);
[~,Cluster]=max(UT,[],1);
if best_clustering ~= 0
    if accurcy_ave(repeat) > accurcy_ave(repeat-1)
        best_clustering = Cluster;
    end
else
    best_clustering = Cluster;
end
%% Show output
load Color_Map
cmap = reshape(best_clustering', [size(double(f_uint8),1) size(double(f_uint8),2)]);
FCM_result = label2rgb(cmap, Color_Map);
figure(4), imshow(FCM_result);title('result image');