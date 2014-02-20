function [ output_args ] = predemo2_3d( input_args )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
addpath('..')
load 'C:\data\image_0.mat'
load 'C:\data\DCT3.mat' %8^3 12^3

IMinnoise= imageM(101:120,101:120,101:120);
%Ps=sum(sum(sum((I-mean(mean(mean(I)))).^2)));
% Ps = sum((IMinnoise(:)-mean(IMinnoise(:))).^2); %signal power
% SNR=0.1; %Ps/Pn = SNR
% Pn = Ps/SNR;
% sigma = sqrt(Pn);
sigma = 20;

bb=8; % block size
% RR=4; % 冗余因素
% K=RR*bb^2; % 字典中的原子数
K=1728; %12^3

%[IoutDCT,output] = denoiseImageDCT(IMinnoise, sigma, K);
%save(['C:\data\IoutDCT_sigma' num2str(sigma) '_blocks1500000' '.mat'],'IoutDCT','output')

[IoutAdaptive,output] = denoiseImageKSVD(IMinnoise, sigma,K);
save(['C:\data\IoutAdaptive_sigma' num2str(sigma) '_blocks65000' '.mat'],'IoutAdaptive','output')

%PSNROut = 20*log10(255/sqrt(mean((IoutDCT(:)-IMinori(:)).^2)));
%disp([datestr(now,0) ' DCT PSNROut=' num2str(PSNROut)])


end

