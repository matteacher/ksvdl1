function [ output_args ] = predemo2_3d( input_args )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
addpath('..')
load 'C:\data\image_0.mat'
testdata_3d = imageM; % (101:120,101:120,101:120);
I= testdata_3d;
%Ps=sum(sum(sum((I-mean(mean(mean(I)))).^2)));%signal power
Ps = sum((I(:)-mean(I(:))).^2);
SNR=0.1; %Ps/Pn = SNR
Pn = Ps/SNR;
sigma = sqrt(Pn);
sigma = 20;

bb=8; % block size
RR=4; % 冗余因素
K=RR*bb^2; % 字典中的原子数
K=2048;

IMinnoise = testdata_3d;
%[IoutDCT,output] = denoiseImageDCT(IMinnoise, sigma, K);
[IoutAdaptive,output] = denoiseImageKSVD(IMinnoise, sigma,K);
%PSNROut = 20*log10(255/sqrt(mean((IoutDCT(:)-IMinori(:)).^2)));
%disp([datestr(now,0) ' DCT PSNROut=' num2str(PSNROut)])


end

