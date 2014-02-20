function [ output_args ] = demo2_3d( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%============================================================
%               demo2 - denoise an image
% this is a run_file the demonstrate how to denoise an image, 
% using dictionaries. The methods implemented here are the same
% one as described in "Image Denoising Via Sparse and Redundant
% representations over Learned Dictionaries", (appeared in the 
% IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
%============================================================

clear

addpath('..')

bb=8; % block size
RR=4; % 冗余因素
K=RR*bb^2; % 字典中的原子数

load 'C:\data\predemo2_3d.mat'

%[IMinori,pp]=imread('images/barbara.bmp');
%[IMinnoise,pp]=imread('images/barbara_25_20.1658.bmp');
%imageName = 'images/barbara_25_20.1658.bmp';
%sigma = 25; 
%PSNRIn = 20.1658;

disp([datestr(now,0) ' sigma=' num2str(sigma)])
disp([datestr(now,0) ' PSNRIn=' num2str(PSNRIn)])

[IoutDCT,output] = denoiseImageDCT(IMinnoise, sigma, K);
PSNROut = 20*log10(255/sqrt(mean((IoutDCT(:)-IMinori(:)).^2)));
disp([datestr(now,0) ' DCT PSNROut=' num2str(PSNROut)])

[IoutGlobal,output] = denoiseImageGlobal(IMinnoise, sigma);
PSNROut = 20*log10(255/sqrt(mean((IoutGlobal(:)-IMinori(:)).^2)));
disp([datestr(now,0) ' Global PSNROut=' num2str(PSNROut)])

[IoutAdaptive,output] = denoiseImageKSVD(IMinnoise, sigma,K);
PSNROut = 20*log10(255/sqrt(mean((IoutAdaptive(:)-IMinori(:)).^2)));
disp([datestr(now,0) ' Adaptive PSNROut=' num2str(PSNROut)])


end

