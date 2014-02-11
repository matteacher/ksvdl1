%============================================================
%               demo2 - denoise an image
% this is a run_file the demonstrate how to denoise an image, 
% using dictionaries. The methods implemented here are the same
% one as described in "Image Denoising Via Sparse and Redundant
% representations over Learned Dictionaries", (appeared in the 
% IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
%============================================================

clear
bb=8; % block size
RR=4; % 冗余因素
K=RR*bb^2; % 字典中的原子数

load 'C:\data\predemo2.mat'

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


% %==========================================================================
% %   P E R F O R M   D E N O I S I N G   U S I N G   O V E R C O M P L E T E 
% %                        D C T    D I C T I O N A R Y
% %==========================================================================
% [IoutDCT,output] = denoiseImageDCT(IMin, sigma, K);
% %imwrite(IoutDCT,'IoutDCT.bmp')
% PSNROut = 20*log10(255/sqrt(mean((IoutDCT(:)-IMin0(:)).^2)));
% figure;
% subplot(1,3,1); imshow(IMin0,[]); title('Original clean image');
% subplot(1,3,2); imshow(IMin,[]); title(strcat(['Noisy image, ',num2str(PSNRIn),'dB']));
% subplot(1,3,3); imshow(IoutDCT,[]); title(strcat(['Clean Image by DCT dictionary, ',num2str(PSNROut),'dB']));
% figure;
% I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb,0,'KSVD');
% title('The DCT dictionary');
% 


% %==========================================================================
% %   P E R F O R M   D E N O I S I N G   U S I N G   G L O B A L 
% %           ( O R   G I V E N )   D I C T I O N A R Y
% %==========================================================================
% [IoutGlobal,output] = denoiseImageGlobal(IMin, sigma);
% 
% PSNROut = 20*log10(255/sqrt(mean((IoutGlobal(:)-IMin0(:)).^2)));
% figure;
% subplot(1,3,1); imshow(IMin0,[]); title('Original clean image');
% subplot(1,3,2); imshow(IMin,[]); title(strcat(['Noisy image, ',num2str(PSNRIn),'dB']));
% subplot(1,3,3); imshow(IoutGlobal,[]); title(strcat(['Clean Image by Global Trained dictionary, ',num2str(PSNROut),'dB']));
% figure;
% I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb);
% title('The dictionary trained on patches from natural images');
% 
% 



%==========================================================================
%   P E R F O R M   D E N O I S I N G   U S I N G   A   D I C T  I O N A R Y
%                  T R A I N E D   O N   N O I S Y   I M A G E
%==========================================================================
% [IoutAdaptive,output] = denoiseImageKSVD(IMin, sigma,K);
% 
% PSNROut = 20*log10(255/sqrt(mean((IoutAdaptive(:)-IMin0(:)).^2)));
% figure;
% subplot(1,3,1); imshow(IMin0,[]); title('Original clean image');
% subplot(1,3,2); imshow(IMin,[]); title(strcat(['Noisy image, ',num2str(PSNRIn),'dB']));
% subplot(1,3,3); imshow(IoutAdaptive,[]); title(strcat(['Clean Image by Adaptive dictionary, ',num2str(PSNROut),'dB']));
% 
% figure;
% I = displayDictionaryElementsAsImage(output.D, floor(sqrt(K)), floor(size(output.D,2)/floor(sqrt(K))),bb,bb);
% title('The dictionary trained on patches from the noisy image');