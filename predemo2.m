sigma = 25; 
pathForImages ='';
imageName = 'images/barbara.png';
[IMinori,pp]=imread(strcat([pathForImages,imageName]));
IMinori=im2double(IMinori);
if (length(size(IMinori))>2)
    IMinori = rgb2gray(IMin0); %彩色图像转成灰度图来处理
end
if (max(IMinori(:))<2) 
    IMinori = IMinori*255; %二值图转成灰度图
end

%imwrite(uint8(IMinori),'images/barbara.bmp')
IMinnoise=IMinori+sigma*randn(size(IMinori)); %添加噪声作为输入图片
PSNRIn = 20*log10(255/sqrt(mean((IMinnoise(:)-IMinori(:)).^2)));  %输入图片的噪声等级

disp([datestr(now,0) ' sigma=' num2str(sigma)])
disp([datestr(now,0) ' PSNRIn=' num2str(PSNRIn)])
%imwrite(uint8(IMinnoise),['images/barbara_' num2str(sigma) '_' num2str(PSNRIn) '.bmp'])

mkdir('C:\data')
save 'C:\data\predemo2.mat' sigma PSNRIn IMinori IMinnoise 