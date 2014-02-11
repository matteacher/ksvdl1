sigma = 25; 
pathForImages ='';
imageName = 'images/barbara.png';
[IMinori,pp]=imread(strcat([pathForImages,imageName]));
IMinori=im2double(IMinori);
if (length(size(IMinori))>2)
    IMinori = rgb2gray(IMin0); %��ɫͼ��ת�ɻҶ�ͼ������
end
if (max(IMinori(:))<2) 
    IMinori = IMinori*255; %��ֵͼת�ɻҶ�ͼ
end

%imwrite(uint8(IMinori),'images/barbara.bmp')
IMinnoise=IMinori+sigma*randn(size(IMinori)); %���������Ϊ����ͼƬ
PSNRIn = 20*log10(255/sqrt(mean((IMinnoise(:)-IMinori(:)).^2)));  %����ͼƬ�������ȼ�

disp([datestr(now,0) ' sigma=' num2str(sigma)])
disp([datestr(now,0) ' PSNRIn=' num2str(PSNRIn)])
%imwrite(uint8(IMinnoise),['images/barbara_' num2str(sigma) '_' num2str(PSNRIn) '.bmp'])

mkdir('C:\data')
save 'C:\data\predemo2.mat' sigma PSNRIn IMinori IMinnoise 