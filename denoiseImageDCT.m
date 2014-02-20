function [IOut,output] = denoiseImageDCT(Image,sigma,K,varargin)
%==========================================================================
%  ���걸��DCT������
%==========================================================================
% function IOut = denoiseImageDCT(Image,sigma,bb,K)
% �ù��걸��DCTϡ���ʾÿһ��block�������� represented parts��ƽ��ֵ.

% INPUT ARGUMENTS : Image - �����Ҷ�ͼ
%                   sigma -�����ȼ�
%                   K -�ֵ��ԭ���źŵĸ���
%    Optional argumeters:              
%                  'blockSize' - ��Ĵ�С��Ĭ����8.
%                  'errorFactor' - ����sigma�õ�����ı�ʾ����Ĭ��1.15
%                  'maxBlocksToConsider' - Ĭ�������250000��
%                  'slidingFactor' - �������룬Ĭ��1.
%                  'waitBarOn' - �Ƿ���ʾ������
% OUTPUT ARGUMENTS : IOut - ������ͼƬͬ��С�ĸɾ�ͼƬ
%                    output ��
%                       D - ���ڽ�����ֵ�
% =========================================================================


Reduce_DC = 1;
[NN1,NN2,NN3] = size(Image);
C = 1.15;
waitBarOn = 1;
maxBlocksToConsider = 1500000;%1500000;
slidingDis = 1;
bb = 8;
z=30;
zs =[];
if(length(varargin)==1)
    varargin=varargin{1};
end
for argI = 1:2:length(varargin)
    if (strcmp(varargin{argI}, 'slidingFactor'))
        slidingDis = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'errorFactor'))
        C = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'maxBlocksToConsider'))
        maxBlocksToConsider = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'blockSize'))
        bb = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'waitBarOn'))
        waitBarOn = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'z'))
        z = varargin{argI+1};
    end
    
    if (strcmp(varargin{argI}, 'zs'))
        zs = varargin{argI+1};
    end
end
errT = C*sigma;
% Create an initial dictionary from the DCT frame
Pn=ceil(K^(1/3)); %K=256  Pn=16


% DCT=zeros(bb,Pn); %8��16��
% for k=0:1:Pn-1,
%     V=cos([0:1:bb-1]'*k*pi/Pn);
%     if k>0, V=V-mean(V); end;
%     DCT(:,k+1)=V/norm(V);
% end;
% DCT=kron(DCT,DCT); %����8x16�Ľ����64x256



while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)
    slidingDis = slidingDis+1;
end

addpath('3d');
[blocks,idx] = my_im2col3d(Image,[bb,bb,bb],slidingDis);%�����64x255025��1x255025


%�ֵ�D 512�� 2048��
%DCT = blocks(:,1:2048);
load('C:\data\DCT3.mat')

% %imwrite(DCT,'myDCT.bmp')
% %imshow(DCT,[]); title('myDCT');
% %64x256  16 256/16 8 8 
% %I = displayDictionaryElementsAsImage(DCT, floor(sqrt(K)), floor(size(DCT,2)/floor(sqrt(K))),bb,bb,0);
% %title('The DCT dictionary');
% 
% [blocks,idx] = my_im2col(Image,[bb,bb],slidingDis);%�����64x255025��1x255025
% 
% %512+64*2=640
% save 'IM3' Image
% Image2 = zeros([640 640]);
% %Image2(:)=255;
% %ind=1:8:(63*8+1);
% totalx=zeros([256 64*64]);
% path = 'dir003/';
% for ii= 9:19
%     for jj = 43:53
%         %ii=6;jj =7;
%         
%         yy=Image((ii-1)*8+1:ii*8,(jj-1)*8+1:jj*8);
%         Image2((ii-1)*10+2:ii*10-1,(jj-1)*10+2:jj*10-1)=yy;
%         scale = 16;
%         imwrite(uint8(yy),[path 'real' num2str(ii) '_' num2str(jj) '.bmp'])
%         imwrite(uint8(imresize(yy,scale,'nearest')),[path num2str(ii) '_' num2str(jj) '.bmp'])
%         yyy=reshape(yy,[64,1]);
%         xx = OMPerr(DCT,yyy,errT);
%         totalx(:,(ii-1)*64+jj)=xx;
%         I = displayDictionaryElementsAsImage(DCT, floor(sqrt(K)), floor(size(DCT,2)/floor(sqrt(K))),bb,bb,0);
%         for i =1:256
%             if(xx(i) ~= 0)
%                 ihang = floor(i/16)+1;
%                 ilie = mod(i,16);
%                 if ilie == 0
%                    ilie =16 
%                 end
%                 dctij = I((ihang-1)*9+2:ihang*9,(ilie-1)*9+2:ilie*9,:);
%                 imwrite(dctij,[path 'real'  num2str(ii) '_' num2str(jj) '_' num2str(i) '_' num2str(xx(i))  '.bmp'])
%                imwrite(imresize(dctij,scale,'nearest'),[path num2str(ii) '_' num2str(jj) '_' num2str(i) '_' num2str(xx(i)) '.bmp'])
%                %imwrite(uint8(reshape(DCT(:,i),[8,8])),[path 'real'  num2str(ii) '_' num2str(jj) '_' num2str(i) '.bmp'])
%                %imwrite(uint8(imresize(reshape(DCT(:,i),[8,8]),scale)),[path num2str(ii) '_' num2str(jj) '_' num2str(i) '.bmp'])
%             end
%         end
%     end
% end
% imwrite(uint8(Image2),'Image2.bmp');
% %reshape(DCT(:,xx),[8,8])
% save 'totalx' totalx
if (waitBarOn)
    counterForWaitBar = size(blocks,2); %255025
    h = waitbar(0,'Denoising In Process ...');
end
% go with jumps of 10000
for jj = 1:10000:size(blocks,2) %��1��255025 ÿ����10000��
    if (waitBarOn)
        waitbar(jj/counterForWaitBar);
    end
    jumpSize = min(jj+10000-1,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       size(blocks,2)); %���ε����һ�������ܳ���255025
    if (Reduce_DC)%1-10000������
        vecOfMeans = mean(blocks(:,jj:jumpSize)); %1-10000�У�10000���е�ƽ��ֵ
        blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1); %��չ��64��10000�У����Զ�Ӧ�м�ȥ��Ӧ�еľ�ֵ
    end
    %Coefs = mexOMPerrIterative(blocks(:,jj:jumpSize),DCT,errT);
    Coefs = OMPerr(DCT,blocks(:,jj:jumpSize),errT);
    if (Reduce_DC)
        blocks(:,jj:jumpSize)= DCT*Coefs + ones(size(blocks,1),1) * vecOfMeans; %��10000��ϵ��xDCT ���϶�Ӧ�еľ�ֵ������10000��
    else
        blocks(:,jj:jumpSize)= DCT*Coefs ;
    end
end

count = 1;
Weight= zeros(NN1,NN2,NN3); %ͼƬ��512x512��    
IMout = zeros(NN1,NN2,NN3); %ͼƬ��512x512��
[rows,cols,zzs] = ind2sub(size(Image)-bb+1,idx); %255025�����������Ӧ��ͼƬ�п��λ�� 
for i  = 1:length(cols)
    col = cols(i); row = rows(i); zz = zzs(i);
    block =reshape(blocks(:,count),[bb,bb,bb]);%����count�л�ԭΪ8x8x8��С
    IMout(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)=IMout(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)+block; %��block���Ӧλ�ã�ԭֵ����
    Weight(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)=Weight(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)+ones([bb bb bb]);
    count = count+1;
end;
if (waitBarOn)
    close(h);
end
%disp(['start for different sigam*lambda ' datestr(now,0)])
%imshow(IMout,[]); title('imout');
if size(zs,2)~= 0
    IOut = zeros([size(IMout) size(zs,2)]);
    for  i=1:size(zs,2)
    IOut(:,:,i)=(Image+1/zs(1,i)*sigma*IMout)./(1+1/zs(1,i)*sigma*Weight);
    end
else
    IOut = (Image+1/z*sigma*IMout)./(1+1/z*sigma*Weight);
end


%imshow(IOut,[]); title('iout');
output.D = DCT;
output.Coefs =Coefs;