function [IOut,output] = denoiseImageDCT(Image,sigma,K,varargin)
%==========================================================================
%  过完备的DCT来降噪
%==========================================================================
% function IOut = denoiseImageDCT(Image,sigma,bb,K)
% 用过完备的DCT稀疏表示每一个block，并计算 represented parts的平均值.

% INPUT ARGUMENTS : Image - 噪声灰度图
%                   sigma -噪声等级
%                   K -字典的原子信号的个数
%    Optional argumeters:              
%                  'blockSize' - 块的大小，默认是8.
%                  'errorFactor' - 乘以sigma得到允许的表示错误，默认1.15
%                  'maxBlocksToConsider' - 默认最大考虑250000块
%                  'slidingFactor' - 滑动距离，默认1.
%                  'waitBarOn' - 是否显示进度条
% OUTPUT ARGUMENTS : IOut - 和输入图片同大小的干净图片
%                    output ：
%                       D - 用于降噪的字典
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


% DCT=zeros(bb,Pn); %8行16列
% for k=0:1:Pn-1,
%     V=cos([0:1:bb-1]'*k*pi/Pn);
%     if k>0, V=V-mean(V); end;
%     DCT(:,k+1)=V/norm(V);
% end;
% DCT=kron(DCT,DCT); %两个8x16的结果是64x256



while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)
    slidingDis = slidingDis+1;
end

addpath('3d');
[blocks,idx] = my_im2col3d(Image,[bb,bb,bb],slidingDis);%结果是64x255025和1x255025


%字典D 512行 2048列
%DCT = blocks(:,1:2048);
load('C:\data\DCT3.mat')

% %imwrite(DCT,'myDCT.bmp')
% %imshow(DCT,[]); title('myDCT');
% %64x256  16 256/16 8 8 
% %I = displayDictionaryElementsAsImage(DCT, floor(sqrt(K)), floor(size(DCT,2)/floor(sqrt(K))),bb,bb,0);
% %title('The DCT dictionary');
% 
% [blocks,idx] = my_im2col(Image,[bb,bb],slidingDis);%结果是64x255025和1x255025
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
for jj = 1:10000:size(blocks,2) %从1到255025 每次跳10000个
    if (waitBarOn)
        waitbar(jj/counterForWaitBar);
    end
    jumpSize = min(jj+10000-1,                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       size(blocks,2)); %本次的最后一个，不能超出255025
    if (Reduce_DC)%1-10000列数据
        vecOfMeans = mean(blocks(:,jj:jumpSize)); %1-10000列，10000个列的平均值
        blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1); %扩展成64行10000列，并对对应列减去对应列的均值
    end
    %Coefs = mexOMPerrIterative(blocks(:,jj:jumpSize),DCT,errT);
    Coefs = OMPerr(DCT,blocks(:,jj:jumpSize),errT);
    if (Reduce_DC)
        blocks(:,jj:jumpSize)= DCT*Coefs + ones(size(blocks,1),1) * vecOfMeans; %将10000列系数xDCT 加上对应列的均值后回填到这10000列
    else
        blocks(:,jj:jumpSize)= DCT*Coefs ;
    end
end

count = 1;
Weight= zeros(NN1,NN2,NN3); %图片是512x512的    
IMout = zeros(NN1,NN2,NN3); %图片是512x512的
[rows,cols,zzs] = ind2sub(size(Image)-bb+1,idx); %255025个随机索引对应在图片中块的位置 
for i  = 1:length(cols)
    col = cols(i); row = rows(i); zz = zzs(i);
    block =reshape(blocks(:,count),[bb,bb,bb]);%将第count列还原为8x8x8大小
    IMout(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)=IMout(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)+block; %将block填到对应位置，原值加上
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