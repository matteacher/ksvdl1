function [IOut,output] = denoiseImageKSVD(Image,sigma,K,varargin)
%用 噪声图片训练出的字典 来降噪
%==========================================================================
%训练出一个字典，然后稀疏表示每一块 and averaging the represented parts.
% 细节在"Image Denoising Via Sparse and Redundant representations over Learned Dictionaries", (IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
% 耗时过程:
%  1. KSVD迭代次数 
%  2. 参与训练的块的最大数量
% INPUT ARGUMENTS : Image - 灰度图
%                   sigma - 噪声等级
%                   K - 字典的原子数
%    Optional arguments:              
%                  'blockSize' - 方块的大小，默认为8.
%                  'errorFactor' - 默认1.15，乘以sigma得到允许的表示错误
%                  'maxBlocksToConsider' - 参与训练的块的最大数量（可以随机选出来）. 默认250000.
%                  'slidingFactor' - 块间滑动距离，默认1.图片太大时，会出于内存考虑增大该值
%                  'numKSVDIters' - KSVD迭代次数 默认是10或5（根据sigma调整）.小一点效果也不错
%                  'maxNumBlocksToTrainOn' - 参与训练的块的最大数量? 默认65000，对于大图可能不够
%                  'displayFlag' - 每次迭代后的显示，平均需要的系数数量
%                  'waitBarOn' - 显示等待进度条
% OUTPUT ARGUMENTS : Iout - 降噪结果
%                    output.D 训练出的字典
% =========================================================================

% 训练字典
reduceDC = 1;
[NN1,NN2,NN3] = size(Image);
waitBarOn = 1;
if (sigma > 5)
    numIterOfKsvd = 10;
else
    numIterOfKsvd = 5;
end
C = 1.15;
maxBlocksToConsider = 260000;
slidingDis = 1;
bb = 8;
z=30;
zs =[];
maxNumBlocksToTrainOn = 65000;
displayFlag = 1;

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
    if (strcmp(varargin{argI}, 'numKSVDIters'))
        numIterOfKsvd = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'blockSize'))
        bb = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'maxNumBlocksToTrainOn'))
        maxNumBlocksToTrainOn = varargin{argI+1};
    end
    if (strcmp(varargin{argI}, 'displayFlag'))
        displayFlag = varargin{argI+1};
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

if (sigma <= 5)
    numIterOfKsvd = 5;
end

% 训练字典

if(prod([NN1,NN2,NN3]-bb+1)> maxNumBlocksToTrainOn)
    randPermutation =  randperm(prod([NN1,NN2,NN3]-bb+1));
    selectedBlocks = randPermutation(1:maxNumBlocksToTrainOn);%从乱序中去除前max个
    blkMatrix = zeros(bb^3,maxNumBlocksToTrainOn); %64行max列的0
    for i = 1:maxNumBlocksToTrainOn
        [row,col,zz] = ind2sub(size(Image)-bb+1,selectedBlocks(i));%一维索引转成多维下标
        currBlock = Image(row:row+bb-1,col:col+bb-1,zz:zz+bb-1);
        blkMatrix(:,i) = currBlock(:);
    end
else
    blkNum= prod([NN1,NN2,NN3]-bb+1);
    blkMatrix = zeros(bb^3,blkNum); %64行max列的0
    for i = 1:blkNum
        [row,col,zz] = ind2sub(size(Image)-bb+1,i);%一维索引转成多维下标
        currBlock = Image(row:row+bb-1,col:col+bb-1,zz:zz+bb-1);
        blkMatrix(:,i) = currBlock(:);
    end
    %blkMatrix = im2col(Image,[bb,bb,bb],'sliding');%分成8*8的子矩阵，每个矩阵作为一列。sliding表示可以重叠，尽可能多
end

param.K = K; %
param.numIteration = numIterOfKsvd ;

param.errorFlag = 1; % 停止条件是达到固定错误，而不是固定非零系数
param.errorGoal = sigma*C; %允许的错误
param.preserveDCAtom = 0; %不保护第一列

Pn=ceil(K^(1/3));%16
% DCT=zeros(bb,Pn);%8行16列
% for k=0:1:Pn-1,%计算DCI的每一列
%     V=cos([0:1:bb-1]'*k*pi/Pn);
%     if k>0, V=V-mean(V); end;
%     DCT(:,k+1)=V/norm(V);
% end;
% DCT=kron(DCT,DCT);%?Kronecker乘法,结果是64*256

%DCT=Image(:,1:2048);
load('C:\data\DCT3.mat')

param.initialDictionary = DCT;%DCT(:,1:param.K );%初始字典为DCT，64*256
param.InitializationMethod =  'GivenMatrix';

if (reduceDC)%每个元素都减去当前列的均值
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
end

if (waitBarOn)
    counterForWaitBar = param.numIteration+1;%迭代次数
    h = waitbar(0,'Denoising In Process ...');
    param.waitBarHandle = h;
    param.counterForWaitBar = counterForWaitBar;
end


param.displayProgress = displayFlag;
[Dictionary,output] = KSVD(blkMatrix,param);
output.D = Dictionary;
save('C:\data\KSVDDic.mat','Dictionary','output')
% output =1;
% load('C:\data\KSVDDic.mat');

if (displayFlag)
    disp('finished Trainning dictionary');
end


% denoise the image using the resulted dictionary
errT = sigma*C;
IMout=zeros(NN1,NN2,NN3);
Weight=zeros(NN1,NN2,NN3);
%blocks = im2col(Image,[NN1,NN2],[bb,bb],'sliding');
while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)
    slidingDis = slidingDis+1;
end
[blocks,idx] = my_im2col(Image,[bb,bb,bb],slidingDis);

if (waitBarOn)
    newCounterForWaitBar = (param.numIteration+1)*size(blocks,2);
end


% go with jumps of 30000
for jj = 1:30000:size(blocks,2)
    if (waitBarOn)
        waitbar(((param.numIteration*size(blocks,2))+jj)/newCounterForWaitBar);
    end
    jumpSize = min(jj+30000-1,size(blocks,2));
    if (reduceDC)
        vecOfMeans = mean(blocks(:,jj:jumpSize));
        blocks(:,jj:jumpSize) = blocks(:,jj:jumpSize) - repmat(vecOfMeans,size(blocks,1),1);
    end
    
    %Coefs = mexOMPerrIterative(blocks(:,jj:jumpSize),Dictionary,errT);
    Coefs = OMPerr(Dictionary,blocks(:,jj:jumpSize),errT);
    

%         myDic  = Dictionary;
%         myData = blocks(:,jj:jumpSize);
%         numOfAtoms = size(myDic,2);
%         numOfSignals = size(myData,2);
%         mydim=size(myDic,1);
%         Coefs = zeros([numOfAtoms,numOfSignals]);
%         for iii = 1:numOfSignals
%             [s, err_mse, iter_time]=greed_omp_qr(myData(:,iii),myDic,numOfAtoms);
%             Coefs(:,iii)=s';
%         end
        
        
        
    if (reduceDC)
        blocks(:,jj:jumpSize)= Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
    else
        blocks(:,jj:jumpSize)= Dictionary*Coefs ;
    end
end

count = 1;
Weight = zeros(NN1,NN2,NN3);
IMout = zeros(NN1,NN2,NN3);
[rows,cols,zzs] = ind2sub(size(Image)-bb+1,idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);zz=zzs(i);        
    block =reshape(blocks(:,count),[bb,bb,bb]);
    IMout(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)=IMout(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)+block;
    Weight(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)=Weight(row:row+bb-1,col:col+bb-1,zz:zz+bb-1)+ones([bb,bb,bb]);
    count = count+1;
end;

if (waitBarOn)
    close(h);
end

if size(zs,2)~= 0
    IOut = zeros([size(IMout) size(zs,2)]);
    for  i=1:size(zs,2)
    IOut(:,:,i)=(Image+1/zs(1,i)*sigma*IMout)./(1+1/zs(1,i)*sigma*Weight);
    end
else
    IOut = (Image+1/z*sigma*IMout)./(1+1/z*sigma*Weight);
end

%output.Coefs =Coefs;
