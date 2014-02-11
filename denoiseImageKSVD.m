function [IOut,output] = denoiseImageKSVD(Image,sigma,K,varargin)
%�� ����ͼƬѵ�������ֵ� ������
%==========================================================================
%ѵ����һ���ֵ䣬Ȼ��ϡ���ʾÿһ�� and averaging the represented parts.
% ϸ����"Image Denoising Via Sparse and Redundant representations over Learned Dictionaries", (IEEE Trans. on Image Processing, Vol. 15, no. 12, December 2006).
% ��ʱ����:
%  1. KSVD�������� 
%  2. ����ѵ���Ŀ���������
% INPUT ARGUMENTS : Image - �Ҷ�ͼ
%                   sigma - �����ȼ�
%                   K - �ֵ��ԭ����
%    Optional arguments:              
%                  'blockSize' - ����Ĵ�С��Ĭ��Ϊ8.
%                  'errorFactor' - Ĭ��1.15������sigma�õ������ı�ʾ����
%                  'maxBlocksToConsider' - ����ѵ���Ŀ������������������ѡ������. Ĭ��250000.
%                  'slidingFactor' - ��们�����룬Ĭ��1.ͼƬ̫��ʱ��������ڴ濼�������ֵ
%                  'numKSVDIters' - KSVD�������� Ĭ����10��5������sigma������.Сһ��Ч��Ҳ����
%                  'maxNumBlocksToTrainOn' - ����ѵ���Ŀ���������? Ĭ��65000�����ڴ�ͼ���ܲ���
%                  'displayFlag' - ÿ�ε��������ʾ��ƽ����Ҫ��ϵ������
%                  'waitBarOn' - ��ʾ�ȴ�������
% OUTPUT ARGUMENTS : Iout - ������
%                    output.D ѵ�������ֵ�
% =========================================================================

% ѵ���ֵ�
reduceDC = 1;
[NN1,NN2] = size(Image);
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

% ѵ���ֵ�

if(prod([NN1,NN2]-bb+1)> maxNumBlocksToTrainOn)
    randPermutation =  randperm(prod([NN1,NN2]-bb+1));
    selectedBlocks = randPermutation(1:maxNumBlocksToTrainOn);%��������ȥ��ǰmax��
    blkMatrix = zeros(bb^2,maxNumBlocksToTrainOn); %64��max�е�0
    for i = 1:maxNumBlocksToTrainOn
        [row,col] = ind2sub(size(Image)-bb+1,selectedBlocks(i));%һά����ת�ɶ�ά�±�
        currBlock = Image(row:row+bb-1,col:col+bb-1);
        blkMatrix(:,i) = currBlock(:);
    end
else
    blkMatrix = im2col(Image,[bb,bb],'sliding');%�ֳ�8*8���Ӿ���ÿ��������Ϊһ�С�sliding��ʾ�����ص��������ܶ�
end

param.K = K; %K=256
param.numIteration = numIterOfKsvd ;

param.errorFlag = 1; % ֹͣ�����Ǵﵽ�̶����󣬶����ǹ̶�����ϵ��
param.errorGoal = sigma*C; %�����Ĵ���
param.preserveDCAtom = 0; %��������һ��

Pn=ceil(sqrt(K));%16
DCT=zeros(bb,Pn);%8��16��
for k=0:1:Pn-1,%����DCI��ÿһ��
    V=cos([0:1:bb-1]'*k*pi/Pn);
    if k>0, V=V-mean(V); end;
    DCT(:,k+1)=V/norm(V);
end;
DCT=kron(DCT,DCT);%?Kronecker�˷�,�����64*256

param.initialDictionary = DCT(:,1:param.K );%��ʼ�ֵ�ΪDCT��64*256
param.InitializationMethod =  'GivenMatrix';

if (reduceDC)%ÿ��Ԫ�ض���ȥ��ǰ�еľ�ֵ
    vecOfMeans = mean(blkMatrix);
    blkMatrix = blkMatrix-ones(size(blkMatrix,1),1)*vecOfMeans;
end

if (waitBarOn)
    counterForWaitBar = param.numIteration+1;%��������
    h = waitbar(0,'Denoising In Process ...');
    param.waitBarHandle = h;
    param.counterForWaitBar = counterForWaitBar;
end


param.displayProgress = displayFlag;
[Dictionary,output] = KSVD(blkMatrix,param);
output.D = Dictionary;

if (displayFlag)
    disp('finished Trainning dictionary');
end


% denoise the image using the resulted dictionary
errT = sigma*C;
IMout=zeros(NN1,NN2);
Weight=zeros(NN1,NN2);
%blocks = im2col(Image,[NN1,NN2],[bb,bb],'sliding');
while (prod(floor((size(Image)-bb)/slidingDis)+1)>maxBlocksToConsider)
    slidingDis = slidingDis+1;
end
[blocks,idx] = my_im2col(Image,[bb,bb],slidingDis);

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
    if (reduceDC)
        blocks(:,jj:jumpSize)= Dictionary*Coefs + ones(size(blocks,1),1) * vecOfMeans;
    else
        blocks(:,jj:jumpSize)= Dictionary*Coefs ;
    end
end

count = 1;
Weight = zeros(NN1,NN2);
IMout = zeros(NN1,NN2);
[rows,cols] = ind2sub(size(Image)-bb+1,idx);
for i  = 1:length(cols)
    col = cols(i); row = rows(i);        
    block =reshape(blocks(:,count),[bb,bb]);
    IMout(row:row+bb-1,col:col+bb-1)=IMout(row:row+bb-1,col:col+bb-1)+block;
    Weight(row:row+bb-1,col:col+bb-1)=Weight(row:row+bb-1,col:col+bb-1)+ones(bb);
    count = count+1;
end;

if (waitBarOn)
    close(h);
end
disp(['start for different sigam*lambda ' datestr(now,0)])
if size(zs,2)~= 0
    IOut = zeros([size(IMout) size(zs,2)]);
    for  i=1:size(zs,2)
    IOut(:,:,i)=(Image+1/zs(1,i)*sigma*IMout)./(1+1/zs(1,i)*sigma*Weight);
    end
else
    IOut = (Image+1/z*sigma*IMout)./(1+1/z*sigma*Weight);
end

%output.Coefs =Coefs;