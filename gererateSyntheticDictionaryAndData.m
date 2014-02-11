function [Dictionary, data, coefs] = gererateSyntheticDictionaryAndData(N, L, dim, K, SNRdB)
%(信号总数1500，每个数据信号3个元素，每个数据20维，每个字典50个元素，噪声等级20)
 
randn('state',sum(100*clock)); %matlab5时代用来控制随机数发生器，推荐用rng('shuffle') 代替。clock是当前时间
rand('state',sum(100*clock));



Dictionary = randn(dim,K); %20行50列，50列20维的数据
Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));%字典的归一化 %默认是按列求和。50个和开方 分之1，转化为50*50的对角矩阵。20*50 * 50*50 =20*50

%20*1500 = 20*50  50*1500
%dim*N   = dim*K   K*N 
% data   =   D  *  X
%初始X的每一列有L个非零值
[data,coefs] = CreateDataFromDictionarySimple(Dictionary, N, L);
%输入参数是20x50的dictionary和1500和3
%输出的xOrig为50行1500列，每一列有3个随机数和47个0
%输出的 D = dictionary*xOrig; 结果是1500列20维的数据

if (SNRdB==0) | (SNRdB == 80) %?
    return
else
    noise = randn(size(data)); %与data同大小的随机噪声
    actualNoise = calcNoiseFromSNR(SNRdB,data, noise); %(20,1500列20维，随机的1500列20维)
    SNR = calcSNR(data, data+actualNoise);
    data =  data + actualNoise*SNR/SNRdB;   %根据随机的噪声计算需要的噪声
end

% 从字典创建1500条数据信号，每个数据信号从字典中取出3列
%(字典，信号总数1500，每个数据信号3个元素)
function [D,xOrig] = CreateDataFromDictionarySimple(dictionary, numElements, numCoef)
maxRangeOfCoef = 1;
resolution = 0.0001;

xOrig = zeros(size(dictionary,2),numElements); %50行1500列的全0
%vecOfValues = -1*maxRangeOfCoef:resolution:maxRangeOfCoef;
%coefs = randsrc(numCoef,numElements,vecOfValues);
coefs = randn(numCoef,numElements)*maxRangeOfCoef; %3行1500列的随机数
xOrig(1:numCoef,:) = coefs; %这3行1500列随机数填到xOrig的前3行，后47行全0

for i=1:size(xOrig,2) %从1到1500列
    xOrig(:,i) = xOrig(randperm(size(xOrig,1)),i); %第i列  将1-50打乱顺序，将第i列的50个元素打乱顺序，就是3个随机值分在了50行中
end
%dictionaryElementIndices = randsrc(numCoef*numElements,1,[1:size(dictionary,2)])   ; 
%matrixOfIndices = repmat([1:numElements],numCoef,1);
%xOrig(sub2ind(size(xOrig),dictionaryElementIndices,matrixOfIndices(:))) = coefs;

%输入参数是20x50的dictionary和1500和30，上面的代码生成了Orig，50行1500列，每一列有3个随机数和47个0
D = dictionary*xOrig; %D为1500列20维的数据


function  actualNoise = calcNoiseFromSNR(TargerSNR, signal, randomNoise) %(20,1500列20维，随机的1500列20维)
signal = signal(:); %从二维的20*1500转化为一维的30000个元素
randomNoiseRow = randomNoise(:);%也转成一维数组
signal_2 = sum(signal.^2); %所有元素的平方和
ActualNoise_2 = signal_2/(10^(TargerSNR/10)); %等级为20对应的实际噪声大小
noise_2 = sum(randomNoiseRow.^2); %所有元素的平方和
ratio = ActualNoise_2./noise_2; %信噪比？
%disp(randomNoiseRow.*repmat(sqrt(ratio),size(randomNoiseRow,1),1) - randomNoiseRow*sqrt(ratio))  %这样计算更简单
actualNoise = randomNoiseRow.*repmat(sqrt(ratio),size(randomNoiseRow,1),1); %将根号信噪比填充为为30000行1列，实际就是给一维噪声数组乘以信噪比
actualNoise = reshape(actualNoise,size(randomNoise));%变回20*1500的大小

function SNR = calcSNR(origSignal, noisySignal)
errorSignal = origSignal-noisySignal;
signal_2 = sum(origSignal.^2);
noise_2 = sum(errorSignal.^2);

SNRValues = 10*log10(signal_2./noise_2);
SNR = mean(SNRValues);
