function [Dictionary, data, coefs] = gererateSyntheticDictionaryAndData(N, L, dim, K, SNRdB)
%(�ź�����1500��ÿ�������ź�3��Ԫ�أ�ÿ������20ά��ÿ���ֵ�50��Ԫ�أ������ȼ�20)
 
randn('state',sum(100*clock)); %matlab5ʱ������������������������Ƽ���rng('shuffle') ���档clock�ǵ�ǰʱ��
rand('state',sum(100*clock));



Dictionary = randn(dim,K); %20��50�У�50��20ά������
Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));%�ֵ�Ĺ�һ�� %Ĭ���ǰ�����͡�50���Ϳ��� ��֮1��ת��Ϊ50*50�ĶԽǾ���20*50 * 50*50 =20*50

%20*1500 = 20*50  50*1500
%dim*N   = dim*K   K*N 
% data   =   D  *  X
%��ʼX��ÿһ����L������ֵ
[data,coefs] = CreateDataFromDictionarySimple(Dictionary, N, L);
%���������20x50��dictionary��1500��3
%�����xOrigΪ50��1500�У�ÿһ����3���������47��0
%����� D = dictionary*xOrig; �����1500��20ά������

if (SNRdB==0) | (SNRdB == 80) %?
    return
else
    noise = randn(size(data)); %��dataͬ��С���������
    actualNoise = calcNoiseFromSNR(SNRdB,data, noise); %(20,1500��20ά�������1500��20ά)
    SNR = calcSNR(data, data+actualNoise);
    data =  data + actualNoise*SNR/SNRdB;   %�������������������Ҫ������
end

% ���ֵ䴴��1500�������źţ�ÿ�������źŴ��ֵ���ȡ��3��
%(�ֵ䣬�ź�����1500��ÿ�������ź�3��Ԫ��)
function [D,xOrig] = CreateDataFromDictionarySimple(dictionary, numElements, numCoef)
maxRangeOfCoef = 1;
resolution = 0.0001;

xOrig = zeros(size(dictionary,2),numElements); %50��1500�е�ȫ0
%vecOfValues = -1*maxRangeOfCoef:resolution:maxRangeOfCoef;
%coefs = randsrc(numCoef,numElements,vecOfValues);
coefs = randn(numCoef,numElements)*maxRangeOfCoef; %3��1500�е������
xOrig(1:numCoef,:) = coefs; %��3��1500��������xOrig��ǰ3�У���47��ȫ0

for i=1:size(xOrig,2) %��1��1500��
    xOrig(:,i) = xOrig(randperm(size(xOrig,1)),i); %��i��  ��1-50����˳�򣬽���i�е�50��Ԫ�ش���˳�򣬾���3�����ֵ������50����
end
%dictionaryElementIndices = randsrc(numCoef*numElements,1,[1:size(dictionary,2)])   ; 
%matrixOfIndices = repmat([1:numElements],numCoef,1);
%xOrig(sub2ind(size(xOrig),dictionaryElementIndices,matrixOfIndices(:))) = coefs;

%���������20x50��dictionary��1500��30������Ĵ���������Orig��50��1500�У�ÿһ����3���������47��0
D = dictionary*xOrig; %DΪ1500��20ά������


function  actualNoise = calcNoiseFromSNR(TargerSNR, signal, randomNoise) %(20,1500��20ά�������1500��20ά)
signal = signal(:); %�Ӷ�ά��20*1500ת��Ϊһά��30000��Ԫ��
randomNoiseRow = randomNoise(:);%Ҳת��һά����
signal_2 = sum(signal.^2); %����Ԫ�ص�ƽ����
ActualNoise_2 = signal_2/(10^(TargerSNR/10)); %�ȼ�Ϊ20��Ӧ��ʵ��������С
noise_2 = sum(randomNoiseRow.^2); %����Ԫ�ص�ƽ����
ratio = ActualNoise_2./noise_2; %����ȣ�
%disp(randomNoiseRow.*repmat(sqrt(ratio),size(randomNoiseRow,1),1) - randomNoiseRow*sqrt(ratio))  %�����������
actualNoise = randomNoiseRow.*repmat(sqrt(ratio),size(randomNoiseRow,1),1); %��������������ΪΪ30000��1�У�ʵ�ʾ��Ǹ�һά����������������
actualNoise = reshape(actualNoise,size(randomNoise));%���20*1500�Ĵ�С

function SNR = calcSNR(origSignal, noisySignal)
errorSignal = origSignal-noisySignal;
signal_2 = sum(origSignal.^2);
noise_2 = sum(errorSignal.^2);

SNRValues = 10*log10(signal_2./noise_2);
SNR = mean(SNRValues);
