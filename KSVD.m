function [Dictionary,output] = KSVD(Data,param)
% DataΪn*N����N��nά�ź�
% �����param  
%    K, ...                    ѵ�������ֵ��ж�����
%    numIteration,...          ��������
%    errorFlag...             %��ʽ
%    (optional) L,...                 % 0��ʽ��OMP������ϵ��������ԭ�ӵĸ�����
%    (optional) errorGoal, ...        % 1��ʽ��ÿ���ź�����ı�����
%    preserveDCAtom...      Ϊ1ʱ���ֵ�ĵ�һ�б��ֲ��䣬��ʵ��ͼ���п��ܻ��õ� 

%    InitializationMethod,...  �ֵ�ĳ�ʼ�������� 'DataElements' (�������źų�ʼ��)�� 'GivenMatrix' (������initialDictionary��ʼ��).
%    (optional) initialDictionary,
%    (optional) TrueDictionary, ...        % ÿ�ζ����㲢��ʾ��ѵ�������ֵ�ʹ��ֵ�Ĳ�ͬ
%    displayProgress, ...      ����ģʽ�� 0��ʽ��ʾƽ����ʾ���RMSE��1��ʽ��ʾƽ��ϵ������
%                                
% =========================================================================
% OUTPUT ARGUMENTS:
%  Dictionary                  n��K�е��ֵ�
%  �����output                      
%    CoefMatrix                  ���յ�ϵ������x (it should hold that Data equals approximately Dictionary*output.CoefMatrix.
%    ratio                       TrueDictionaryʱ��holds��ΪnumIteration ����������ÿ�ε�����detection ratios 
%    totalerr                    ÿ�ε������ܵı�ʾ����0��ʽ������ģʽ
%    numCoef                    ��ΪnumIteration ������������ÿ�ε��� ƽ����Ҫϵ�������ź� 1��ʽ������ģʽ
% =========================================================================

if (~isfield(param,'displayProgress'))
    param.displayProgress = 0;
end
totalerr(1) = 99999;
if (isfield(param,'errorFlag')==0)
    param.errorFlag = 0;
end

if (isfield(param,'TrueDictionary'))
    displayErrorWithTrueDictionary = 1;
    ErrorBetweenDictionaries = zeros(param.numIteration+1,1); %��TrueDic�Ĵ���Ҳ��ʾ����
    ratio = zeros(param.numIteration+1,1);
else
    displayErrorWithTrueDictionary = 0;
	ratio = 0;
end
if (param.preserveDCAtom>0)
    FixedDictionaryElement(1:size(Data,1),1) = 1/sqrt(size(Data,1));
else
    FixedDictionaryElement = [];
end
% coefficient calculation method is OMP with fixed number of coefficients

if (size(Data,2) < param.K)
    disp('Size of data is smaller than the dictionary size. Trivial solution...');
    Dictionary = Data(:,1:size(Data,2));
    return;
elseif (strcmp(param.InitializationMethod,'DataElements')) % �ֵ�ĳ�ʼ�������� 'DataElements' (�������źų�ʼ��)�� 'GivenMatrix' (������initialDictionary��ʼ��).
    Dictionary(:,1:param.K-param.preserveDCAtom) = Data(:,1:param.K-param.preserveDCAtom);
elseif (strcmp(param.InitializationMethod,'GivenMatrix'))
    Dictionary(:,1:param.K-param.preserveDCAtom) = param.initialDictionary(:,1:param.K-param.preserveDCAtom);
end
% reduce the components in Dictionary that are spanned by the fixed
% elements
if (param.preserveDCAtom)
    tmpMat = FixedDictionaryElement \ Dictionary;
    Dictionary = Dictionary - FixedDictionaryElement*tmpMat;
end
%normalize the dictionary.
Dictionary = Dictionary*diag(1./sqrt(sum(Dictionary.*Dictionary)));
Dictionary = Dictionary.*repmat(sign(Dictionary(1,:)),size(Dictionary,1),1); % multiply in the sign of the first element.
totalErr = zeros(1,param.numIteration);

% the K-SVD algorithm starts here.

for iterNum = 1:param.numIteration
    % find the coefficients
    CoefMatrix = l1([FixedDictionaryElement,Dictionary],Data, param.errorGoal);
    param.L = 1;
%             myDic  = [FixedDictionaryElement,Dictionary];
%         numOfAtoms = size(myDic,2);
%         numOfSignals = size(Data,2);
%         mydim=size(myDic,1);
%         CoefMatrix = zeros([numOfAtoms,numOfSignals]);
%         parfor iii = 1:numOfSignals
%             [s, err_mse, iter_time]=greed_omp_qr(Data(:,iii),myDic,numOfAtoms);
%             CoefMatrix(:,iii)=s';
%         end
%         param.L = 1;

    
%     if (0)%(param.errorFlag==0)
%         %CoefMatrix = mexOMPIterative2(Data, [FixedDictionaryElement,Dictionary],param.L);
%         CoefMatrix = OMP([FixedDictionaryElement,Dictionary],Data, param.L);
%         
% 
%         
%     else 
%         %CoefMatrix = mexOMPerrIterative(Data, [FixedDictionaryElement,Dictionary],param.errorGoal);
%         CoefMatrix = OMPerr([FixedDictionaryElement,Dictionary],Data, param.errorGoal);
%         param.L = 1;
%     end
    
    replacedVectorCounter = 0;
	rPerm = randperm(size(Dictionary,2));%�ֵ���20*50��  ������50������
    for j = rPerm
        [betterDictionaryElement,CoefMatrix,addedNewVector] = I_findBetterDictionaryElement(Data,...
            [FixedDictionaryElement,Dictionary],j+size(FixedDictionaryElement,2),...
            CoefMatrix ,param.L);%�������ź�Y���ֵ�D������j��ϵ��X������L
        Dictionary(:,j) = betterDictionaryElement;
        if (param.preserveDCAtom)
            tmpCoef = FixedDictionaryElement\betterDictionaryElement;
            Dictionary(:,j) = betterDictionaryElement - FixedDictionaryElement*tmpCoef;
            Dictionary(:,j) = Dictionary(:,j)./sqrt(Dictionary(:,j)'*Dictionary(:,j));
        end
        replacedVectorCounter = replacedVectorCounter+addedNewVector;
    end

    if (iterNum>1 & param.displayProgress)
        if (param.errorFlag==0)
            output.totalerr(iterNum-1) = sqrt(sum(sum((Data-[FixedDictionaryElement,Dictionary]*CoefMatrix).^2))/prod(size(Data)));
            disp(['Iteration   ',num2str(iterNum),'   Total error is: ',num2str(output.totalerr(iterNum-1))]);
        else
            output.numCoef(iterNum-1) = length(find(CoefMatrix))/size(Data,2);
            disp(['Iteration   ',num2str(iterNum),'   Average number of coefficients: ',num2str(output.numCoef(iterNum-1))]);
        end
    end
    if (displayErrorWithTrueDictionary ) 
        [ratio(iterNum+1),ErrorBetweenDictionaries(iterNum+1)] = I_findDistanseBetweenDictionaries(param.TrueDictionary,Dictionary);
        disp(strcat(['Iteration  ', num2str(iterNum),' ratio of restored elements: ',num2str(ratio(iterNum+1))]));
        output.ratio = ratio;
    end
    Dictionary = I_clearDictionary(Dictionary,CoefMatrix(size(FixedDictionaryElement,2)+1:end,:),Data);
    
    if (isfield(param,'waitBarHandle'))
        waitbar(iterNum/param.counterForWaitBar);
    end
end

output.CoefMatrix = CoefMatrix;
Dictionary = [FixedDictionaryElement,Dictionary];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  findBetterDictionaryElement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [betterDictionaryElement,CoefMatrix,NewVectorAdded] = I_findBetterDictionaryElement(Data,Dictionary,j,CoefMatrix,numCoefUsed)
if (length(who('numCoefUsed'))==0) %�ж��Ƿ���ڴ˲���
    numCoefUsed = 1;
end
relevantDataIndices = find(CoefMatrix(j,:)); %�ҳ�ϵ��X�ĵ�j�з���ֵ���±�
if (length(relevantDataIndices)<1) %��ǰ��ȫΪ0
    ErrorMat = Data-Dictionary*CoefMatrix;%�����ʾ���
    ErrorNormVec = sum(ErrorMat.^2);%����ÿһ�е������
    [d,i] = max(ErrorNormVec);
    betterDictionaryElement = Data(:,i);%ErrorMat(:,i); %
    betterDictionaryElement = betterDictionaryElement./sqrt(betterDictionaryElement'*betterDictionaryElement);
    betterDictionaryElement = betterDictionaryElement.*sign(betterDictionaryElement(1));
    CoefMatrix(j,:) = 0;
    NewVectorAdded = 1;
    return;
end

NewVectorAdded = 0;
tmpCoefMatrix = CoefMatrix(:,relevantDataIndices); 
tmpCoefMatrix(j,:) = 0;% the coeffitients of the element we now improve are not relevant.
errors =(Data(:,relevantDataIndices) - Dictionary*tmpCoefMatrix); % vector of errors that we want to minimize with the new element
% % the better dictionary element and the values of beta are found using svd.
% % This is because we would like to minimize || errors - beta*element ||_F^2. 
% % that is, to approximate the matrix 'errors' with a one-rank matrix. This
% % is done using the largest singular value.
[betterDictionaryElement,singularValue,betaVector] = svds(errors,1); %�õ����õ��ֵ�Ԫ�ط��ڵ�j��
CoefMatrix(j,relevantDataIndices) = singularValue*betaVector';% *signOfFirstElem %���µ�j��ϵ��

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  findDistanseBetweenDictionaries
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ratio,totalDistances] = I_findDistanseBetweenDictionaries(original,new)
% first, all the column in oiginal starts with positive values.
catchCounter = 0;
totalDistances = 0;
for i = 1:size(new,2)
    new(:,i) = sign(new(1,i))*new(:,i);
end
for i = 1:size(original,2)
    d = sign(original(1,i))*original(:,i);
    distances =sum ( (new-repmat(d,1,size(new,2))).^2);
    [minValue,index] = min(distances);
    errorOfElement = 1-abs(new(:,index)'*d);
    totalDistances = totalDistances+errorOfElement;
    catchCounter = catchCounter+(errorOfElement<0.01);
end
ratio = 100*catchCounter/size(original,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  I_clearDictionary
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Dictionary = I_clearDictionary(Dictionary,CoefMatrix,Data)
T2 = 0.99;
T1 = 3;
K=size(Dictionary,2);
Er=sum((Data-Dictionary*CoefMatrix).^2,1); % remove identical atoms
G=Dictionary'*Dictionary; G = G-diag(diag(G));
for jj=1:1:K,
    if max(G(jj,:))>T2 | length(find(abs(CoefMatrix(jj,:))>1e-7))<=T1 ,
        [val,pos]=max(Er);
        Er(pos(1))=0;
        Dictionary(:,jj)=Data(:,pos(1))/norm(Data(:,pos(1)));
        G=Dictionary'*Dictionary; G = G-diag(diag(G));
    end;
end;

