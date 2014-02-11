function [X]=OMP(D,Y,L); 
%=============================================
%       D - 字典，列已经归一化
%       Y -目标信号
%       L - 最大非零值个数
%       X - 系数矩阵
%=============================================
[d,N]=size(Y);
[d,K]=size(D);
for i=1:1:N,
    y=Y(:,i);%一次处理一列信号
    indx=zeros(L,1);%最多3个非零值X
    
    
    residual=y;%残差
    a=[];
    for j=1:1:L, 
        proj=abs(D'*residual);
        [maxVal,pos]=max(proj);pos=pos(1); %找到使得内积最大的L列并记录
        indx(j)=pos;
        
        a=pinv(D(:,indx(1:j)))*y; %a是D的前j列为字典时，y的系数
        residual=y-D(:,indx(1:j))*a;%求出残差
        if sum(residual.^2) < 1e-6 %残差够小或者找到L个系数就停止
            break;
        end
    end;
    temp=zeros(K,1);
    temp(indx(1:j))=a; %当前列信号对应的系数的系数在a中，非零位置在indx中
    X(:,i)=sparse(temp);
end;
return;
