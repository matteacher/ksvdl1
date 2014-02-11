function [X]=OMP(D,Y,L); 
%=============================================
%       D - �ֵ䣬���Ѿ���һ��
%       Y -Ŀ���ź�
%       L - ������ֵ����
%       X - ϵ������
%=============================================
[d,N]=size(Y);
[d,K]=size(D);
for i=1:1:N,
    y=Y(:,i);%һ�δ���һ���ź�
    indx=zeros(L,1);%���3������ֵX
    
    
    residual=y;%�в�
    a=[];
    for j=1:1:L, 
        proj=abs(D'*residual);
        [maxVal,pos]=max(proj);pos=pos(1); %�ҵ�ʹ���ڻ�����L�в���¼
        indx(j)=pos;
        
        a=pinv(D(:,indx(1:j)))*y; %a��D��ǰj��Ϊ�ֵ�ʱ��y��ϵ��
        residual=y-D(:,indx(1:j))*a;%����в�
        if sum(residual.^2) < 1e-6 %�вС�����ҵ�L��ϵ����ֹͣ
            break;
        end
    end;
    temp=zeros(K,1);
    temp(indx(1:j))=a; %��ǰ���źŶ�Ӧ��ϵ����ϵ����a�У�����λ����indx��
    X(:,i)=sparse(temp);
end;
return;
