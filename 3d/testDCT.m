function [ DCT ] = testDCT( input_args )

% bb=8;
% Pn=16;
% DCT=zeros(bb,Pn); %8行16列
% for k=0:1:Pn-1,
%     V=cos([0:1:bb-1]'*k*pi/Pn);
%     if k>0, V=V-mean(V); end;
%     DCT(:,k+1)=V/norm(V);
% end;
% DCT=kron(DCT,DCT); %两个8x16的结果是64x256

% DCT coefficient function

 
M=8;N=8;
MM=16;NN=16;

figure;

number=1;
subplot(MM,NN,1)
hold on;
I=[];
for u=1:1:MM
    for v=1:1:NN
        for i=1:1:M
            for j=1:1:N
                f(i,j)=cos(pi*(i+0.5)*(u-1)/M)   * cos(pi*(j+0.5)*(v-1)/N);
            end
        end
        I=mat2gray(f);
        subplot(MM,NN,number);
        imshow(I);
        number=number+1;
    end
end

end

