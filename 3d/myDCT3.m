function [ DCT ] = myDCT3( input_args )


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

 
M=88;
MM=12;

figure;

number=1;
%subplot(MM,NN,1)
%hold on;
DCT=zeros([M*M*M,MM*MM*MM]);
for u=1:1:MM
    for v=1:1:MM
        for w=1:1:MM
            for i=1:1:M
                for j=1:1:M
                    for k = 1:M
                        f(i,j,k)=cos(pi*(i+0.5)*(u-1)/M)   * cos(pi*(j+0.5)*(v-1)/M) * cos(pi*(k+0.5)*(w-1)/M);
                    end
                end
            end
        end
        
        I=mat2gray(f);
        DCT(:,number) = I(:);
        %subplot(MM,NN,number);
        %imshow(I);
        number=number+1;
    end
end

end



