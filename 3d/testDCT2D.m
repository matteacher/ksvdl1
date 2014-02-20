function [ myDCT ] = testDCT2D( input_args )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
dim = 2;
bb =8;
myDCT = zeros([bb bb]);
for i = 1:bb
    for j = 1:bb
        sum=0;
        for ii = 0:i-1
           for jj = 0:j-1
              sum=sum+cos(pi*(2*ii+1))/16*i * cos(pi*(2*jj+1))/16*j; 
           end
        end
        myDCT(i,j) = sum;
    end
end

