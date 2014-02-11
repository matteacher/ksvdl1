classdef denosingData
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        imageName
        oriImage
        sigma
        noiseImage
        PSNRIn
        
        PSNROutDCT
        cleanImageDCT
        DicDCT
        PSNROutGlobal
        cleanImageGlobal
        DicGlobal
        PSNROutKSVD
        cleanImageKSVD
        DicKSVD
    end
    
    methods
        
        function obj = denosingData()

        end
        
        function obj=loadImage(obj,imagename)
           obj.imageName = imagename; 
           %obj.oriImage = obj.loadImage();
           [IMin0,pp]=imread(obj.imageName);
            IMin0=im2double(IMin0);
            if (length(size(IMin0))>2)
                IMin0 = rgb2gray(IMin0); %彩色图像转成灰度图来处理
            end
            if (max(IMin0(:))<2) 
                IMin0 = IMin0*255; %二值图转成灰度图
            end 
            obj.oriImage= IMin0;
        end
        
        function obj=makeNoise(obj,sigma)
           obj.sigma=sigma;
           obj.noiseImage = obj.oriImage+sigma*randn(size(obj.oriImage)); 
           obj.PSNRIn = 20*log10(255/sqrt(mean((obj.noiseImage(:)-obj.oriImage(:)).^2)));
        end
        
        
        function obj=DCT(obj,K,varargin)
            [obj.cleanImageDCT,output] = denoiseImageDCT(obj.noiseImage,obj.sigma,K,[varargin]);
            obj.DicDCT = output.D;
            if(length(size(obj.cleanImageDCT))>2)
                obj.PSNROutDCT = zeros([1,size(obj.cleanImageDCT,3)]);
                for i = 1:size(obj.cleanImageDCT,3)
                    IOut  = obj.cleanImageDCT(:,:,i);
                    obj.PSNROutDCT(i) = 20*log10(255/sqrt(mean((IOut(i)-obj.oriImage(:)).^2))); 
                end
            else
                obj.PSNROutDCT = 20*log10(255/sqrt(mean((obj.cleanImageDCT(:)-obj.oriImage(:)).^2))); 
            end
            
        end
        
        function obj=Global(obj,varargin)
            [obj.cleanImageGlobal,output] = denoiseImageGlobal(obj.noiseImage,obj.sigma,[varargin]);
            obj.DicGlobal = output.D;
            
            if(length(size(obj.cleanImageGlobal))>2)
                obj.PSNROutGlobal = zeros([1,size(obj.cleanImageGlobal,3)]);
                for i = 1:size(obj.cleanImageGlobal,3)
                    IOut  = obj.cleanImageGlobal(:,:,i);
                    obj.PSNROutGlobal(i) = 20*log10(255/sqrt(mean((IOut(i)-obj.oriImage(:)).^2))); 
                end
            else
                obj.PSNROutGlobal = 20*log10(255/sqrt(mean((obj.cleanImageGlobal(:)-obj.oriImage(:)).^2))); 
            end
            
            
        end
        
        function obj=KSVD(obj,K,varargin)
            [obj.cleanImageKSVD,output] = denoiseImageKSVD(obj.noiseImage,obj.sigma,K,[varargin]);
            obj.DicKSVD= output.D;
            
             if(length(size(obj.cleanImageKSVD))>2)
                obj.PSNROutKSVD = zeros([1,size(obj.cleanImageKSVD,3)]);
                for i = 1:size(obj.cleanImageKSVD,3)
                    IOut  = obj.cleanImageKSVD(:,:,i);
                    obj.PSNROutKSVD(i) = 20*log10(255/sqrt(mean((IOut(i)-obj.oriImage(:)).^2))); 
                end
            else
                obj.PSNROutKSVD = 20*log10(255/sqrt(mean((obj.cleanImageKSVD(:)-obj.oriImage(:)).^2))); 
            end
            
        end
        
        function obj=d3(obj,K,varargin)
             obj=obj.DCT(K,varargin);
             obj=obj.Global(varargin);
             obj=obj.KSVD(K,varargin);
        end
        
        function obj=showImages2(obj)
            figure();
            subplot(231)
            imshow(obj.oriImage,[]);
            title(['oriImage' ]);
            
            subplot(232)
            imshow(obj.noiseImage,[]);
            title(['noiseImage' ]);

            subplot(234)
            imshow(obj.cleanImageDCT,[]);
            title(['cleanImageDCT' ]);
            
            subplot(235)
            imshow(obj.cleanImageGlobal,[]);
            title(['cleanImageGlobal' ]);
            
            subplot(236)
            imshow(obj.cleanImageKSVD,[]);
            title(['cleanImageKSVD' ]);
        end
        
        function obj=showImages(obj)
            figure();
            subplot(331)
            imshow(obj.oriImage,[]);
            title(['oriImage' ]);
            
            subplot(332)
            imshow(obj.noiseImage,[]);
            title(['noiseImage' ]);
            

            
            subplot(334)
            imshow(obj.cleanImageDCT,[]);
            title(['cleanImageDCT' ]);
            
            subplot(335)
            imshow(obj.cleanImageGlobal,[]);
            title(['cleanImageGlobal' ]);
            
            subplot(336)
            imshow(obj.cleanImageKSVD,[]);
            title(['cleanImageKSVD' ]);
            
            subplot(337)
            imshow(abs(obj.noiseImage-obj.cleanImageDCT));
            title(['cleanVSnoiseDCT' ]);
            
            subplot(338)
            imshow(abs(obj.noiseImage-obj.cleanImageGlobal));
            title(['cleanVSnoiseGlobal' ]);
            
            subplot(339)
            imshow(abs(obj.noiseImage-obj.cleanImageKSVD));
            title(['cleanVSnoiseKSVD' ]);
        end
        
        function obj=sb(obj,bb,x,y,m,n)
            ss=regexp(obj.imageName, '\.', 'split')
            obj.imageName = [ ss{1} '_sb.bmp']
            obj.oriImage = someblock(obj.oriImage,bb,x,y,m,n);
            obj.noiseImage = someblock(obj.noiseImage ,bb,x,y,m,n);
            obj.cleanImageDCT = someblock(obj.cleanImageDCT ,bb,x,y,m,n);
            obj.cleanImageGlobal = someblock(obj.cleanImageGlobal ,bb,x,y,m,n);
            obj.cleanImageKSVD = someblock(obj.cleanImageKSVD ,bb,x,y,m,n);
        end
        
        function obj=i2b(obj,bb,r)
            image2blocks(obj.oriImage,bb,r,[obj.imageName '_blocks/' 'oriImage/'] )
            image2blocks(obj.noiseImage,bb,r,[obj.imageName '_blocks/' 'noiseImage/'] )
            image2blocks(obj.cleanImageDCT,bb,r,[obj.imageName '_blocks/' 'cleanImageDCT/'] )
            image2blocks(obj.cleanImageGlobal,bb,r,[obj.imageName '_blocks/' 'cleanImageGlobal/'] )
            image2blocks(obj.cleanImageKSVD,bb,r,[obj.imageName '_blocks/' 'cleanImageKSVD/'] )
        end
        

        function obj=saveImage(obj,r)
             imwrite(uint8(imresize(obj.noiseImage,r,'nearest')),[obj.imageName '_noiseImage_' num2str(obj.PSNRIn) '.bmp']);
             imwrite(uint8(imresize(obj.cleanImageDCT,r,'nearest')),[obj.imageName '_cleanImageDCT_' num2str(obj.PSNROutDCT) '.bmp']);
             imwrite(uint8(imresize(obj.cleanImageGlobal,r,'nearest')),[obj.imageName '_cleanImageGlobal_' num2str(obj.PSNROutGlobal) '.bmp']);
             imwrite(uint8(imresize(obj.cleanImageKSVD,r,'nearest')),[obj.imageName '_cleanImageKSVD_' num2str(obj.PSNROutKSVD) '.bmp']);
             imwrite(uint8(imresize(abs(obj.noiseImage-obj.cleanImageDCT),r,'nearest')),[obj.imageName '_cleanVSnoiseDCT_'  '.bmp']);
             imwrite(uint8(imresize(abs(obj.noiseImage-obj.cleanImageGlobal),r,'nearest')),[obj.imageName '_cleannVSnoiseGlobal_'  '.bmp']);
             imwrite(uint8(imresize(abs(obj.noiseImage-obj.cleanImageKSVD),r,'nearest')),[obj.imageName '_cleannVSnoiseKSVD_'  '.bmp']);
        end

    end
    
end
% di = demo001()
% disb = di.sb(8,2,2,5,5)
% disb.i2b(8,10)
