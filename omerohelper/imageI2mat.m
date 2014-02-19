function [ imageM ] = imageI2mat( session,imageID )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    pixels = imageID.getPrimaryPixels(); % Same thing image1.getPixels(0)
    sizeZ = pixels.getSizeZ().getValue(); % The number of z-sections. 220
    sizeT = pixels.getSizeT().getValue(); % The number of timepoints. 1
    sizeC = pixels.getSizeC().getValue(); % The number of channels. 1
    sizeX = pixels.getSizeX().getValue(); % The number of pixels along the X-axis. 220
    sizeY = pixels.getSizeY().getValue(); % The number of pixels along the Y-axis. 220
    
    imageM = zeros([sizeX sizeY sizeZ]);
    
    for j = 0:sizeZ-1 
        imageM(:,:,j+1) = getPlane(session,imageID,j,0,0);  %plane = getPlane(session, image, z, c, t) 
    end

end

