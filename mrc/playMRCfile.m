function a = playMRCfile(fname)
%% Setup
% Create a temporary working folder to store the image sequence. 
% readMRCfile readMRCfile (fname)
[fid,message]=fopen(fname,'r');
if fid == -1
    error('can''t open file');
    a = -1;
    return;
end
nx=fread(fid,1,'long');
ny=fread(fid,1,'long');
nz=fread(fid,1,'long');
type= fread(fid,1,'long');

%fprintf(1,'nx= %d ny= %d nz= %d type= %d', nx, ny,nz,type);
% Seek to start of data
status=fseek(fid,1024,-1);
% Shorts
if type== 1
    a=fread(fid,nx*ny*nz,'int16');
end
%floats
if type == 2
    a=fread(fid,nx*ny*nz,'float32');
end
fclose( fid);
if type ~= 0
    a= reshape(a, [nx ny nz]);
end
if nz == 1
    a= reshape(a, [nx ny]);
end

workingDir = tempname;
mkdir(workingDir);
mkdir(workingDir,'images');

for ii = 1:nz
    img = uint8(mat2gray(a(:,:,ii))*255);
    imwrite(img,fullfile(workingDir,'images',sprintf('img%d.jpg',ii)));
end

imageNames = dir(fullfile(workingDir,'images','*.jpg'));
imageNames = {imageNames.name}';

imageStrings = regexp([imageNames{:}],'(\d*)','match');
imageNumbers = str2double(imageStrings);

[~,sortedIndices] = sort(imageNumbers);
sortedImageNames = imageNames(sortedIndices);

outputVideo = VideoWriter(fullfile(workingDir,'shuttle_out.avi'));
outputVideo.FrameRate = 5; % the video frame rate
open(outputVideo);

for ii = 1:length(sortedImageNames)
    img = imread(fullfile(workingDir,'images',sortedImageNames{ii}));
    
    writeVideo(outputVideo,img);
end

close(outputVideo);

shuttleAvi = VideoReader(fullfile(workingDir,'shuttle_out.avi'));

mov(shuttleAvi.NumberOfFrames) = struct('cdata',[],'colormap',[]);
for ii = 1:shuttleAvi.NumberOfFrames
    mov(ii) = im2frame(read(shuttleAvi,ii));
end

set(gcf,'position', [150 150 shuttleAvi.Width shuttleAvi.Height])
set(gca,'units','pixels');
set(gca,'position',[0 0 shuttleAvi.Width shuttleAvi.Height])

image(mov(1).cdata,'Parent',gca);
axis off;

movie(mov,1,shuttleAvi.FrameRate);