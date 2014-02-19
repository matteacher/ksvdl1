%����Omero Matlab��
addpath 'C:\OMEROmatlab'
loadOmero;

%����Server
client = loadOmero('localhost', 4064);
session = client.createSession('root', 'spc');
client.enableKeepAlive(60);

dataset = getDatasets(session,1);  %IDΪ1�����ݼ�
for i = 0:0  %dataset.linkedImageList.size -1 
    imageID = dataset.linkedImageList.get(i);
    imageM = imageI2mat(session,imageID);
    mkdir('C:\data')
    save(['C:\data\image_' num2str(i) '.mat'],'imageM');    
end

client.closeSession();