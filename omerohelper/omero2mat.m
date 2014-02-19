%导入Omero Matlab库
addpath 'C:\OMEROmatlab'
loadOmero;

%连接Server
client = loadOmero('localhost', 4064);
session = client.createSession('root', 'spc');
client.enableKeepAlive(60);

dataset = getDatasets(session,1);  %ID为1的数据集
for i = 0:0  %dataset.linkedImageList.size -1 
    imageID = dataset.linkedImageList.get(i);
    imageM = imageI2mat(session,imageID);
    mkdir('C:\data')
    save(['C:\data\image_' num2str(i) '.mat'],'imageM');    
end

client.closeSession();