function status= writeMRCfile(a, fname,type)
% type is float or short
[fid,message]=fopen(fname,'w');
if fid == -1
    error('can''t open file');
    a= -1;
    return;
end
ndims=size(size(a));
ndim=ndims(2);
dims(1:256)=1;
if ndims > 3
    ferror('too many dimensions');
end
b=size(a);
dims(1:ndim)=b(1:ndim);
if type=='float'
    dims(4)=2;
    itype='float32';
end
if type =='short'
    dims(4)=1;
    itype='uint16';
end
status=fwrite(fid,dims(1:256),'int32');
status=fwrite(fid,a,itype);
fclose(fid);