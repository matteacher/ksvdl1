addpath('3d')
try
    predemo2_3d
catch
   
    system('C:\matlab\zhy\notify.py')
     rethrow(lasterror)
end
system('C:\matlab\zhy\notify.py')