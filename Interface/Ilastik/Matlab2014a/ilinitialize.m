function ilinitialize(varargin)
%
% ilinitialize()
%
% description:
%     set up the interface to ilasitk via the python interface
%
% input:
%    hintpath  (optional) path to ilastik

%ldpath = getenv('LD_LIBRARY_PATH');
%setenv('LD_LIBRARY_PATH', ''); 

cpath = fileparts(mfilename('fullpath'));
ipath = ilpath(varargin{:});

% set ilastik path
fid = fopen(fullfile(cpath, 'IlastikConfigTemplate.py'),'r');
f = fread(fid,'*char')';
fclose(fid);
 
f = strrep(f, 'ILASTIK_PATH', ipath);

fid  = fopen(fullfile(cpath, 'IlastikConfig.py'),'w');
fprintf(fid,'%s',f);
fclose(fid);

% set current path in python 
py('set', 'currentpath', cpath);
py('eval', 'import os,sys')
py('eval', 'sys.path.append(currentpath)')

% set so that the matlab h5 lib version is not used

if count(py.sys.path, cpath) == 0
    insert(py.sys.path,int32(0), cpath);
end

% initialize the Ilastik classifer
py('eval', 'from IlastikInterface import *')
py('eval', 'ilc = IlastikClassifier()')

end



