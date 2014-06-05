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


% initialize the Ilastik classifer
py('eval', 'from IlastikInterface import *')
py('eval', 'ilc = IlastikClassifier()')



