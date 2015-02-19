function compileStr2doubleq()
%
% compileStr2doubleq()
%
% description:
%      compiles all routines for the Polygon package
%

oldpath = pwd;
path = fileparts(mfilename('fullpath'));

cd(path)

try
   
   clc
   mex str2doubleq.cpp
   
catch
   cd(oldpath)
   error('compileStr2double: failed!');
end

cd(oldpath)

end
