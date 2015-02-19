function compileVignetting()
%
% compileVignetting()
%
% description:
%    compiles all files necessary for vignetting algorithms
%

oldpath = pwd;
path = fileparts(mfilename('fullpath'));

cd(path)

try
   clc
   mex mexResponseTransform.cpp

catch
   cd(oldpath)
   error('compileVignetting: error while compiling code in ./Photmetry/Vignetting');
end

cd(oldpath)

end

