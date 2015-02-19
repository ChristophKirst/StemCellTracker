function compileSegmentation()
%
% compileSegmentation()
%
% description:
%    compiles all files necessary for segmentation algorithms
%
%

oldpath = pwd;
path = fileparts(mfilename('fullpath'));

cd(path)

try

   clc
   mex segmentByPropagationMEX.cpp
   mex segmentByPropagation3DMEX.cpp

   mex seedPropagationMEX.cpp
   mex seedPropagation3DMEX.cpp

catch
   cd(oldpath)
   error('compileSegmentation: error while compiling code in ./Segmentation');
end

cd(oldpath)

end




