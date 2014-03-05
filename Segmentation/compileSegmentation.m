function compileSegmentation()
%
% compileSegmentation()
%
% description:
%    compiles all files necessary for segmentation algorithms
%
%

clc
cd Segmentation

try

mex segmentByPropagationMEX.cpp
mex segmentByPropagation3DMEX.cpp

mex seedPropagationMEX.cpp
mex seedPropagation3DMEX.cpp

catch
end

cd ..

end


