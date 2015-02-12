function compileVignetting()
%
% compileVignetting()
%
% description:
%    compiles all files necessary for vignetting algorithms
%
%

clc
cd Photmetry/Vignetting

try

mex mexResponseTransform.cpp

catch
   cd ../..
   error('compileVignetting: error while compiling code in ./Photmetry/Vignetting');
end

cd ..

end


