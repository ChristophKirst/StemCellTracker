% compile py.cpp

clc

if isdir('/tmp/matpy')
   rmdir('/tmp/matpy', 's');
end
delete('./Interface/Python/py.mex*')


if isunix() && ~ismac
   % matlabs mex funciotn is buggy so we compile directly
   mkdir('/tmp/matpy');
   %compile
   [stat, cmdout] = system('/usr/bin/g++ -c -D_GNU_SOURCE -DMATLAB_MEX_FILE -I/usr/lib64/python2.7/site-packages/numpy/core/include/ -I"/usr/local/MATLAB/R2014a/extern/include" -I"/usr/local/MATLAB/R2014a/simulink/include" -ansi -fexceptions -fPIC -fno-omit-frame-pointer -pthread -O -DNDEBUG ./Interface/Python/py.cpp -o /tmp/matpy/py.o')
   %link
   [stat, cmdout] = system('/usr/bin/g++ -pthread -Wl,--no-undefined  -shared -O -Wl,--version-script,"/usr/local/MATLAB/R2014a/extern/lib/glnxa64/mexFunction.map" /tmp/matpy/py.o  -lpython2.7   -L/usr/lib64   -Wl,-rpath-link,/usr/local/MATLAB/R2014a/bin/glnxa64 -L"/usr/local/MATLAB/R2014a/bin/glnxa64" -lmx -lmex -lmat -lm -lstdc++ -o ./Interface/Python/py.mexa64')

   rmdir('/tmp/matpy', 's')
else
   
  % under mac / windows try if this works -> not tested yet
  mex py.cpp -lpyhton2.7
end
