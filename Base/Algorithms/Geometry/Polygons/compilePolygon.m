function compilePolygon()
%
% compilePolygon()
%
% description:
%      compiles all routines for the Polygon package
%

oldpath = pwd;
path = fileparts(mfilename('fullpath'));

cd(path)

try
   
   mex mexPolygonBuffer.cpp clipper.cpp
   mex mexPolygonExecute.cpp clipper.cpp
   
catch
   cd(oldpath)
   error('compilePolygon: failed!');
end

cd(oldpath)

end

