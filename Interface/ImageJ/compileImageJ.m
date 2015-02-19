function compileImageJ(force)
%
% compileImageJ(force)
%
% description:
%    compiles MImageJ.java
%
% input:
%    force   force compilation
%+

% check if we need to compile

mijname = fullfile(fileparts(which(mfilename)), 'MImageJ.java');
target = fullfile(fileparts(which(mfilename)), 'MImageJ.class');


if (nargin < 1 || ~force) && exist(target, 'file')   
   fprintf('compileImageJ: %s already exists!\n', target) 
   return;
end


% compiles code for ImageJ interface
[ipath, ijar] = ijpath();

system(['export CLASSPATH=' ipath]);
res = system(['javac -cp ' ijar ' ' mijname]);

if (res ~=0)
   fprintf('compileImageJ: compilation of %s failed!\n', mijname)
else
   fprintf('compileImageJ: compilation of %s succeeded!\n', mijname)
end


end


