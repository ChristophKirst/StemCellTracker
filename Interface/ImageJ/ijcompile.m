function ijcompile(force)
%
% ijcompile(force)
%
% description:
%    compiles MImageJ.java
%
% input:
%    force   force compilation
%

% check if we need to compile

mijname = fullfile(fileparts(which(mfilename)), 'MImageJ.java');


if (nargin < 1 || ~force) && exist(mijname, 'file')   
   fprintf('ijcompile: %s already exists!\n', mijname) 
   return;
end


% compiles code for ImageJ interface
[ipath, ijar] = ijpath();

system(['export CLASSPATH=' ipath]);
res = system(['javac -cp ' ijar ' ' mijname]);

if (res ~=0)
   fprintf('ijcompile: compilation of %s failed!\n', nijname)
else
   fprintf('ijcompile: compilation of %s succeeded!\n', mijname)
end


end


