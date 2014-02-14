function success = javaaddjar(jarname, search)
%
% javaaddjar(directory, search)
%
% description:
%    adds jar file jarname to java class path
%
% input:
%    jarname    name of jar file or path to jar files
%    search     (optional) if 'all' add all jars in jarname and all its subdirectories
%
% output:
%   success     true if at least one jar file was added
%
% See also: javaaddpath

success = 0;
cpath = javaclasspath('-all');

if nargin < 2
   search = '';
end

if iscell(jarname)
   for i = 1: length(jarname)
      success = succes || javaaddjar(jarname{i}, search);
   end
end

if ~isequal(search, 'all')  % search for single file
   if exist(jarname) == 2    %#ok<EXIST>
      success = true;
      javaaddjarcheck(jarname, cpath);
   end
   
else   % search recursive
     
   % browse to subdirectories
   if isdir(jarname)
      directory = jarname;
   else
      directory = fileparts(jarname);
   end
   
   fns = dir(jarname);
   for i = 1:length(fns)
      if fns(i).isdir && ~strcmp(fns(i).name, '.') && ~strcmp(fns(i).name, '..')
         success = success || javaaddjar(fullfile(directory, fns(i).name));
      end
   end
   
   % add jar files in current directory
   fns = dir(fullfile(directory, '*.jar'));
   for i = 1:length(fns)
      if ~fns(i).isdir
         success = true;
         javaaddjar(fullfile(directory, fns(i).name));
      end
   end
   
end

end
    

 
%%% helper
function javaaddjarcheck(jarfile, cpath)
   %apath = absolutepath(jarfile);
   %[~, fname, ext] = fileparts(jarfile);
   %jarfile = [fname ext];
   
   if isempty(cell2mat(regexp(cpath, [jarfile '$'])))
      %javaaddpath(fullfile(apath, jarname), '-end');
      javaaddpath(jarfile, '-end');
   end
end