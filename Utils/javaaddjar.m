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
      success = javaaddjar(jarname{i}, search) || success;
   end
end

if ~isequal(search, 'all')  % search for single file
   if exist(jarname) == 2    %#ok<EXIST>
      javaaddjarcheck(jarname, cpath);
      success = true;
   else
      warning('javaaddjar: could not add %s', jarname)
   end
   
else   % search recursive
     
   % browse to subdirectories
   if isdir(jarname)
      directory = jarname;
   else
      directory = fileparts(jarname);
   end
   
   fns = dir(directory);
   %fprintf('dir: length =%g, jarname= %s\n', length(fns), jarname);
   
   for i = 1:length(fns)
      %fprintf('dir: n=%s dir=%g  full=%s\n', fns(i).name, fns(i).isdir, fullfile(directory, fns(i).name)) 
      if fns(i).isdir && ~strcmp(fns(i).name, '.') && ~strcmp(fns(i).name, '..')
         %fprintf('new subdir: %s\n', fullfile(directory, fns(i).name));
         success = javaaddjar(fullfile(directory, fns(i).name), 'all') || success;
      end
   end
   
   % add jar files in current directory
   fns = dir(fullfile(directory, '*.jar'));
   %fprintf('file: length =%g, directory = %s\n', length(fns), directory);
   
   for i = 1:length(fns)
      %fprintf('file: n=%s dir=%g\n', fns(i).name, fns(i).isdir)
      if ~fns(i).isdir
         javaaddjarcheck(fullfile(directory, fns(i).name), cpath);
         success = true;
      end
   end
   
%    fns = dir(fullfile(directory, '*.class'));
%    for i = 1:length(fns)
%       if ~fns(i).isdir
%          success = true;
%          javaaddjar(fullfile(directory, fns(i).name));
%       end
%    end
      
end

end


%%% helper
function javaaddjarcheck(jarfile, cpath)
   apath = absolutepath(jarfile);
   [~, fname, ext] = fileparts(jarfile);
   jarfile = fullfile(apath, [fname ext]);
   
   if isempty(cell2mat(regexp(cpath, [fname ext '$'])))
      %javaaddpath(fullfile(apath, jarname), '-end');
      fprintf('javaaddjarr: %s\n', jarfile);
      javaaddpath(jarfile, '-end');
   end
end