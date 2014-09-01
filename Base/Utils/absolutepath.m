function apath = absolutepath(name)
%
% description:
%    absolute path to the file/folder name
%
% input:
%    name    file or folder name
%
% output:
%    apath   path to file/folder name

if ismac() || ispc() % windows / mac
   
   curdir = cd;
   
   if isdir(name)
      pname = name;
   else
      pname = fileparts(name);
   end
   
   try
      cd(pname);
   catch
      apath = [];
      return
   end
   
   apath = pwd;
   cd(curdir); 
    
elseif isunix()
   
   [status, apath] = system(['readlink -f ' name]);
   
   if status
      error('absolutepath: could not determine absolute path for ''%s!''', name)
   end   
   
   %apath = strtrim(apath);
   apath = strsplit(apath, '\n');
   apath = apath{1};

   if ~isdir([apath filesep])
      apath = fileparts(apath);
      if ~isdir(apath)
         error('absolutepath: could not determine absolute path for ''%s!''', name)
      end
   end

else 
    error('absolutepath: unsupported os')
end

end
   