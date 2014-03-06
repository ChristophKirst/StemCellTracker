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
   apath = strtrim(apath);
   if status
      error('absolutepath: could not determine absolute path!')
   end

   if ~isdir([apath filesep])
      apath = fileparts(apath);
      if ~isdir(apath)
         error('absolutepath: could not determine absolute path!')
      end
   end

else 
    error('absolutepath: unsupported os')
end

end
   