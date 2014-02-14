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

pname = fileparts(name);
if isempty(pname)
   pname = name;
end

w = what(pname);

if ~isempty(w)
   apath = w.path;
else
   apath = [];
end

end
   