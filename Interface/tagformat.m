function [tfrmt, tnames, tags] = tagformat(fname, tagnames)
%
% tfrmt = tagformat(fname)
% tfrmt = tagformat(fname, tagnames)
% [tfrmt, tnames, tags] = tagformat(...)
%
% description:
%        tries to infer tagged name format from files specified by fname 
%        image names must have same length of chars
%        example: img_T001_Z01.tif, img_T001_Z02.tif img_T002_Z01.tif, img_T002_Z02.tif -> img_T<tag1, 3>_Z<tag2, 2>.tif
%
% input:
%      fname     directory or variable filename (e.g.  img_T*_Z*.tif)
%      tagnames  (optional) use these names for the tags in order of appearance
%
% output:
%      tfrmt     taggedfile format
%      tnames    (optional) tag names
%      tags      (optional) indices in each tag dimension
%
% See also: tags2name, tagformat2tagnames, name2tags

if ~ischar(fname)
   error('tagformat: first argument not valid file or directory name');
end

if ~isdir(fname)
   [dirname, fname] = fileparts(fname);
else
   dirname = fname;
   fname = '*';
end

if ~isdir(dirname)
   error('tagformat: %s not a valid directory', dirname);
end

if nargin < 2
   tagnames = {};
else
   if ~iscellstr(tagnames) 
      error('tagformat: expect cell of tag names as second argument');
   end
end


fns = dir(fullfile(dirname, fname));
fns([fns.isdir]) = [];
fns = {fns.name};

if isempty(fns)
   error('tagformat: no files in %s', fullfile(dirname, fname));
end

%filenames need to be of same size
if any(diff(cellfun(@length, fns)))
   error('tagformat: filenames not of same length');
end

%get common chars
fnall = fns{1};
for i = 2:length(fns)
   fnall(fnall ~= fns{i}) = '*';
end

%indentify tags
[sp, tags] = strsplit(fnall, '*');

for i = 1:length(tags)
   if length(tagnames) < i
      tagnames{end+1} = ['tag', num2str(i)]; %#ok<AGROW>
   end
end

tfrmt = sp{1};
for i = 2:length(sp)
      %check for trailing zeros
      while tfrmt(end) == '0'
         tfrmt(end) = [];
         tags{i-1} = [tags{i-1} '*'];
      end
   
      tfrmt = [tfrmt, '<' tagnames{i-1}, ',' num2str(length(tags{i-1})), '>' sp{i}]; %#ok<AGROW>
end

if nargout < 2
   return
end

tnames = tagnames;

if nargout < 3
   return
end

% check for the indices
tags = name2tags(tfrmt, fns, tagnames);
%mintags = min(mintags, tags);
%maxtags = max(maxtags, tags);

%for i = size(tags, 1):-1:1
%   idx{i} = unique(tags(i,:));
%end   

end




