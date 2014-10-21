function [texpr, tnames, tags] = tagExpression(fname, varargin)
%
% texpr = tagExpression(fname)
% texpr = tagExpression(fname, param)
% [texpr, tnames, tags] = tagExpression(...)
%
% description:
%        tries to infer tagged name format from files specified by fname 
%        image names must have same length of chars
%        example: img_T001_Z01.tif, img_T001_Z02.tif img_T002_Z01.tif, img_T002_Z02.tif -> img_T<tag1, 3>_Z<tag2, 2>.tif
%
% input:
%      fname     directory, variable filename or file list (e.g.  img_T*_Z*.tif)
%      param     (optional) parameter struct with entries:
%                .tagnames  use these names for the tags in order of appearance ([]= unsorted)
%                .reduce    reduce tags if they show 100% correlation (true)
%
% output:
%      texpr     taggedfile format
%      tnames    (optional) tag names
%      tags      (optional) tag struct array with entries .name = index
%
% See also: tagExpressionToRegularExpression

if ischar(fname)
   fns = dirr(fname);
   if isempty(fns)
      error('tagExpression: no files found for %s!', fname);
   end
else
   fns = fname;
end

if ~iscellstr(fns) || isempty(fns)
   error('tagExpression: first argument not valid or empty!');
end

%filenames need to be of same size for this code to work
if any(diff(cellfun(@length, fns)))
   error('tagExpression: filenames not of same length');
end


param = parseParameter(varargin{:});

tnames = getParameter(param, 'tagnames', {});
if ischar(tnames)
   tnames = {tnames};
end
if ~iscellstr(tnames) 
      error('tagExpression: expect cell of tag names for tagnames parameter');
end

reduce = getParameter(param, 'reduce', true);


%get common chars
fnall = fns{1};
for i = 2:length(fns)
   fnall(fnall ~= fns{i}) = '*';
end

%indentify tags
[sp, tags] = strsplit(fnall, '*');

tagnames = cellfunc(@(x) ['tag' num2str(x)], num2cell(1:length(tags)));

texpr = sp{1};
s = length(texpr) + 1;

for i = 2:length(sp)
      
      e = s + length(tags{i-1}) - 1;   
      %check for trailing zeros
      while texpr(end) == '0'
         texpr(end) = [];
         tags{i-1} = [tags{i-1} '*'];
         s = s - 1;
      end
   
      % get all tags and check for string vs number
      tn = strfun(@(x) x(s:e), fns);
      n = true;
      try 
         cellfun(@str2num, tn); 
      catch
         n = false;
      end
      
      if n
         texpr = [texpr, '<' tagnames{i-1}, ',' num2str(length(tags{i-1})), '>' sp{i}]; %#ok<AGROW>
      else
         texpr = [texpr, '<' tagnames{i-1}, ',s,' num2str(length(tags{i-1})), '>' sp{i}]; %#ok<AGROW>
      end
      
      s = e + 1 + length(sp{i});
end


% reduce 
if reduce
   tags = tagExpressionToTags(texpr, fns);
   [tags, texpr] = tagsReduce(tags, texpr);
   tagnames = fieldnames(tags);
end

% rename
if length(tagnames) < length(tnames)
   tnames = tnames(1:length(tagnames));
else
   for i = length(tnames)+1:length(tagnames)
      tnames{i} = ['tag' num2str(i)];
   end
end
texpr = tagExpressionRename(texpr, tagnames, tnames);

% calculate tags
if nargout > 2
   tags = tagExpressionToTags(texpr, fns);
end

end




