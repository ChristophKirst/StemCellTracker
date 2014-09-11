function [tnames, tagsplit, taginfo] = tagexpr2tagnames(texpr)
%
% [tnames, tagsplit, taginfo] = tagexpr2tagnames(texpr)
%
% description:
%      return tag information for the tag expression texpr
%
% input:
%    texpr     tagexpr string, e.e. <name,#chars,type>
%
% output:
%    tnames    tag names
%    tsplit    text strings between tags
%    taginfo   info on tags (type, length etc)
%
% See also: tagexpr

tnames = regexp(texpr, '<(?<name>.*?)>', 'names');
tnames = {tnames.name};

%[tids, tide] =  regexp(texpr, '<(?<name>.*?)>');
[torig, tagsplit] = regexp(texpr, '<(?<name>.*?)>', 'match', 'split');


nnames = length(tnames);
twidth = cell(1,nnames);
ttype = cell(1, nnames);

for i = 1:nnames
   nms = strtrim(strsplit(tnames{i}, ','));
   
   if length(nms) < 2
      nms{2} = '0';
   end
   
   if any(ismember({'s', 'd'}, nms{2}))
      if length(nms) > 2
         nms([3,2]) = nms(2:3);
      else
         nms{3} = nms{2};
         nms{2} = '0';
      end
   end
   twidth{i} = nms{2};
   
   if length(nms) < 3
      nms{3} = 'd';
   end
   ttype{i} = nms{3};
   tnames{i} = nms{1};
end


% check syntax
for i = 1:nnames
   if ~any(ismember({'s', 'd'}, ttype{i}))
      warning('tagexpr2tagnames: syntax error: type %s is neither d(igit) nor s(tring)', ttype{i})
      ttype{i} = 'd';
   end
   
   tw = str2double(twidth{i});
   if isnan(tw)
      warning('tagexpr2tagnames: syntax error: tag width %s not a integer number', twidth{i})
      twidth{i} = 0;
   else
      twidth{i} = tw;
   end
end



% handle same tag appearing more than once

tnamesall = tnames;
tnames = {};
taginfo = struct;


j = 1;
for i = 1:nnames
   k = find(strcmp(tnamesall, tnamesall{i}));
   if k(1) == i
      tnames = [tnames, tnamesall(i)]; %#ok<AGROW>
      
      if nargout > 2
         taginfo(j).name = tnamesall{i};
         %taginfo(j).istart = tids(k);
         %taginfo(j).iend = tide(k);
         taginfo(j).tag = torig(k);
         taginfo(j).pos = k;
         taginfo(j).width = [twidth{k}];

         % type should be the same
         taginfo(j).type = ttype{k(1)};
         
         if any(~strcmp(ttype(k), ttype(k(1))))
            error('tagexpr2tagnames: inconsistent types for tagname %s with multiple appearances',  tnamesall{i})
         end
         
         % for strings the length must be the same
         
         if ttype{k(1)} == 's' 
            tw = [twidth{k}];
            if any(tw ~= tw(1))
               error('tagexpr2tagnames: string tag with tagname %s has multiple appearances with different lengths!',  tnamesall{i});
            end
         end
         
      end
      
      j = j + 1;

   end
end


end







