function [tnames, tagsplit, taginfo] = tagformat2tagnames(tfrmt)
%
% [tnames, tagsplit, taginfo] = tagformat2tagnames(tfrmt)
%
% description:
%      return tag information for the tagformat
%
% input:
%    tfrmt     tagformat string, e.e. <name,#chars,type>
%
% output:
%    tnames    tag names
%    tsplit    text strings between tags
%    taginfo   info on tags (type, length etc)

tnames = regexp(tfrmt, '<(?<name>.*?)>', 'names');
tnames = {tnames.name};

[tids, tide] =  regexp(tfrmt, '<(?<name>.*?)>');
[torig, tagsplit] = regexp(tfrmt, '<(?<name>.*?)>', 'match', 'split');


nnames = length(tnames);
twidth = cell(1,nnames);
ttype = cell(1, nnames);

for i = 1:nnames
   nms = strtrim(strsplit(tnames{i}, ','));
   
   if length(nms) < 2
      nms{2} = '0';
   end
   
   if length(nms) == 2 && any(ismember({'s', 'd'}, nms{2}))
      nms{3} = nms{2};
      nms{2} = '0';
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
      warning('tagformat2tagnames: syntax error: type %s is neither d(igit) nor s(tring)', ttype{i})
      ttype{i} = 'd';
   end
   
   tw = str2double(twidth{i});
   if isnan(tw)
      warning('tagformat2tagnames: syntax error: tag width %s not a integer number', twidth{i})
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
         taginfo(j).istart = tids(k);
         taginfo(j).iend = tide(k);
         taginfo(j).tag = torig(k);
         taginfo(j).pos = k;
         taginfo(j).width = [twidth{k}];

         % type should be the same
         taginfo(j).type = ttype{k(1)};
         
         if any(~strcmp(ttype(k), ttype(k(1))))
            error('tagforamt2tagnames: inconsistent types for tagname %s with multiple tag appearances',  tnamesall{i})
         end
         
      end
      
      j = j + 1;

   end
end


end







