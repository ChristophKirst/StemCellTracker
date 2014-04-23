function [tnames, twidth, ttype, torig] = tagformat2tagnames(tfrmt)

% return tag information

tnames = regexp(tfrmt, '<(?<name>.*?)>', 'names');
tnames = {tnames.name};
nnames = length(tnames);

if nargout > 3
   torig = tnames;
   for i = 1:nnames
      torig{i} = ['<' torig{i} '>'];
   end
end

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
   twidth{i} = str2double(nms{2});
   
   if length(nms) < 3
      nms{3} = 'd';
   end
   ttype{i} = nms{3};
   tnames{i} = nms{1};
end

end







