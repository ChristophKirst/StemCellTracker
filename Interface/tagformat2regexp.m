function tfrmt = tagformat2regexp(tfrmt)
%
% returns regular expression for a tag format

[tnames, tsplit, tinfo] = tagformat2tagnames(tfrmt);

% replace . with \.
tsplit = strfun(@(x) strrep(x, '.', '\.'), tsplit);
res = repmat({''}, length(tsplit)-1, 1);

for i = 1:length(tnames)
   for k = 1:length(tinfo(i).tag)
      if tinfo(i).type == 'd'
         if tinfo(i).width(k) == 0
            re = ['(?<' tnames{i} '>\d*?)'];
         else
            re = ['(?<' tnames{i} '>\d{' num2str(tinfo(i).width(k)) '})'];
         end
      else
         if tinfo(i).width(k) == 0
            re = ['(?<' tnames{i} '>\w*?)'];
         else
            re = ['(?<' tnames{i} '>\w{' num2str(tinfo(i).width(k)) '})'];
         end
      end
      
      res{tinfo(i).pos(k)} = re;
   end  
end

tfrmt = tsplit{1};
for i = 1:length(res);
   tfrmt = [tfrmt,  res{i}, tsplit{i+1}]; %#ok<AGROW>
end

end