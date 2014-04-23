function tfrmt = tagformat2regexp(tfrmt)
%
% returns regular expression for a tag format

[tnames, twidth, ttype, torig] = tagformat2tagnames(tfrmt);

% replace . with \.
tfrmt = strrep(tfrmt, '.', '\.');

for i = 1:length(tnames)
   if ttype{i} == 'd'
      if twidth{i} == 0
         re = ['(?<' tnames{i} '>\d*?)'];
      else
         re = ['(?<' tnames{i} '>\d{' num2str(twidth{i}) '})'];
      end
   else
      if twidth{i} == 0
         re = ['(?<' tnames{i} '>\w*?)'];
      else
         re = ['(?<' tnames{i} '>\w{' num2str(twidth{i}) '})'];
      end
   end  
   tfrmt = strrep(tfrmt, torig{i}, re);
end

end