function fn = tagformat2filename(tfrmt)
%
% returns filename pattern that is in accordance with the tag format

[~, tsplit, ~] = tagformat2tagnames(tfrmt);
fn = tsplit{1};
for i = 2:length(tsplit)
   fn = [fn, '*', tsplit{i}]; %#ok<AGROW>
end

end