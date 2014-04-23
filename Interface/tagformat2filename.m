function fn = tagformat2filename(tfrmt)
%
% returns filname that are in accordance with the tag format

[tnames, ~, ~, torig] = tagformat2tagnames(tfrmt);
fn = tfrmt;
for i = 1:length(tnames)
   fn = strrep(fn, torig{i}, '*');
end

end