function fns = tagformat2files(tfrmt)
%
% returns filnames that are in accordance with the tag format

re = tagformat2regexp(tfrmt);
fns = tagformat2filename(tfrmt);

fns = dirr(fns);
fns = fns(cellfun(@(x) ~isempty(regexp(x, re, 'once')), fns, 'UniformOutput', true));

end