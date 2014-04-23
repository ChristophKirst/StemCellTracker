function fns = tagformat2files(tfrmt)
%
% returns filnames that are in accordance with the tag format

re = tagformat2regexp(tfrmt);
fns = tagformat2filename(tfrmt);

fns = dir(fns);
fns = fns(~[fns.isdir]);
fns = {fns.name};
fpath = fileparts(tfrmt);
for i = 1:length(fns)
   fns{i} = fullfile(fpath, fns{i});
end

fns = fns(cellfun(@(x) ~isempty(regexp(x, re, 'once')), fns, 'UniformOutput', true));

end