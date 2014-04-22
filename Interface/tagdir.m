function fns = tagdir(tfrmt)
%
% fns = tagdir(tfrmt)
%
% description:
%     list all files accroding to the tag format trmt


%generate regular expression
[tnames, tagw] = tagformat2tagnames(tfrmt);
re = tfrmt;
for i = 1:length(tnames)
   if tagw(i) == 0 
      repl = regexp(tfrmt, ['<\s*?' tnames{i} '\s*?>'],'match');
      if length(repl) ~= 1
         error('tagdir: cannot infer tag values')
      end
      re = strrep(re, repl{1}, ['(?<' tnames{i} '>\d+?)']);
   else
      repl = regexp(tfrmt,  ['<\s*?' tnames{i}  '\s*?,\s*?' num2str(tagw(i)) '\s*?>'], 'match');
      if length(repl) ~= 1
         error('tagdir: cannot infer tag values')
      end
      re = strrep(re, repl{1}, ['(?<' tnames{i} '>\d{' num2str(tagw(i)) '})']);
   end
end

fprintf('tagdir: seraching for filenames: %s\n', re);

%file names
fns = fileparts(tfrmt);
if isempty(fns)
   fns = '.';
end
fns = dir(fns);
fns = fns(~[fns.isdir]);
fns = {fns.name};

fpath = fileparts(tfrmt);
for i = 1:length(fns)
   fns{i} = fullfile(fpath, fns{i});
end

fns = fns(cellfun(@(x) ~isempty(regexp(x, re, 'once')), fns, 'UniformOutput', true));

end



