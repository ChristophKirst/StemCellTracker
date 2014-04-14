function cmd = fileformat2filename(fformat, imagespec, formatnames)
%
% cmd = fileformat2filename(fformat, imagespec)
%
% description:
%    generates a filename from a a file format string and image specifications imagespec
%    replacements are done as follows:
%    for a field with name xxx in the struct imagespec the substring 
%       <xxx>    is replaced by the value, 
%       <xxx,k>  is replaced by the value xxx using k digits with trailling zeros
%    for a vector or cell imagespec the substring  
%       <#n>     is replaced by the nth element of  imagespec
%       <#n,k>   is replaced by the value xxx using k digits with trailling zeros
%
% input:
%    fformat     string of the command that opens the image
%    imagespec   struct/vector or cell with the image specifications
%    formatname  (optional) list of names
%
% See also: imageRead num2str0

cmd = fformat;
if nargin < 2 || isempty(imagespec)
   return 
end

if nargin > 2 || ~isempty(formatnames)
   if length(imagespec) ~= length(formatnames) 
      error('fileformat2filename: inconsistent size of image specifications and format names!')
   end
   if ~iscellstr(formatnames)
      error('fileformat2filename: inconsistent size of image specifications and format names!')
   end
   
   if isstruct(imagespec)
      error('fileformat2filename: image specifications should be cell or vector if formatnames are specified!')
   end
             
   if ~iscell(imagespec)
      imagespec = num2cell(imagespec);
   end
   
   imagespec = {formatnames{:}; imagespec{:}};
   imagespec = struct(imagespec{:});
end
   
if isstruct(imagespec)
   fnames = fieldnames(imagespec);
   for i = 1:length(fnames)
      cmd = strrep(cmd, ['<' fnames{i} '>'], num2str(imagespec.(fnames{i})));
      k = regexp(cmd, ['<' fnames{i} '\s*,\s*(?<k>\d+)>'], 'names');
      if ~isempty(k)
         k = {k.k};  
         for l = 1:length(k)
            repl = num2str0(imagespec.(fnames{i}), str2double(k{l}));
            cmd = regexprep(cmd,  ['<' fnames{i}  '\s*,\s*' k{l} '>'], repl);
         end 
      end
   end
     
else
   if ~iscell(imagespec)
      imagespec = num2cell(imagespec);
   end

   n = length(imagespec);
   for i = 1:n
      cmd = strrep(cmd, ['<#' num2str(i) '>'], num2str(imagespec{i}));

      k = regexp(cmd, ['<#' num2str(i) '\s*,\s*(?<k>\d+)>'], 'names');
      if ~isempty(k)
         k = {k.k}; 
         for l = 1:length(k)
            repl = num2str0(imagespec{i}, str2double(k{l}));
            cmd = regexprep(cmd,  ['<#' num2str(i) '\s*,\s*' k{l} '>'], repl);
         end 
      end
   end
end

end
   
   



