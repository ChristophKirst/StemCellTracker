function sort2subfolder(dirn, pattern, method)
%
% sort2subfolder(dirn, pattern, method)
%
% description:
%    takes images and sortes them into subfolders deterimined by the pattern
%    use regular expression to define patterns and subfolder names:
%    pattern examples:
%    '\w*p(?<position>\d+)\w*.tif' -> folders will be positionXX and files containg pXX will be moved in there
%    '\w*t(?<time>\d+)_p(?<tile>\d+)\w*' -> folders will be time01_tile01 etc and files with t01_p01 will be moved there, etc
%    some standard examples
%    '\w*p(?<p>\d+)\w*.tif'
%
% input:
%    dirn          name of directory
%    pattern       regular expression defining the pattern for subfolders
%    method        (optional) 'move' 'copy' 'test' ('copy')
%
% See also: regexp

if nargin < 3
   method = 'copy';
end

switch method
   case 'copy'
      mode = 1;
   case 'move'
      mode = 2;
   case 'test'
      mode = 3;
   otherwise
      error('sort2subfolder: method must be copy, move or test');
end
   
% get files in folder
fns = dir(dirn);
fns = fns(~[fns.isdir]);
fns = {fns.name};

if isempty(fns)
   warning('sort2subfolder: no filnames to sort in %s', dirn);
   return
end

% get files matching pattern
match = regexp(fns, pattern, 'names');
matchid = ~cellfun(@isempty, match);
fns = fns(matchid);
match = match(matchid);

% move files
dirnames = {};
for f = 1:length(fns)
   dn = dirname(match{f}, dirn);
   if isempty(find(strcmp(dirnames, dn), 1))
      dirnames = [dirnames, dn]; %#ok<AGROW>
      if ~isdir(dn)
         mkdir(dn)
      end
   end
   
   source = fullfile(dirn, fns{f});
   dest = fullfile(dn, fns{f});
   
   switch mode
      case 1
         fprintf('move2subfolder: copying %s -> %s\n', source, dest);
         copyfile(source,dest);
      case 2
         fprintf('move2subfolder: moving %s -> %s\n', source, dest);
         movefile(source, dest);
      case 3
         fprintf('move2subfolder: %s -> %s\n', source, dest);
   end
   
end

end


% sub functions
function dn = dirname(match, path)
   fn = fieldnames(match);
   dn = '';
   for i = 1:length(fn)
      dn = [dn fn{i} match.(fn{i}) '_']; %#ok<AGROW>
   end
   dn = dn(1:end-1);
   dn = fullfile(path, dn);
end



