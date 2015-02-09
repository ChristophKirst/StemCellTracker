function [hpath, varargout] = hipath(hintpath)
%
% [hpath, tools] = hipath(hintpath)
%
% description:
%    tries to locate Hugin command line tools
%
% input:
%    hintpath  (optional) path to Hugin
%
% output:
%    hpath    path to valid Hugin tools
%    tools    map of tools to tool commands  {pto_gen, pto_var, enblend, nona, autooptimiser, vig_optimize}
%
% See also: hiinitialize

if ispc()
   error('hipath: hugin tools for windows not supportedhi!')
end


found = false;

% use hint
if nargin == 1
   [hpath, tools] = findhi(hintpath);
   if ~isempty(hpath)
      found = true;
   end
end

%edjucated guessing
if ~found
   if isunix()
      [succ, hintpath] = system('which hugin');
      hintpath = fileparts(hintpath);
      if succ == 0
         [hpath, tools] = findhi(hintpath);
         if ~isempty(hpath)
            found = true;
         end
      end
   end
end

if ~found
    error('hipath: cannot find valid Hugin installation, try passing a path to hipath')
end 

if nargout > 1
   varargout{1} = tools; 
end

end



% helper

%check if the path is a valid Hugin path with the needed tools
function tools = checkpath(hipath)
   tools = containers.Map();
   if ~isdir(hipath)
      return
   end

   bins = {'hugin', 'pto_gen', 'pto_var', 'cpfind', 'enblend', 'nona', 'autooptimiser', 'vig_optimize'};

   for b = 1:length(bins)
      fn = fullfile(hipath, bins{b});
      if ~isfile(fn)
         tools = {};
         fprintf('hipath: cannot find %s\n!', bins{b});
         return
      else
         tools(bins{b}) = fn;
      end
   end
end



% find imagej in ipath
function [rpath, tools] = findhi(hpath)  

   rpath = [];
   
   if ~isdir(hpath)
      hpath = fileparts(hpath);
   end
        
   tools = checkpath(hpath);
   if ~isempty(tools)
      rpath = absolutepath(hpath);
      return;
   end
   
   % get absolute path and check all directories below
   hpath = absolutepath(hpath);
   %ipath = fullfile(ipath, 'ImageJ');

   hpathcomp = strsplit(hpath, filesep);
   for l = length(hpathcomp):-1:1
      hpath = fullfile(hpathcomp{1:l});
      tools = checkpath(hpath);
      if ~isempty(tools)
         rpath = absolutepath(hpath);
         return;
      end
   end

end