function names = dirr(foldername, varargin)
%
% names = dirr(foldername, param)
%
% description:
%     returns absolute filenames of all fildes in the folder/file specification foldername
%
% input:
%     foldername    name of folder and files that can contain * for places holders in filenames as well as directories
%     param         parameter structure with entries
%                   .directories    include directory names in result (false)
%                   .files          include files in result (true)
%                   .recursive      propagate into sub directories after final folder is matched
%                   .dots           include '.' and '..' in dirctories (false)
%
% output:
%     fnames        all filenames that match foldername


% split into subdirectories

param = parseParameter(varargin{:});

includedirs = getParameter(param, 'directories', false);
includefiles = getParameter(param, 'files', true);
recursive = getParameter(param, 'recursive', false);
dots = getParameter(param, 'dots', false);
%apath = getParameter(param, 'absolutepath', false);


names = {};

if ~includedirs && ~includefiles
   return
end
if isempty(foldername)
   return
end


folders = strsplit(foldername, filesep);

%if recursive
%   folders = [folders, {'*'}];
%end

% walk through directories

if strcmp(foldername(1), filesep)
   folders = folders(2:end);
   folders{1} = [filesep, folders{1}];
end
   
dirs = {''};
for f = 1:length(folders)-1
   % check for place holder
   if isempty(strfind(folders{f},'*'))
      for d = 1:length(dirs)
         %check if folder exists
         dirname = fullfile(dirs{d}, folders{f});
         if isdir(dirname)
            dirs{d} = dirname;
         else
            dirs{d} = '';
         end
      end
      dirs = dirs(cellfun(@(x) ~isempty(x), dirs));
   else
      newdirs = {};
      for d = 1:length(dirs)
         dirname = fullfile(dirs{d}, folders{f});
         newsubdirs = dir(dirname);
         newsubdirs = newsubdirs([newsubdirs.isdir]);
         newsubdirs = {newsubdirs.name};  
         
         if ~dots
            newsubdirs = newsubdirs(cellfun('isempty',regexp(newsubdirs,'\<\.', 'match')));
         end
         
         newsubdirs = cellfun(@(x) fullfile(dirs{d}, x), newsubdirs, 'UniformOutput', false);
         newdirs = [newdirs, newsubdirs]; %#ok<AGROW>
      end
      dirs = newdirs;
   end
   
   if isempty(dirs)
      return
   end
end


% find files
if recursive
   names = findfiles(dirs, folders{end});
else
   for d = 1:length(dirs)
      na = fullfile(dirs{d}, folders{end});
      
      if isempty(strfind(folders{end},'*'))
         if includefiles && isfile(na)
            names = [names, {na}]; %#ok<AGROW>
         end
         if includedirs && isdir(na)
            names = [names, {na}]; %#ok<AGROW>
         end
      else   
         newfiles = dir(na);
         
         if ~includefiles % get rid of file if not wanted
            newfiles = newfiles([newfiles.isdir]);
         end
         if ~includedirs % get rid of dirs if not wanted
            newfiles = newfiles(~[newfiles.isdir]);
         end
         newfiles = {newfiles.name};
         if ~dots % get rid of dots
            newfiles = newfiles(cellfun('isempty',regexp(newfiles,'\<\.', 'match')));
         end
         newfiles = cellfun(@(x) fullfile(dirs{d}, x), newfiles, 'UniformOutput', false);
         names = [names, newfiles]; %#ok<AGROW>
      end
   end
end

end
   

% recursive file finder
function newfiles = findfiles(dirs, filepattern)

   newfiles = {};

   for d = 1:length(dirs)
      ndirs = dir(fullfile(dirs{d}, '*'));
      ndirs = ndirs([ndirs.isdir]);
      
      nfiles = dir(fullfile(dirs{d}, filepattern));
      nfiles = nfiles(~[nfiles.isdir]);
         
      ndirs = {ndirs.name};
      ndirs = ndirs(3:end);     % remove '.' and '..'

      nfiles = {nfiles.name};
      
      ndirs = cellfun(@(x) fullfile(dirs{d}, x), ndirs, 'UniformOutput', false);
      nfiles = cellfun(@(x) fullfile(dirs{d},x), nfiles, 'UniformOutput', false);
            
      newfiles = [newfiles, nfiles]; %#ok<AGROW>


      nfiles  = findfiles(ndirs, filepattern);
      newfiles = [newfiles, nfiles]; %#ok<AGROW>

   end
end


