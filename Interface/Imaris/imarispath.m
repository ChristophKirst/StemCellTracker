function [ipath, server, exe, jar] = imarispath(hintpath)
%
% ipath = imarispath(ipath)
%
% description:
%    returns a valid Imaris installation path, [] otherwise
%
% input:
%    hintpath   (optional) path to Imaris installation or ImarisLib.jar file
%
% output:
%   ipath        path to ImarisLib.jar or [] if not found
%   server       full name of ImarisServerIce.exe
%   exe          full name of Imaris.exe
%   jar          full name of ImarisLib.jar
%
% See also: imarisinitialize

  
if nargin > 0
   % see if we can find Imaris using the information in hintpath    
   ipath = findrootpath(hintpath);
   
   if ~isempty(ipath)
      return
   end
end

% is it setup in the javaclasspath already
jpc = javaclasspath('-all');
ipath = strfind(jpc, 'ImarisLib.jar');
ipath = ~cellfun(@isempty, ipath)
for i = 1:length(ipath)
   
	ipathi = findrootpath(jpc{find(ipath, 1, 'first')});
   if ~isempty(ipathi)
      ipath = ipathi;
      return
   end
end
  
% environment variable
ipath = getenv('IMARISPATH'); 
if ~isempty(ipath)
   ipath = findroootpath(ipath);
   if ~isempty(ipath)
      return
   end
end

% some guessing
if ispc()
   ipath = 'C:\Program Files\Bitplane\';
elseif ismac()
   ipath = '/Applications';
end

ipath = findroootpath(ipath);
if ~isempty(ipath)
   return
end

% TODO:
% check running processes
% if ispc()
%    
%    [~, result] = system(...
%       'tasklist /NH /FI "IMAGENAME eq ImarisServerIce.exe"');
%    if strfind(result, 'ImarisServerIce.exe')
%       isRunning = 1;
%       return;
%    end
%    
% elseif ismac()
%    
%    [~, result] = system('ps aux | grep ImarisServerIce');
%    if strfind(result, this.mImarisServerExePath)
%       isRunning = 1;
%       return;
%    end

end


%%% helper

% find installation root
function rpath = findrootpath(ipath)  
   if ~isdir(ipath)
      ipath = fileparts(ipath);
   end
        
   rpath = checkpath(ipath);
   if ~isempty(rpath)
      return;
   end
   
   % add Imaris and then browse all possible parent directories
   ipath = absolutepath(ipath);
   ipath = fullfile(ipath, 'Imaris');

   ipathcomp = stringsplit(ipath, filesep);
   for l = length(ipathcomp):-1:1
      ipath = fullfile(ipathcomp{1:l});
      rpath = checkpath(ipath);
      if ~isempty(rpath)
         return;
      end
   end

   rpath = [];
end

function cpath = checkpath(imarisPath)  % check if path or path/Imaris* contains a valid Imaris path

   dirs = dir(fullfile(imarisPath, 'Imaris*'));
   dirs = dirs([dirs.isdir]);
   dirs = {dirs.name};
   for i = 1:length(dirs)
      dirs{i} = fullfile(imarisPath, dirs{i});
   end
   dirs{end+1} = imarisPath;
   
   for i = length(dirs):-1:1 % starting from top will give newest version

      cpath = dirs{i};
      
      if ~isdir(cpath)
         break
      end
      
      % specify files assuming imarisPath is installation root dir
      if ispc()
         exe = fullfile(cpath, 'Imaris.exe');
         server = fullfile(imarisPath, 'ImarisServerIce.exe');
         libjar = fullfile(cpath, 'XT', 'matlab', 'ImarisLib.jar');
      elseif ismac()
         exe = fullfile(cpath, 'Contents', 'MacOS', 'Imaris');
         server = fullfile(cpath, 'Contents', 'MacOS', 'ImarisServerIce');
         libjar = fullfile(cpath, 'Contents', 'SharedSupport', 'XT', 'matlab', 'ImarisLib.jar');
      end
      
      if ~exist(exe, 'file')
         %fprintf('imarisinitialize: could not find the Imaris executable.');
         break;
      end
      
      if ~exist(server, 'file')
         %fprintf('imarisinitialize: could not find the ImarisServer executable.';
         break;
      end

      if ~exist(libjar, 'file')
         %errorMessage = 'imarisinitialize: could not find the ImarisLib jar file.';
         break;
      end
      
      % consistent Imaris installation in cpath    
      return;
      
   end
   
   cpath = [];
end
