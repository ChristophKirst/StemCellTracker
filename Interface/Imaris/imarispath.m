function [ipath, exe, server, jar] = imarispath(hintpath)
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

exe = [];
server = [];
jar = [];
  
if nargin > 0
   % see if we can find Imaris using the information in hintpath    
   ipath = findrootpath(hintpath);
   
   if ~isempty(ipath)
      [exe, server, jar] = imarisfilenames(ipath);
      return
   end
end

% is it setup in the javaclasspath already
jpc = javaclasspath('-all');
ipath = strfind(jpc, 'ImarisLib.jar');
ipath = find(~cellfun(@isempty, ipath));

for i = 1:length(ipath)
   
	ipathi = findrootpath(jpc{find(ipath, 1, 'first')});
   if ~isempty(ipathi)
      ipath = ipathi;
      [exe, server, jar] = imarisfilenames(ipath);
      return
   end
end
  
% environment variable
ipath = getenv('IMARISPATH'); 
if ~isempty(ipath)
   ipath = findroootpath(ipath);
   if ~isempty(ipath)
      [exe, server, jar] = imarisfilenames(ipath);
      return
   end
end

% some guessing
if ispc()
   ipath = 'C:\Program Files\Bitplane\';
elseif ismac()
   ipath = '/Applications';
else
   error('imarispath: imaris requires Windows or Mac');
end

ipath = findrootpath(ipath);
if ~isempty(ipath)
   [exe, server, jar] = imarisfilenames(ipath);
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

function [exe, server, jar] = imarisfilenames(ipath)

   if ispc()
      exe = fullfile(ipath, 'Imaris.exe');
      server = fullfile(ipath, 'ImarisServerIce.exe');
      jar = fullfile(ipath, 'XT', 'matlab', 'ImarisLib.jar');
   elseif ismac()
      exe = fullfile(ipath, 'Contents', 'MacOS', 'Imaris');
      server = fullfile(ipath, 'Contents', 'MacOS', 'ImarisServerIce');
      jar = fullfile(ipath, 'Contents', 'SharedSupport', 'XT', 'matlab', 'ImarisLib.jar');
   else
      error('imarisfilenames: Imaris only runs on Windows or Mac OS')
   end

end


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
   %ipath = fullfile(ipath, 'Imaris');

   ipathcomp = strsplit(ipath, filesep);
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
         continue
      end
      
      % specify files assuming imarisPath is installation root dir
      [exe, server, libjar] = imarisfilenames(cpath);
      
      if ~exist(exe, 'file')
         %fprintf('imarisinitialize: could not find the Imaris executable %s.\n', exe);
         continue;
      end
      
      if ~exist(server, 'file')
         %fprintf('imarisinitialize: could not find the ImarisServer executable %s.\n', server);
         continue;
      end

      if ~exist(libjar, 'file')
         %fprintf('imarisinitialize: could not find the ImarisLib jar file %s.\n', libjar);
         continue;
      end
      
      % consistent Imaris installation in cpath    
      return;
      
   end
   
   cpath = [];
end
