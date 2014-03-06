function [ipath, varargout] = ijpath(hintpath)
%
% [ipath, ijjar, viewer3djar] = ijpath(hintpath)
%
% description:
%    tries to locate ImageJ using an optional hint to a path
%
% input:
%    hintpath  (optional) path to ImageJ or Fiji ([] = autodetection)
%
% output:
%    ipath       path to valid ImageJ or Fiji 
%    ijjar       main ij jar file
%    viewer3djar 3D viewer jar file
%
% See also: ijinitialize

found = 0;

% use hint
if nargin == 1
   [ipath, ijjar, viewer3djar] = checkpath(hintpath);
   if ~isempty(ipath)
      found = 1;
   end
end

%edjucated guessing
if ~found
   if isunix()
      hintpath = '/usr/share/imagej/';
   elseif ispc()
      hintpath = 'C:\Program Files\ImageJ';
   elseif ismac()
       hintpath = '/Applications/ImageJ/ImageJ64.app/Contents/Resources/Java/';
   else
      error('ijinitialize: ImageJ runs on Linux, Max or Windows only!');
   end
   
   [ipath, ijjar, viewer3djar] = findij(hintpath);
   if ~isempty(ipath)
      found = 1;
   end
end

%edjucated guessing 2
if ~found && isunix()
   hintpath = '~/programs/ImageJ/';  
   [ipath, ijjar, viewer3djar] = findij(hintpath);
   if ~isempty(ipath)
      found = 1;
   end
end


% guessing for Fiji
if ~found
   if isunix()
      hintpath = '/usr/share/Fiji.app/';
   elseif ispc()
      hintpath = 'C:\Program Files\Fiji.app';
   elseif ismac()
      error('TODO: Fiji path for mac');
   else
      error('ijinitialize: ImageJ runs on Linux, Max or Windows only!');
   end
   
   [ipath, ijjar, viewer3djar] = findij(hintpath);

   if isempty(ipath)
      error('ijpath: cannot find valid ImageJ installation.')
   end
end


if nargout > 1
   varargout{1} = ijjar;
end
if nargout > 2
   varargout{2} = viewer3djar;
end


end




% helper

function [ijjar, viewer3djar] = ijfilepatterns(ipath, version)

   %Fiji or ImageJ
   switch version
      case 'ImageJ'
         ijjar = fullfile(ipath, 'i*j.jar');
         viewer3djar =  fullfile(ipath, 'plugins', '3D', 'ImageJ_3D_Viewer.jar');
      case 'Fiji'
         ijjar = fullfile(ipath, 'jars', 'ij-*.jar');
         viewer3djar = fullfile(ipath, 'plugins', '3D_Viewer.jar');
   
      otherwise
         error('ijfilenamepatterns:  version must be ImageJ or Fiji')
   end

end



%check if the path is a valid ImageJ path based on ij.jar and 3D_Viewer.jar and return their location
function [ijar, v3djar] = checkpath(ijpath)
 
   ijar = [];
   v3djar = [];
   if ~isdir(ijpath)
      return
   end

   versions = {'Fiji', 'ImageJ'};
   for v = 1:length(versions)
      
      [ijjar, viewer3djar] = ijfilepatterns(ijpath, versions{v});
   
      ijfns = dir(ijjar);
      if isempty(ijfns)
         continue
      end
      
      vfns = dir(viewer3djar);
      if ~isempty(vfns)
         ijar = fullfile(absolutepath(ijjar), ijfns(1).name);
         v3djar = fullfile(absolutepath(ijpath), vfns(1).name);
      end
   end
   
end



% find imagej in ipath
function [rpath, ijar, v3djar] = findij(ipath)  

   rpath = [];

   if ~isdir(ipath)
      ipath = fileparts(ipath);
   end
        
   [ijar, v3djar] = checkpath(ipath);
   if ~isempty(ijar)
      rpath = ipath;
      return;
   end
   
   % get absolute path and check all directories below
   ipath = absolutepath(ipath);
   %ipath = fullfile(ipath, 'ImageJ');

   ipathcomp = stringsplit(ipath, filesep);
   for l = length(ipathcomp):-1:1
      ipath = fullfile(ipathcomp{1:l});
      [ijar, v3djar] = checkpath(ipath);
      if ~isempty(ijar)
         rpath = ipath;
         return;
      end
   end

end