function mij = ijinitialize(ijpath)
%
% mij = ijinitialize(ijpath)
%
% description:
%    initializes ImageJ interface by making installing the jars
%
% input:
%    ijpath  (optional) path to ImageJ or Fiji ([] = autodetection)
%
% output:
%    mij     ij.ImageJ instance of imagej
%
% See also: ijstart


if nargin < 1
   if isunix()
      ijpath = '/usr/share/imagej/';
   elseif ispc()
      ijpath = 'C:\Program Files\ImageJ';
   elseif ismac()
      error('TODO: path for mac');
   end 
end

if ~checkpath(ijpath)
   wrning('ijinitialize: ijpath might not be correct!');
end

javaaddjar(ijpath, 'all');

if ~exist('ij.ImageJ', 'class')
   error('ijinitialize: failed try to specify correct ijpath!');
end


mij = ij.ImageJ([], 2);


%%% old Fiji version
h% dirn = pwd;
% 
% if nargin < 1
%    ijpath = '/home/ckirst/programs/Fiji/scripts/';
% end
% 
% addpath(ijpath);
% 
% Miji(false); %run startup script from Fiji
% 
% cd(dirn)


end

% helper

function b = checkpath(ijpath)
   % in the path there should be a imagej.jar or ij.jar or a plugins or a jar folder
   b = 0;
   
   fns = dir(ijpath);
   b = b | any(strcmp({fns.name}, 'imagej.jar');
   b = b | any(strcmp({fns.name}, 'ij.jar');
   b = b | any(~cellfun(@isempty, strfind({fns.name}, 'ij')));
   
   fns = fns([fns.isdir]);
   
   b = b | any(strcmp({fns.name}, 'plugins');
   b = b | any(strcmp({fns.name}, 'jars');
   
end