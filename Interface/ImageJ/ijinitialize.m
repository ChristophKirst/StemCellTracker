function [mij, mimagej] = ijinitialize(varargin)
%
% mij = ijinitialize(ijpath)
%
% description:
%    initializes ImageJ interface by installing the java jars
%
% input:
%    ijpath  (optional) hint to path to ImageJ or Fiji ([] = autodetection)
%
% output:
%    mij     ij.ImageJ instance of imagej
%
% See also: ijstart

if ~exist('ij.ImageJ', 'class')
   % add imagej to java path
   [ipath, ijar, view3djar] = ijpath(varargin{:}); %#ok<ASGLU>
   
   if ~javacheckclasspath(ijar)
        fprintf('ijinitialize: installing %s\n', ijar);
        javaaddpath(ijar, '-end');
   end
   if ~javacheckclasspath(view3djar)
        fprintf('ijinitialize: installing %s\n', view3djar);
        javaaddpath(view3djar, '-end');
   end
   
   %javaaddjar(ipath, 'all');
end

if ~exist('ij.ImageJ', 'class')
   error('ijinitialize: failed try to specify correct ijpath!');
end

% add MImageJ to jave path
mpath = fileparts(which(mfilename));
if ~javacheckclasspath(mpath)
   javaaddpath(fullfile(mpath, 'MImageJ.class') , '-end');
   javaaddpath(mpath, '-end');
end

mij = ijinstance();
if isempty(mij)
   try
      mij = ij.ImageJ([], 2);
   catch %#ok<CTCH>
      error('ijinitialize: error while initializing ImageJ classes')
   end
   %ijinfo();
end

if nargout > 1
   mimagej = MImageJ();
end

fprintf('ijinitialize: ImageJ interface installed using: %s\n', ijpath(varargin{:}));


end

%%% old Fiji version
% dirn = pwd;
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

