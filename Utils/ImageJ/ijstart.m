function mij = ijstart(ijpath)
%
% ijstart(ijpath)
%
% description:
%    starts ImageJ interface
%
% input:
%    ijpath    path to ImageJ or Fiji
%



%%% Fiji version
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

% return

%%% version for mij

if nargin < 1
   % if ispc() make itr work for pc, mac and windows 
       
   ijpath = '/usr/share/imagej/';
end

javaaddjar(ijpath, 'all');

mij = ij.ImageJ([], 2);

end




