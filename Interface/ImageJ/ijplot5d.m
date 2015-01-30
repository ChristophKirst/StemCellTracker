function ijplot5d(img, varargin)
%
% ijplot5d(img, param)
%
% description:
%   display 5d image in ImageJ, assumes image of format XYZCT
%
% input:
%   img     image(x,y,z,c,t)
%   param   (optional) parameter struct with
%           .title   title of image
%
% Notes: to develop further commands to feed to ImageJ turn on the recorder while doing 
%        things interactively...Plugins -> Macros -> Recorder
%
% modified from EDS 12/2014

param = parseParameter(varargin);
title = getParameter(param, 'title', []);

simg = size(img);
if length(simg) ~= 5
   error('ijplot5d: inconsistent image dimensions: %s\n', var2char(simg));
end

% transfrom image to QXY where Q = C Z T
nz = size(img,3);
nc = size(img,4);
nt = size(img,5);

img = imfrmtReformat(img, 'XYZCT', 'XYCZT');
img = reshape(img, [size(img,1), size(img,2), nz*nc*nt]);
img = imfrmtReformat(img, 'XYZ', 'ZXY');

size(img)

if isempty(title) 
   title = '';
end
imp = MImageJ.createImage(title, 'lpq', img); % fast creation of image data

imp.show();
imp.updateAndDraw();

%if simg(3) == 3
%    MIJ.run('Stack to RGB');
%else % check for nz=1 or nt=1

cmd = sprintf('order=xyczt(default) channels=%d slices=%d frames=%d display=Composite', nc, nz, nt);
MImageJ.run('Stack to Hyperstack...', cmd);

%end
 
