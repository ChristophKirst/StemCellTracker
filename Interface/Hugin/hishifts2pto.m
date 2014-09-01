function pto = hishifts2pto(shifts, isizes, varargin)
%
% pto = hishifts2pto(isizes, shifts)
%
% description:
%   convert imgs and shifts to a corresponding pto struct
%
% input:
%   isizes  image sizes in pixel units
%   shifts  image shifts in pixel units 
%   param   parameter struct with entries
%           hFoV   horizotal field of view of first image in degree (1)
%
% output:
%   pto     struct with relevant pto file parameters
%
% See also: hipto2shifts

deg2rad = 2*pi/360;
rad2deg = 1/deg2rad;

param = parseParameter(varargin{:});

nimgs = length(shifts(:));
if nimgs < 1
   error('hishifts2pto: no images and shifts!');
end


pto = struct();

is1 = isizes{1};
pto(1).w = is1(2);
pto(1).h = is1(1);
pto(1).TrZ = 0;

% reference for field of view
hFoV = getParameter(param, 'hFoV', 1) * deg2rad; 
vFoV = 2 * atan(tan(hFoV/2)* pto(1).h / pto(1).w);

zref = pto(1).w / 2 / tan(hFoV/2);

% pixel to panorama coordinate system
phfov = 2 * sin(hFoV/2);
pixel2pano = phfov / pto(1).w;

%fprintf('p2p: %g\n',pixel2pano)

% ref point of first image in pano coords
xy0 = - (1 - pto(1).TrZ) * sin([hFoV, vFoV]/2);

for i = nimgs:-1:1
   
   %image size w,h
   is = isizes{i};
   pto(i).w = is(2);
   pto(i).h = is(1);
   
   % field of view
   hFoV = 2 * atan2(pto(i).w /2, zref);
   vFoV = 2 * atan(tan(hFoV/2)* pto(i).h / pto(i).w);
   
   pto(i).v = hFoV * rad2deg; % hFoV in degree

   % transformations 
   pto(i).TrZ = 0;
   
   psh = pixel2pano * shifts{i};  %shifts in panorama units
   psh = circshift(psh, [0,1]);
   
   xyi = - (1 - pto(i).TrZ) * sin([hFoV, vFoV]/2); % origins of centered image i
   txy = xy0 + psh - xyi;  %translation of camera origin

   pto(i).TrX = txy(1);
   pto(i).TrY = txy(2);

end

end






