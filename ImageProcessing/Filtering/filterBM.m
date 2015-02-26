function [img, psf] = filterBM(img, varargin)
%
% img = filterBM(img, param)
%
% description:
%   Block Matrix collaborative filtering using external BM3D filter package
% 
% input:
%   img   image
%   param parameter struct with entries
%         .sigma    Std. dev. of the noise (corresponding to intensities in the range [0,255] even if the range of img is [0,1])
%         .profile  'np' --> Normal Profile, 'lc' --> Fast Profile, 'high' --> High Profile (high quality) 'vn' --> automatically enabled for high noise when sigma > 40
%         .print    0 --> do not print output information, 1 --> print information and plot figures
%   .colorspace:    'opp' --> use opponent colorspace, 'yCbCr' --> use yCbCr colorspace
%
% See also: BM3D, CBM3D

param   = parseParameter(varargin);

sigma   = getParameter(param, 'sigma', 50);
profle  = getParameter(param, 'profile', 'np');
prnt    = getParameter(param, 'print', 0);
cspace  = getParameter(param, 'colorspace', 'opp');
ref     = getParameter(param, 'reference', 1);
depth   = getParameter(param, 'depth', 3);

if isempty(ref)
   ref = 1;
end

% check image format

frmt = imfrmtFormat(img);

switch frmt
   case 'XY'
      [psf, img] = BM3D(ref, img, sigma, profle, prnt); 
   case 'XYC'
      [psf, img] = CBM3D(ref, img, sigma, profle, prnt, cspace);
   case 'XYZ'
      %if ref == 1
         [psf, img] = VBM3D(img, sigma, 0, depth, prnt, profle);
      %else
      %   [psf, img] = VBM3D(img, sigma, prnt, profle, ref);
      %end
   case 'XYZC'
      img = imfrmtReformat(img, 'XYZC', 'XYCZ');   
      %if ref == 1
         [psf, img] = CVBM3D(img, sigma, depth, prnt, profle);
      %else
      %   [psf, img] = CVBM3D(img, sigma, prnt, profle, ref);
      %end     
   case 'XYCZ'
      %img = imfrmtReformat(img, 'XYZC', 'XYCZ');   
      if ref == 1
         [psf, img] = CVBM3D(img, sigma, depth, prnt, profle);
      else
         [psf, img] = CVBM3D(img, sigma, depth, prnt, profle, ref);
      end      
      
   otherwise
      error('filterBM: image format %s not supported', frmt)
end
   
end








