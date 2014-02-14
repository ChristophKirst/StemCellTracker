function imarisopen(varargin)
%
% imarisopen(filename, opts)
% imarisopen(imaris, filename, opts)
%
% description:
%    opens a file in imaris
%
% input:
%    filename   file to open
%    imaris     (optional) Imaris Application instance
%    opts       (optional) options ('')
%      
% usage:
%    Formats:    reader="Imaris3",
%                Available Formats: Imaris5, Imaris3, Imaris, AndorIQ, Andor, DeltaVision, Biorad, IPLab, 
%                IPLabMac, Gatan, CXD, SlideBook, MRC, LeicaLif, LeicaSingle, LeicaSeries, LeicaVista, 
%                MicroManager, MetamorphSTK, MetamorphND, ICS, NikonND2, OlympusCellR, OlympusOIB, OlympusOIF,
%                Olympus, OlympusVSI, OmeTiff, OmeXml, OpenlabLiff, OpenlabRaw, PerkinElmer2, Prairie,
%                Till, AxioVision, Lsm510, Lsm410, BmpSeries, TiffSeries, ZeissCZI.
%    Cropping:   croplimitsmin="x0 y0 z0 c0 t0"
%                croplimitsmax="x y z c t"   (0 = automatic detection of size)
%    Resampling: resample="rx ry rz rc rt"
%                LoadDataSet="eDataSetYes | eDataSetNo | eDataSetWithDialog" 
%                size of the loaded dataset is "(croplimitsmax - croplimitsmin) / resample" along each dimension
% 
% See also: imarissave

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1
   error('imarisopen: need to sepcify filename')
else
   filename = varargin{1};
end
if nargin < 2
   opts = '';
else
   opts = varargin{2};
end
   
imaris.FileOpen(filename, opts);

end   