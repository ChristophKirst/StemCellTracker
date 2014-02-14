function type = imaristype(varargin)
% 
% type = imaristype(object)
%
% description:
%   inferst the tzpe of the object
% 
% input:
%   object  any surpass scene object
%
% output:
%   type    one of:
%               'Cells'
%               'ClippingPlane'
%               'Dataset'
%               'Filaments'
%               'Frame'
%               'LightSource'
%               'MeasurementPoints'
%               'Spots'
%               'Surfaces'
%               'SurpassCamera'
%               'Volume'
%               [] for unknwon
%
% See also: isimaristype, imariscast

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1
   error('imaristype: need object as input parameter');
else
   object = varargin{1};
end

factory = imaris.GetFactory();
 
type = [];

if isempty(object)
    return
   
elseif factory.IsCells(object)
   type = 'Cells';
   
elseif factory.IsClippingPlane(object)
   type = 'ClippingPlane';
   
elseif factory.IsDataSet(object)
   type = 'DataSet';
   
elseif factory.IsFilaments(object)
   type = 'Filaments';
   
elseif factory.IsFrame(object)
   type = 'Frame';
   
elseif factory.IsLightSource(object)
   type = 'LightSource';
   
elseif factory.IsMeasurementPoints(object)
   type = 'MeasurementPoints';
   
elseif factory.IsSpots(object)
   type = 'Spots';
   
elseif factory.IsSurfaces(object)
   type = 'Surfaces';
   
elseif factory.IsSurpassCamera(object)
   type = 'SurpassCamera';
   
elseif factory.IsVolume(object)
   type = 'Volume';
   
elseif factory.IsFactory(object)
   type = 'Factory';
   
end

end
   
   

