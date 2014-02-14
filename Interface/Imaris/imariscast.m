function derived = imariscast(varargin)
%
% derived = imariscast(object, type)
%
% description:
%   cast object to its derived Imaris type
% 
% input:
%   object   Imaris::object object
%   type     (optional) try tocast to this type ([] = automatic)
% 
% output:
% 
%   derived  one of the Imaris::object subclasses:
%                                 - Imaris::IClippingPlane
%                                 - Imaris::IDataContainer
%                                 - Imaris::IFilaments
%                                 - Imaris::IFrame
%                                 - Imaris::IDataSet
%                                 - Imaris::IICells
%                                 - Imaris::ILightSource
%                                 - Imaris::IMeasurementPoints
%                                 - Imaris::ISpots
%                                 - Imaris::ISurfaces
%                                 - Imaris::IVolume
%                                 - Imaris::ISurpassCamera
%                                 - Imaris::IImageProcessing
%                                 - Imaris::IFactory
%            or [] if it fails;
%
% See also: isimaristype, imaristype


[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1
   error('imariscast: need object as input parameter');
else
   object = varargin{1};
end


if isempty(object)
   derived = [];
   return
end

factory = imaris.GetFactory();
 
if factory.IsLightSource(object)
    derived = factory.ToLightSource(object);

elseif factory.IsFrame(object)
    derived = factory.ToFrame(object);

elseif factory.IsVolume(object)
    derived = factory.ToVolume(object);
    
elseif factory.IsSpots(object)
    derived = factory.ToSpots(object);
   
elseif factory.IsSurfaces(object)
    derived = factory.ToSurfaces(object);

elseif factory.IsDataSet(object)
    derived = factory.ToDataSet(object);

elseif factory.IsSurpassCamera(object)
    derived = factory.ToSurpassCamera(object);

elseif factory.IsFilaments(object)
    derived = factory.ToFilaments(object);

elseif factory.IsClippingPlane(object)
    derived = factory.ToClippingPlane(object);

elseif factory.IsApplication(object)
    derived = factory.ToApplication(object);

elseif factory.IsMeasurementPoints(object)
    derived = factory.ToMeasurementPoints(object);

elseif factory.IsDataContainer(object)
    derived = factory.ToDataContainer(object);

elseif factory.IsCells(object)
    derived = factory.ToCells(object);

elseif factory.IsFactory(object)
    derived = factory.ToFactory(object);

elseif factory.IsImageProcessing(object)
    derived = factory.ToImageProcessing(object);   
else
    derived = object;
end

if nargin > 1
   type = varargin{1};

   if ~isimaristype(derived, type)
      error('imariscast: failed to convert to type %s', type);
   end
end

end
   


