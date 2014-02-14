function type = imaristype(varargin)
% 
% type = imaristype(iobject)
%
% description:
%   inferst the tzpe of the iobject
% 
% input:
%   iobject  any surpass scene iobject
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

[~, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1
   error('imaristype: need iobject as input parameter');
else
   iobject = varargin{1};
end

if isempty(iobject)
    return
end

type = class(iobject);

types = {'Cells', 'ClippingPlane', 'DataSet', 'Filaments', 'Frame',...
   'LightSource', 'MeasurementPoints', 'Spots', 'Surfaces', 'SurpassCamera', 'Volume','Factory'};

itypes = find(cellfun(@(s) strcmp(['Imaris.I' s 'PrxHelper'], type), types),1,'first');

if ~isempty(itypes)
   type = types(itypes);
else
   type = [];
end
   
end
   
   

