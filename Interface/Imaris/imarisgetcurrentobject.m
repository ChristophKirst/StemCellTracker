function selection = imarisgetcurrentobject(varargin)
% 
% selection = imarisgetcurrentobject(type)
%
% description: 
%   returns current surpass selection if of type type
% 
% input:
%   type   (optional) check for type
% 
% output:
%   selection : autocasted selected object or error if not correct type

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1
   type = [];
else
   type = varargin{1};
end

selection = imariscast(imaris.GetSurpassSelection());

if nargin == 2
    if ~isimaristype(selection, type)
        error('imarisgetcurrentobject: current object is not of type %s', type)
    end
end

end
