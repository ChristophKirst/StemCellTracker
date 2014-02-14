function selection = imarisgetcurrentobject(varargin)
% 
% selection = imarisgetcurrentobject(type)
%
% description: 
%   returns current surpass selection and optionally checks for type
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
if isempty(selection)
   fprintf('imarisgetcurrentobject: no object selected in Imaris Surpass Scene!\n');
end

if nargin == 1
    if isempty(selection) || ~isimaristype(selection, type)

       error('imarisgetcurrentobject: current object is not of type %s', type)
    end
end

end
