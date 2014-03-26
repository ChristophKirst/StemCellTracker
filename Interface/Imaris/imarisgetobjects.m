function  objects = imarisgetobjects(varargin)
%
% object = imarisgetobjects(type, recursive)
%
% description:
%   returns cell array of all objects in surpass scene of type type.
% 
% input:
%   type         (optional) type of object ([] = all objects)
%   recursive    (optional)  scan recursively (true)
% 
% output:
%   objects      cell array of references to objects
%
% See also: imarisgetobject, imarisgetcurrentobject

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin < 1
   type = [];
else
   type = varargin{1};
end

if nargin < 2
    recursive = 1;
else
   recursive = varargin{2};
end

objects = {};

scene = imaris.GetSurpassScene();
if isempty(scene)
    return
end

nchildren = scene.GetNumberOfChildren();
if nchildren == 0
    return
end

currentChildNumber = 0;

if isempty(type)
    getChildrenAtLevel(scene);
else
    getFilteredChildrenAtLevel(scene);
end


% Recursive function to scan children of a given container
   function getChildrenAtLevel(container)
      
      % This function scans the children recursively
      for i = 1 : container.GetNumberOfChildren()
         
         % Get current child
         child = container.GetChild(i - 1);
         
         % Is this a folder? If it is, call this function recursively
         if imaris.GetFactory().IsDataContainer(child)
            if recursive == 1
               getChildrenAtLevel(imariscast(imaris, child));
            end
         else
            currentChildNumber = currentChildNumber + 1;
            objects{currentChildNumber} = imariscast(imaris, child);
         end
         
      end
      
   end


% Recursive function to scan children of a given container
   function getFilteredChildrenAtLevel(container)
      
      % This function scans the children recursively
      for i = 1 : container.GetNumberOfChildren()
         
         % Get current child
         child = container.GetChild(i - 1);
         
         % Is this a folder? If it is, call this function recursively
         if imaris.GetFactory().IsDataContainer(child)
            if recursive == 1
               getFilteredChildrenAtLevel(imariscast(imaris, child));
            end
         else
            currentChild = imariscast(imaris, child);
            if isimaristype(imaris, currentChild, type)
               currentChildNumber = currentChildNumber + 1;
               objects{currentChildNumber} = currentChild;
            end
         end
         
      end
      
   end

end
