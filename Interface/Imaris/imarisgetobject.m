function object = imarisgetobject(varargin)
%
% object = imarisobject(objectname, type, recursive)
%
% description:
%   returns object with name objectname in surpass scene, [] otherwise
% 
% input:
%   objectname   name of object to look for
%   type         (optional) object must be of type
%   recursive    (optional) scan recursively (true)
% 
% output:
%   object       the founc object or [] otherwise
%
% See also: imarisgetobjects, imarisgetcurrentobject

[imaris, varargin, nargin] = imarisvarargin(varargin);

if nargin > 3
    error('imarisgetobjects: 0 to 3 input parameters expected.');
end

if nargin < 1
   object = imarisgetcurrentobject();
   return
end

objectname = varargin{1};

if nargin < 2
   type = [];
else
   type = varargin{2};
end

if nargin < 3
    recursive = 1;
else
   recursive = varargin{3};
end

object = [];

scene = imaris.GetSurpassScene();
if isempty(scene)
    return
end

nchildren = scene.GetNumberOfChildren();
if nchildren == 0
    return
end

found = false;
if isempty(type)
    getChildrenAtLevel(scene);
else
    getFilteredChildrenAtLevel(scene);
end


% recursive scans
   function getChildrenAtLevel(container)

      for i = 1 : container.GetNumberOfChildren()
         
         if found; return; end
         
         child = container.GetChild(i - 1);
         
         cName = child.GetName();
         if strcmp(cName, objectname)
            object = imariscast(imaris, child);
            found  = true;
            return
         end
         
         if recursive == 1 && imaris.GetFactory().IsDataContainer(child) 
            getChildrenAtLevel(imariscast(imaris, child));
         end
     
      end

   end

 
   function getFilteredChildrenAtLevel(container)
      
      for i = 1 : container.GetNumberOfChildren()
    
         if found; return; end
         
         child = container.GetChild(i - 1);

         cName = child.GetName();
         currentChild = imariscast(imaris, child);

         if strcmp(cName, objectname) && isimaristype(imaris, currentChild, type)
            found = true;
            object = currentChild;
            return
         end
    
         if recursive == 1 && imaris.GetFactory().IsDataContainer(child)
            getFilteredChildrenAtLevel(imariscast(imaris, child));
         else

         end
         
      end
      
   end

end


