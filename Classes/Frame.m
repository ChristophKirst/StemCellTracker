classdef Frame < handle
%
% Frame is a class that conains data about objects from single time point
%

   properties
      filename = {};   % the file names of the source image
      objects  = [];   % objects in image (e.g. nuclei)  

   end
   
   methods
      function obj = Frame(filename, objects)
      %
      % Image(filename, objects)
      %
      % input:
      %    filename     source image file name
      %    objects      idenified objects in image
      %
         if nargin > 0
            obj.filename = filename;
         end
         if nargin > 1
            obj.objects = objects;
         end
      end
      
      function f = copy(obj)
      % 
      % f = copy(obj)
      %
      % description:
      %    deep copy of the frame and its objects
      %
         f = Frame(obj.filename, obj.objects.copy);
      end
      
 
      function d = dim(obj)
            d = obj(1).objects(1).dim;
      end

      function data = toArray(obj)
      %
      % data = toArray(obj)
      %
      % convert data of all objects in picture to array
      %  
         data = obj.objects.toArray;
      end
           
      function t = time(obj)
      %
      % t = time(obj)
      %
      % output:
      %   t    time of the frame
      %
         if length(obj) == 1
            t = obj.objects(1).time;
         else
            t = cellfun(@(x) x(1).time, {obj.objects});
         end
      end
      
      
      function xyz = r(obj)
      %
      % xyz = r(obj)
      %
      % output:
      %   xyz    coordinates of the objects in the image os column vectors
      %
         if length(obj) > 1 % for array of images
            xyz = cellfun(@(x) [ x.r ], { obj.objects }, 'UniformOutput', false);
         else               % single image
            xyz = [ obj.objects.r ];
         end   
      end

      function vol = volume(obj)
      %
      % vol = volume(obj)
      %
      % output:
      %   vol    volumes of the objects in the image
      %
         if length(obj) > 1 % for array of images
            vol = cellfun(@(x) [ x.volume ], { obj.objects }, 'UniformOutput', false);
         else               % single image
            vol = [ obj.objects.volume ];
         end   
      end
      
      function i = intensity(obj)
      %
      % i = intensity(obj)
      %
      % output:
      %   i    intensities of the objects in the image
      %
         if length(obj) > 1 % for array of images
            i = cellfun(@(x) [ x.intensity ], { obj.objects }, 'UniformOutput', false);
         else               % single image
            i = [ obj.objects.intensity ];
         end   
      end

      function t = type(obj)
      %
      % t = type(obj)
      %
      % output:
      % t    type data of the objects in the image as column vectors
      %
         if length(obj) > 1 % for array of images
            t = cellfun(@(x) [ x.type ], { obj.objects }, 'UniformOutput', false);
         else               % single image
            t = [ obj.objects.type ];
         end   
      end

      function i = id(obj)
      %
      % i = id(obj)
      %
      % output:
      %   i    coordinates of the objects in the image
      %
         if length(obj) > 1 % for array of images
            i = cellfun(@(x) [ x.id ], { obj.objects }, 'UniformOutput', false);
         else               % single image
            i = [ obj.objects.id ];
         end   
      end
      

      function obj = transformCoordinates(obj, R, T, C)
      %
      % obj = transformCoordinates(obj, R, T, C)
      %
      % applies rotation R, scaling C  and translation T to coordinates of objects r
      %  
         obj.objects = obj.objects.transformCoordinates(R,T,C);
         
      end

   end
   
end
