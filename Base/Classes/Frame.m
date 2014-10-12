classdef Frame < handle
%
% Frame is a class that conains data about objects from a single time point
%
% See also: Object, DataObject, TimeSeries

   properties
      t = 0;           % absolute time of frame
      objects  = [];   % segmented objects in image (e.g. array of Object, Cell, classes...)
      
      experiment = []; % reference to Experiment class
      %timeseries = []; % optional reference to the time series this frame belongs to
      source     = []; % image source
   end
   
   methods
      
      function obj = Frame(varargin)
         %
         % Frame()
         % Frame(frame)
         % Frame(...,fieldname, fieldvalue,...)
         %

         if nargin == 1 && isa(varargin{1}, 'Frame') %% copy constructor
            obj = copy(varargin{1});
         else
            obj.fromParameter(varargin);
         end
      end
      
      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, [], varargin);
      end
      

      function newobj = copy(obj)
      % 
      % f = copy(obj)
      %
      % description:
      %    deep copy of the frame and its objects
      %
         nobjs = length(obj);
         newobj(nobjs) = Frame();
         for k = 1:nobjs
            newobj(k).t          = obj(k).t;
            newobj(k).objects    = obj(k).objects.copy;
            newobj(k).experiment = obj(k).experiment; % shallo copy
            newobj(k).source     = obj(k).source;
            %newobj(k).timeseries = obj(k).timeseries;
         end 
      end
      
      
      function data = toArray(obj)
      %
      % data = toArray(obj)
      %
      % description:
      %    convert data of all objects to array
      %  
         data = obj.objects.toArray;
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % basic functionality
      
      function d = dim(obj)
         %
         % d = dim()
         %
         % spatial dimension of object's position
         %
         d = obj(1).objects(1).dim;
      end
      
      
      function n = nobjects(obj)
      %
      % n = nobjects(obj)
      %
      % output:
      %   n    number of objects in this frame
      %
         n = length(obj.objects);
      end


      function t = time(obj)
      %
      % t = time(obj)
      %
      % output:
      %   t    time of the frame
      %
         if ~isempty(obj.t)
            t = obj.t;
         elseif length(obj) == 1
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
      %   xyz    coordinates of the objects in the image as column vectors
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
      
      
      function img = data(obj, varargin)
      %
      % img = readData()
      %
      % returns the image data for this frame
      %
         img = obj.source.data(varargin);   
      end
      
      
      %TODO:
      function img = readData(obj) % only for compability with this version: todo remove
      %
      % img = readData()
      %
      % returns the image data of this frame
      %
         img = obj.getImage();   
      end
      
 

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % utils
      
      function obj = transformCoordinates(obj, R, T, C)
      %
      % obj = transformCoordinates(obj, R, T, C)
      %
      % applies rotation R, scaling C  and translation T to coordinates of objects r
      %  
         obj.objects = obj.objects.transformCoordinates(R,T,C);
         
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % visualization 

      % image
      function imglab = labeledImage(obj)
         imglab = obj.objects.labeledImage();
      end
      
      function plot(obj)
         implot(obj.data);
      end

   end
   
end
