 classdef Colony < matlab.mixin.Copyable
    
    properties 
       source     = [];   % image source pointing to the data       
       roi        = [];   % region information about the the colony (ROI class)

       objects    = [];   % object data, i.e. array of ObjectData classes representing the segmentation result
    end
    
    methods
    
      function obj = Colony(varargin)
         %
         % Colony(source)
         % Colony(source)
         % Colony(source, roi)
         % Colony(source, roi, objects)
         % Colony(..., fieldname, fieldvalue, ...)
         %

         if nargin == 0
            return
         elseif nargin >= 1
            if isa(varargin{1}, 'Colony') %% copy constructor
               obj = copy(varargin{1});
            elseif isa(varargin{1}, 'ImageSource')
               obj.isource = varargin{1};
               if nargin > 1 && isa(varargin{2}, 'ROI')
                  obj.fromImageSourceAligendAndROI(varargin{1}, varargin{2});
               end
               if nargin > 2 && isa(varargin{3}, 'Object')
                  obj.iobjects = varargin{3};
               end

            else
               obj.fromParameter(varargin);
            end
         end
      end

      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, [], varargin);
      end
      
      function obj = fromImageSourceAligendAndROI(obj, ia, roi)
         obj.source = ia.copy();
         r = roi.copy();
         %[~, r] = obj.source.reduceToROI(r);
         obj.roi = r;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % basics
      
      function img = data(obj)
         img = obj.source.dataExtract(obj.roi.boundingBox);
      end
      
      function img = extract(obj)
         img = obj.source.dataExtract(obj.roi);
      end
      
      function img = mask(obj)
         %todo: does no fit to returned data as that gets extraced !!
         img = obj.roi.mask(obj.source.dataSize);
         img = obj.roi.boundingBox.dataExtract(img);
      end
 
      function img = labeledImage(obj)
         img = obj.objects.labeledImage();
      end
      
      
%       function save(obj, filename)
%          save(filename, 'obj');
%       end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % analysis
      
      % add routines here
 
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % visualization

      function img = plotLabeledImage(obj)
         % overlay image and data
         imgo = obj.data();
         imgo = imoverlaylabel(imgo,  obj.objects.labeledImage());    
         imgo = implot(imgo);
         
         if varargout > 0
            img = imgo;
         end
      end
    
    end
    
    
 end