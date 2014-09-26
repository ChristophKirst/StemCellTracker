 classdef Colony < matlab.mixin.Copyable
    
    properties 
       isource     = [];   % image source pointing to the data       
       iroi        = [];   % region information about the the colony (ROI class)

       iobjects    = [];   % object data, i.e. array of ObjectData classes representing the segmentation result
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
               for i = 1:2:nargin % constructor from arguments
                  if ~ischar(varargin{i})
                     error('%s: invalid constructor input, expects char at position %g',class(obj), i);
                  end
                  if isprop(obj, lower(varargin{i}))
                     obj.(lower(varargin{i})) = varargin{i+1};
                  else
                     warning('%s: unknown property name: %s ', class(obj), varargin{i})
                  end
               end
            end
         end
      end

      function obj = fromImageSourceAligendAndROI(obj, ia, roi)
         obj.isource = ia.copy();
         r = roi.copy();
         [~, r] = obj.isource.reduceToROI(r);
         obj.iroi = r;
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % basic functionality
      
      function img = data(obj)
         img = obj.isource.extractdata(obj.iroi.boundingbox);
      end
      
      function dat = objects(obj)
         dat = obj.iobjects;
      end
      
%       function save2disk(obj, filename)
%          save(filename, 'obj');
%       end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % analysis
      
      
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % 

      
      function img = plotDataAndImage(obj)
         % overlay image and data
         
         imgo = obj.data();
         imgo = imoverlaylabel(imgo,  obj.iobjects.labeledImage());
         
         imgo = implot(imgo);
         
         if varargout > 0
            img = imgo;
         end
         
      end
    
    end
    
    
 end