classdef Image < matlab.mixin.Copyable
   %
   % Image class
   %
   % description:
   %     class representing image data together with information about format, colors etc..
   %
   % See also: ImageSource

   properties
      data = [];       % image data if loaded
      source = [];     % image source (ImageSource class)
      format = '';     % permutation of a subset of 'xyzcta'

      scale = [];      % spatial scales
      
      size = [];       % image size
       
      colors = {};     % coloring of the individual channels 'rgb' 'gray' for RGB and gray scale images
      
      tags = [];       % special labeling (TagMap)
         
      lazy = true;     % if true postpone loading of the data until needed if image is given by a ImageSource
      caching = true;  % if keep loaded data in memory, other wise clear data
      
      info = [];       % meta information
   end
   
   
   
   methods
      function obj = Image(varargin)  % basic constructor
         %
         % Image()
         % Image(...,fieldname, fieldvalue,...)
         %
         if nargin == 1 
            if isa(varargin{1}, 'Image') %% copy constructor
               obj = copy(varargin{1});
            elseif isnumeric(varargin{1})
               obj.data = varargin{1};
            elseif isa(varargin{1}, 'ImageSource')
               obj.source = varargin{1};
            else
               error('%s: not valid arguments for constructor', class(obj));
            end
         else
            for i = 1:2:nargin % constructor from arguments
               if ~ischar(varargin{i})
                  error('%s: invalid constructor input, expects char at position %g',class(obj), i);
               end
               if isprop(obj, varargin{i})
                  obj.(varargin{i}) = varargin{i+1};
               else
                  warning('%s: unknown property name: %s ', class(obj), varargin{i})
               end
            end
   
            obj.initialize();           
         end
      end
      
      
      function intialize(obj)
         
         
         
         if isempty(obj.format)
            obj.format = imformat(obj.data)
         
         
      end
      
      
      
      function 
      
      
      
      function plot
         
         
      end
      
      
      
   end
   
   
end