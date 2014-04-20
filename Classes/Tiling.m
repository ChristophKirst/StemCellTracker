classdef Tiling < handler
%
% Tiling is class to handle, load, stitch and organize tiled images
%
   properties
      images = [];         % image indices as column vectors specifying the image directly, or image names
      isize  = [];         % spatial size of each image
      
      shifts = [];         % shifts of the images to obtain the tiling, each shift is a row vector in 2d or 3d

      experiment = [];     % reference to experiment class (would only need FileHandler reference ?)
   end
     
   properties (Dependent)
      tdim = size(shifts, 2);
      idim = numel(isize);
      nimages = 
   end
   
   methods
      function obj = Tiling(varargin) 
         % constructor
         for i = 1:2:nargin
            obj.(varargin{i}) = varargin{i+1};
         end
      end
      
      
      function n = nimages(obj)
         % number of images
         n = size(obj.images,2);
      end
      
      function imgs = readImages(obj)
         % return cell array of images in tiling
         imgs = obj.experiment.ReadImages(obj.images);
      end
      
      function img = stichedImage(obj, varargin)
         % returns the stiched image
         img = obj.ReadImages();
         img = stichImages(img, obj.shifts, varargin{:});
      end  
   end
   
end