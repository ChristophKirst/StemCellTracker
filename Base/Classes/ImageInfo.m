classdef ImageInfo < matlab.mixin.Copyable
   
   % class that manages information about image, including basic measure such as size, format, coloring and 
   
   properties
      isize   = [];               % size of the image
      isizePQLCT = [0,0,0,0,0];   % size of the image for full pqlct dimensions
      
      iformat = '';               % format of the image when obtained with the getRawData routine
      
      iclass = '';                % class of the image
      
      icolor = {};                % colors for the channels
      
      iseries = [];               % (optional) series number if image is from a larger file 
      inimages = [];              % (optional) number of total 2d images

      imetadata = [];             % (optional) meta data
      
      iscale = [];                % (optional) spatial scale of image in pixel per spatial unit
      iunit = '';                 % (optional) spatial unit
      
      iname  = [];                % (optional) name of the image 
      
      ilabel = [];                % translation between lable (e.g. 'dapi') and a subset of the data (e.g. 'c' -> 1)
   end
   
   methods
      function p = sizeP(obj)
         p = obj.isizePQLCT(1);
      end
      function q = sizeQ(obj)
         q = obj.isizePQLCT(2);
      end
      function l = sizeL(obj)
         l = obj.isizePQLCT(3);
      end
      function c = sizeC(obj)
         c = obj.isizePQLCT(4);
      end
      function t = sizeT(obj)
         t = obj.isizePQLCT(5);
      end

      
      function obj = pqlctsizeFromFormatAndSize(obj)
         frmt = 'pqlct';  
         pos = arrayfun(@(x) find(x == frmt,1), obj.iformat);
         obj.isizePQLCT = ones(1, 5);
         obj.isizePQLCT(pos) = obj.isize;
      end
      
   end
end