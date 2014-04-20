classdef Segment < handle
%
% Segment class handling outline of objects etc
%
% See also: Object, Trajectory, Frame

   properties 
      nucleus   = [];  % nuclear pixel stored in compressed format  (coordinates are w.r.t. to image in which the segmentation was done, accesible via the class specified in frame)
      cytoplasm = [];  % cytoplasmic pixel stored in compressed format (coordinates are w.r.t. to image in which the segmentation was done, accesible via the class specified in frame)
   end
   
   %properties (Dependent)
   %   length         % length of data entries obtained when using toArray
   %end
   
   methods
      
      
      
      function pixelIdxList()
         
      end
      
      

      
      
      
      
   end % methods
   
   
   
   
end % classdef
      
      