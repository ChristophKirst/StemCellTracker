classdef Cell < Object
%
% Cell class extending Object class in storing additional statistics / analysis results
%
% See also: Object, Trajectory, Frame
   properties
      nucleus = [];        % nuclear pixel stored in compressed format (coordinates are w.r.t. to image in which the segmentation was done, accesible via the class specified in frame)
      cytoplasm = [];      % cytoplasmic pixel stored in compressed format (coordinates are w.r.t. to image in which the segmentation was done, accesible via the class specified in frame)
      
      statistics = [];     % any statistics for the cell, array of values corresponding to statistics as in the StatisticsNames list of the associated experiment class, column vector to use [Cell.statistics]
      
      experiment = [];     % reference to Experiment class
      frame = [];          % frame to which this Cell belongs too
   end
   
   %properties (Dependent)
   %   length         % length of data entries obtained when using toArray
   %end
   
   methods

      
      
      
      
   end % methods
   
   
   
   
end % classdef
      
      