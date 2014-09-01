classdef TimeSeries < handle
%
% TimeSeries class for storing basic time series data 
%
% See also: Experiment, Frame, Object

   properties
      frames = [];       % individual time frames of the time series (array of Frame class or any class with field 'objects' and methods as in Frame)
      trajectories = []; % trajectory data (vector of trajectories that match the data in frames)
      
      experiment = [];   % reference to Experiment class
   end
   
   methods
           
      function obj = TimeSeries(varargin)
         %
         % TimeSeries()
         % TimeSeries(timeseries)
         % TimeSeries(...,fieldname, fieldvalue,...)
         %

         if nargin == 1 && isa(varargin{1}, 'TimeSeries') %% copy constructor
            obj = copy(varargin{1});
         else
            for i = 1:2:nargin % constructor from arguments
               if ~ischar(varargin{i})
                  error('%s: invalid constructor input, expects char at position %g',class(obj), i);
               end
               if isprop(obj, lower(varargin{i}))
                  obj.(lower(varargin{i})) = varargin{i+1};
               else
                  warning('%s: unknown property name: %s ', class(obj), lower(varargin{i}))
               end
            end
         end
      end

      function newobj = copy(obj)
      % 
      % f = copy(obj)
      %
      % description:
      %    deep copy of the frame and its objects
      %
         newobj = TimeSeries();
         newobj.frames = obj.frames.copy;
         newobj.trajectories = obj.trajectories.copy;
         newobj.experiment = obj.experiment; %shallow copy
      end

   end % methods
   
end