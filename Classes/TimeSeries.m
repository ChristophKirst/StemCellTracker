classdef TimeSeries < handle
%
% TimeSeries class for storing basic time series data 
%
% See also: Plate, Colony, Cell, Frame, Object

   properties
      frames = [];       % individual time frames of the time series (array of Frame or Colony or any class with field 'objects' and methods as in Frame)
      trajectories = []; % trajectory data (vector of trajectories that match the data in frames)
      
      experiment = [];          % reference to Experiment class
   end
   
   methods
      function obj = TimeSeries(varargin)  % basic constructor
         %
         % TimeSeries()
         % TimeSeries(..., fieldname, fieldvalue, ...)
         %
         
         for i = 1:2:nargin
            obj.(varargin{i}) = varargin{i+1};
         end
         
      end
      
      %%% basic routines for extracting time series data using objects and trajectoris

      %function data = 
      

   end % methods
   
end