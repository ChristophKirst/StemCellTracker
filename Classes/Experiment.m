classdef Experiment < FileHandler
%
% Experiment class for storing info, filenames and results of an experiment 
%
% See also: FileHandler, Object, Frame, TimeSeries, Tiling

   properties
      name@char        = '';        % Name of experiment
      date             = '';        % Date of experiment
      description@char = '';        % Description of the Experiment / notes
 
      result  = [];                 % reference to result data of segmentation / analysis
 
      microscope@char = [];         % name of microscope with which images are taken
      resolution      = [];         % resolution of the images in um / Pixel in the different spatial image dimensions
      metadata        = [];         % space to store some meta data
   end
      
   methods
      
      %%% Constructor
      
      function obj = Experiment(varargin)  % basic constructor
         %
         % Experiment()
         % Experiment(..., fieldname, fieldvalue, ...)
         %
         
         for i = 1:2:nargin
            if ~ischar(varargin{i})
               error('Experiment: invalid constructor input, expects char at position %g', i);
            end
            if isprop(obj, varargin{i})
               obj.(varargin{i}) = varargin{i+1};
            else
               warnin('Experiment: unknown property name: %s', varargin{i})
            end        
         end
         
      end
        
           
      %%% info
      
      function txt = infoString(obj)
         % 
         % txt = infoString()
         %
         % descrition:
         %    returns text with information about the experiment
         %
         
         
         txt = ['Experiment : ' obj.Name '\n'];
         if ~isempty(obj.Date)
            txt = [txt, 'Date       : ', obj.Date, '\n'];
         end
         
         if ~isempty(obj.Description)
            txt = [txt, 'Description: ', obj.Description, '\n'];
         end
         
         if ~isempty(obj.Microscope)
            txt = [txt, 'Microscope : ', obj.Microscope , '\n'];
         end
         
         txt = [txt, infoString@FileHandler];
                 
         %r = obj.Result;
         %inf = whos('r');
         %txt = [txt, '\nMemory: ' num2str(inf.bytes/1024) ' kB\n'];
      end

      function info(obj)
         % 
         % Info()
         %
         % descrition:
         %    prints text with information about the experiment
         %

         fprintf(obj.infoString());
        
      end
      
   end % methods
end % classdef
      
      



   

