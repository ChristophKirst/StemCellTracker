classdef Experiment < FileHandler
%
% Experiment class for storing info, filenames and results of an experiment 
%
% See also: FileHandler, Object, Frame, TimeSeries

   properties
      name@char        = '';        % Name of experiment
      date             = '';        % Date of experiment
      description@char = '';        % Description of the Experiment / notes
 
      result  = [];                 % reference to result data of segmentation / analysis
 
      microscope      = '';         % name of microscope with which images are taken
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
         
         obj.date = datestr(now);
         
         for i = 1:2:nargin
            if ~ischar(varargin{i})
               error('Experiment: invalid constructor input, expects char at position %g', i);
            end
            if isprop(obj, varargin{i})
               obj.(varargin{i}) = varargin{i+1};
            else
               warning('Experiment: unknown property: %s', varargin{i})
            end        
         end
         
         obj.initialize();
         
      end
      
      
      function fname = saveExperiment(obj, fname)
         if nargin < 2
            fname = obj.name;
         end
         if isempty(fname)
            fname = ['experiment_' datestr(now, 'yyyy-mm-dd') '.mat'];
         end

         fname = fullfile(obj.ResultDirectory, fname);
         data = obj; %#ok<NASGU>
         
         save(fname, 'data'); 
      end
      
      %%% info
      
      function txt = infoString(obj)
         % 
         % txt = infoString()
         %
         % descrition:
         %    returns text with information about the experiment
         %
         
         txt = ['\nExperiment : ' obj.name '\n'];
         if ~isempty(obj.date)
            txt = [txt, 'Date       : ', obj.date, '\n'];
         end
         
         if ~isempty(obj.description)
            txt = [txt, 'Description: ', obj.description, '\n'];
         end
         
         if ~isempty(obj.microscope)
            txt = [txt, 'Microscope : ', obj.microscope , '\n'];
         end
         
         txt = [txt, obj.infoString@FileHandler];
                 
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
      
      



   

