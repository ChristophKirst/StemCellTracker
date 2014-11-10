classdef FileHandler < matlab.mixin.Copyable
   % FileHandler is a class that organizes file acces and saving
   %
   % See also: Experiment
   
   properties
      % file info
      fbasedirectory      = '';       % all directories are w.r.t to this base directory
      fdatadirectory      = '';       % directory name in which image data is stored relative to fbasedirectory
      fresultdirectory    = '';       % directory for saving results relative to fbasedirectory
   end
 
   methods
      
      function obj = FileHandler(varargin)  % simple constructor
         %
         % FileHandler()
         % FileHandler(...,fieldname, fieldvalue,...)
         %
         if nargin == 1 && isa(varargin{1}, 'FileHandler') %% copy constructor
            obj = copy(varargin{1});
         else
            obj = obj.fromParameter(varargin);
         end
      end
      
      function initialize(obj)
         %
         % initialize(obj)
         %
         % description: initializes the FileHandler
            
         if ~isdir(obj.fbasedirectory)
            warning('FileHandler: basedirectory %s is not a directory!', obj.fbasedirectory);
         end
         
         if ~isdir(obj.fdatadirectory)
            warning('FileHandler: datadirectory %s is not a directory!', obj.fdatadirectory);
         end
      end
      
      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, 'f', varargin);
         obj.initialize();
      end
    
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % directories 
      
      function bdir = baseDirectory(obj)
         bdir = obj.fbasedirectory;
      end
      
      
      function idir = dataDirectory(obj, varargin)
         % description:
         % returns the image data directorty for the experiment
         %
         % output:
         %     ddir    image data directory
         
         idir = fullfile(obj.fbasedirectory, obj.fdatadirectory);
         idir = tagExpressionToString(idir, varargin);       
      end
      
      function ddir = resultDirectory(obj, varargin)
         % description:
         % returns the results directorty for the experiment
         %
         % output:
         %     ddir    result directory
         
         ddir = fullfile(obj.fbasedirectory, obj.fresultdirectory);
         ddir = tagExpressionToString(ddir, varargin);
      end
      
      function ddir = dirBaseDirectory(obj, varargin)
         ddir = dir(obj.baseDirectory(varargin));    
      end
      
      function ddir = dirDataDirectory(obj, varargin)
         ddir = dir(obj.dataDirectory(varargin));    
      end

      function ddir = dirResultDirectory(obj, varargin)
         ddir = dir(obj.resultDirectory(varargin));    
      end
      
      function fn = dataFile(obj, fname, varargin)
         fn = tagExpressionToString(fullfile(obj.dataDirectory, fname), varargin);
      end
      
      function fn = resultFile(obj, fname, varargin)
         fn = tagExpressionToString(fullfile(obj.resultDirectory, fname), varargin);
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % io 

      function fn = saveData(obj, fname, data, varargin) %#ok<INUSL>
         %
         % fn = saveData(obj, fname, varargin)
         %
         % description:
         %     save data in the result folder creating subdirectories if necessary
         %
         % input:
         %     fname    the file name relative to the ResultDirectory
         %
         % output:
         %     fn       file name where results were saved
         %
         % note: 
         %     only saves one data object at a time
         %
         % See also: FileHandler.loadData
         
         fn = fullfile(obj.resultDirectory(varargin), fname);
         
         % create subfolder if necessary
         fnpath = fileparts(fn);
         if ~isdir(fnpath)
            mkdir(fnpath)
         end
         
         %save
         save(fn, 'data', varargin{:});
      end
      
      
      function result = loadData(obj, fname, varargin)
         %
         % result = loadData(fname, varargin)
         %
         % description:
         %      load data from the result folder
         %
         % input:
         %      fname the file name relative to the ResultDirectory
         %
         % output:
         %      result the results
         %
         % See also: FileHandler.saveData
         
         fn = fullfile(obj.resultDirectory(varargin), fname);
         result = load(fn, varargin{:});
         result = result.data;
      end     
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % info 
      
      function txt = infoString(obj)
         %
         % txt = infoString()
         %
         % descrition:
         %       returns text with information about the experiment
         %
         % See also: FileHandler.info
    
         txt = '';
         txt = [txt, 'Directories:\n'];
         txt = [txt, 'Base Directory  : ', obj.baseDirectory, '\n'];
         txt = [txt, 'Data Directory  : ', obj.dataDirectory, '\n'];
         txt = [txt, 'Result Directory: ', obj.resultDirectory, '\n'];  
      end
      
      function printInfo(obj)
         %
         % descrition:
         %      prints text with information about the experiment
         %
         % See also: FileHandler.infoString
         
         fprintf(obj.infoString());
      end

   end % methods
   
end