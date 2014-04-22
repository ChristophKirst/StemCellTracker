classdef FileHandler < handle
   % FileHandler is a class that handles file access to raw image data and results
   %
   % note:
   %      image data acces is divided into two routines: 
   %             - readImage   reads the data necessary for segmentation, e.g. the full image in each Frame
   %             - loadImage   allows for reading also sub sets of images, and is intended for user insection, testing 
   %                           
   
   properties
      % file info
      BaseDirectory       = '.';      % all directories are w.r.t to this base directory
      ImageDirectoryName  = '';       % directory name in which the image data is stored relative to BaseDirectory
      ResultDirectoryName = '';       % directory for saving the results relative to BaseDirectory

      ReadImageCommandFormat = 'imload(<file>)'; % command to open the image data on which segmentation is performed, typically the image data in one Frame, eg: 'imload(<directory>\<file>)'
      ReadImageFileFormat    = '';   % format of individual file names, example: 'IMG_T<time,3>_C<channel,2>_P<position,2>.TIF'      
      ReadImageTagNames      = {};   % images are obtinaed by replacing the tags in ImageFileFormat and ReadImageCommandFormat with numbers, the ordering of these numbers stored in Frame or derived classes is specified in this cell array
      ReadImageTagRange      = {};   % cell array with ranges of the tag indices, each range is specifiede by {min, max}, {min, max, delta} or array of indices    
   end
   
   properties (Dependent)
      ImageDirectory;        % absolute image directory
      ResultDirectory;       % absolute result directory  
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
            
            obj.initialize();           
         end
      end
      
      
      function newobj = copy(obj)
         %
         % c = copy(obj)
         %
         % description:
         %     deep copy of the object
         
         nobjs = length(obj);
         newobj(nobjs) = FileHandler();
         props = properties(newobj);
         for k = 1:nobjs
            for i = 1:length(props)
               newobj(k).(props{i}) = obj(k).(props{i});
            end
         end
      end
      
      function initialize(obj)
         %
         % initialize(obj)
         %
         % description: initializes missing info and checks consistency
            
         if ~isdir(obj.BaseDirectory)
            warning('BaseDirectory %s is not a directory!', obj.BaseDirectory);
         end
         
         if ~isdir(obj.ImageDirectory)
            warning('ImageDirectory %s is not a directory!', obj.ImageDirectory);
         end
         
         % check tag names 
         cmd = obj.ReadImageCommand([]);
         tnames = tagformat2tagnames(cmd);
         
         if isempty(obj.ReadImageTagNames)
            obj.ReadImageTagNames = tnames;
            fprintf('FileHandler; generated tag names: %s', var2char(tnames));
         end
         
         extratags = setdiff(tnames, obj.ReadImageTagNames);
         if ~isempty(extratags)
            warning('There are more tags in the format than specified, adding extra tags to the list!');
            obj.ReadImageTagNames = [obj.ReadImageTagNames, extratags];
         end
         if ~isempty(obj.ReadImageTagNames, tnames)
            warning('There are more tags in ReadImageTagNames than used in the format!');   
         end
         
         % tag ranges   
         if isempty(obj.ReadImageTagRange)
            tfrmt = obj.fileName();
            fns = tagdir(obj.fileName());
            tags = name2tags(tfrmt, fns, obj.ReadImageTagNames);
            
            for i = 1:size(tags, 1)
               obj.ReadImageTagRange{i} = {min(tags(i,:)), max(tags(i,:))};
               fprintf('ReadImageTagRange; generated tagrange for %s: %s', obj.ReadImageTagNames{i},  var2char(obj.ReadImageTagRange{i}));
            end
         end
      end
    

      % directories 
      function ddir = get.ImageDirectory(obj)
         % description:
         % returns the image data directorty for the experiment
         %
         % output:
         % ddir image data directory
         
         ddir = fullfile(obj.BaseDirectory, obj.ImageDirectoryName);
      end
      
      function ddir = get.ResultDirectory(obj)
         % description:
         % returns the results directorty for the experiment
         %
         % output:
         % ddir result directory
         
         ddir = fullfile(obj.BaseDirectory, obj.ResultDirectoryName);
      end

      
      
      %%% ReadImage functionality 
 
      function fn = fileName(obj, tags)
         % description:
         %     returns the filename specified by tags
         %
         % output:
         %     fn    filename
         %
         % Experiment.readImageCommand, Experiment.readImage, tags2name, tagformat
         
         fn = obj.ReadImageTagFormat;
         
         if nargin > 1
            fn = tags2name(fn, tags, obj.ReadImageTagNames);
         end
         
         fn = fullfile(obj.ImageDirectory, fn);
      end
 
      function cmd = ReadImageCommand(obj, tags)
         % description:
         % returns the command to open an image specified by tags
         %
         % output:
         %     cmd the comand string the opens the image using eval(cmd)
         %
         % See also:  Experiment.readImage, tags2name, tagformat
         
         cmd = obj.ReadImageCommandFormat;
         cmd = strrep(cmd, '<file>', obj.fileName);
         
         if nargin < 2
            return
         else
            if isempty(obj.ReadImageTagNames) || isstruct(tags)
               cmd = tags2name(cmd, tags);
            else
               cmd = tags2name(cmd, tags, obj.ReadImageTagNames);
            end
         end
      end
      
      function data = readImage(obj, varargin)
         % description:
         %     returns image data as generated by readDataCommand
         %
         % output:
         %     data    image data
         %
         % See also: Experiment.readImageCommand, tags2name, tagformat
         
         data = eval(obj.ReadImageCommand(varargin{:}));     
      end

      function b = checkTags(obj, tags)
         %
         % b = checkTags(tags)
         %
         % description: checks tags to lie in specified tag range
         
         b = 1;
         for i = 1:length(obj.ReadImageTagNames)
            range = obj.ReadImageTagRange{i};
            if iscell(range)
               if length(range) == 2
                  if tags(i) < range{1} || tags(i) > range{2}
                     b = 0;
                     return
                  end
               else 
                  range = range{1}:range{3}:range{2};
                  if ~any(range == tags(i));
                     b = 0; 
                     return
                  end
               end
            else
               if ~any(range == tags(i));
                  b = 0;
                  return
               end
            end
         end
      end
      
      
      %%% Saving
      function fn = saveResult(obj, fname, data, varargin) %#ok<INUSL>
         %
         % fn = saveResult(obj, fname, varargin)
         %
         % description:
         % save a result in the result folder creating subdirectories if necessary
         %
         % input:
         % fname the file name relative to the ResultDirectory
         % varargin arguments to save
         %
         % output:
         % fn the file name where results were saved
         %
         % See also: loadResults
         
         fn = fullfile(obj.ResultDirectory, fname);
         
         % create subfolder if necessary
         fnpath = fileparts(fn);
         if ~isdir(fnpath)
            mkdir(fnpath)
         end
         
         %save
         save(fn, 'data', varargin{:});
      end
      
      
      function result = loadResult(obj, fname, varargin)
         %
         % result = LoadResult(fname, varargin)
         %
         % description:
         % load a result in the result folder
         %
         % input:
         % fname the file name relative to the ResultDirectory
         % varargin arguments to load
         %
         % output:
         % result the results
         %
         % See also: saveResults
         
         fn = fullfile(obj.ResultDirectory, fname);
         result = load(fn, varargin{:});
         result = result.data;
      end     
      
      %%% Info
      function txt = infoString(obj)
         %
         % txt = InfoString()
         %
         % descrition:
         % returns text with information about the experiment
         %
    
         txt = '';
         txt = [txt, 'Directories:\n'];
         txt = [txt, 'BaseDirectory : ', obj.BaseDirectory, '\n'];
         txt = [txt, 'ImageDirectory : ', obj.ImageDirectory, '\n'];
         txt = [txt, 'ResultDirectory: ', obj.ResultDirectory, '\n'];
         
         %r = obj.Result;
         %inf = whos('r');
         %txt = [txt, '\nMemory: ' num2str(inf.bytes/1024) ' kB\n'];
      end
      
      function info(obj)
         %
         % Info()
         %
         % descrition:
         % prints text with information about the experiment
         %
         
         fprintf(obj.infoString());
      end

   end % methods
   
end