classdef FileHandler < handle
   % FileHandler is a class that handles file access to raw image data and results
   %
   % note:
   %      image data acces is divided into two routines: 
   %             - readImage   reads the data necessary for segmentation, e.g. the full image in each Frame
   %             - loadImage   allows for reading also sub sets of images, and is intended for user insection, testing 
   %
   % See also: Experiment
   
   properties
      % file info
      BaseDirectory       = '.';      % all directories are w.r.t to this base directory
      ImageDirectoryName  = '';       % directory name in which the image data is stored relative to BaseDirectory
      ResultDirectoryName = '';       % directory for saving the results relative to BaseDirectory

      ReadImageCommandFormat = 'imload(<file>)'; % command to open the image data on which segmentation is performed, typically the image data in one Frame, eg: 'imload(<directory>\<file>)'
      ReadImageFileFormat    = '';   % format of individual file names, example: 'IMG_T<time,3>_C<channel,2>_P<position,2>.TIF'      
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
               if isprop(obj, varargin{i})
                  obj.(varargin{i}) = varargin{i+1};
               else
                  warning('%s: unknown property name: %s ', class(obj), varargin{i})
               end
            end
     
            %obj.ReadImageCommandFormat 
            %obj.ReadImageFileFormat
            %obj.ReadImageTagNames
            %obj.ReadImageTagRange  
             
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
            warning('FileHandler: BaseDirectory %s is not a directory!', obj.BaseDirectory);
         end
         
         if ~isdir(obj.ImageDirectory)
            warning('FileHandler: ImageDirectory %s is not a directory!', obj.ImageDirectory);
         end
         
         % check if files exists specified by the 
         if ~isempty(obj.ReadImageFileFormat)
            tfrmt = obj.fileName();
            fns = tagformat2files(tfrmt);
            
            if isempty(fns)
               warning('FileHandler: ReadImageFileFormat %s, does not point to any files!', obj.ReadImageFileFormat);
            else
               fprintf('FileHandler: %g, files of format %s found!\n', length(fns), tfrmt);
               fprintf('FileHandler: first file is %s\n', fns{1});
            end
         end

      end
    

      % directories 
      function ddir = get.ImageDirectory(obj)
         % description:
         % returns the image data directorty for the experiment
         %
         % output:
         %     ddir    image data directory
         
         ddir = fullfile(obj.BaseDirectory, obj.ImageDirectoryName);
      end
      
      function ddir = get.ResultDirectory(obj)
         % description:
         % returns the results directorty for the experiment
         %
         % output:
         %     ddir    result directory
         
         ddir = fullfile(obj.BaseDirectory, obj.ResultDirectoryName);
      end

      
      
      %%% ReadImage functionality 
 
      function fn = fileName(obj, varargin)
         %
         % fn = fileName(tagspec)
         %
         % description:
         %     returns the filename specified by the tags
         %
         % input: 
         %     tagspec  (optional) tag specification as struct or parameter list
         %
         % output:
         %     fn       filename
         %
         % See also: FileHandler.readImageCommand, FileHandler.readImage, tags2name, tagformat
         
         fn = obj.ReadImageFileFormat;
         fn = tags2name(fn, varargin{:});
         fn = fullfile(obj.ImageDirectory, fn);
      end
      
      function fns = files(obj, varargin)
         %
         % fn = files(tagspec)
         %
         % description:
         %     returns the files specified by the tags
         %
         % input: 
         %     tagspec  (optional) tag specification as struct or parameter list
         %
         % output:
         %     fn    filename
         %
         % See also: FileHandler.fileTags, FileHandler.readImage, tags2name, tagformat    
         
         fns = tagformat2files(obj.fileName(varargin{:})); 
      end
      
      
      function tags = fileTags(obj, varargin) 
         %
         % fn = fileTags(tagspec)
         %
         % description:
         %     returns the list of all tags not specified in tagspec
         %
         % input: 
         %     tagspec  (optional) tag specification as struct or parameter list
         %
         % output:
         %     fn    filename
         %
         % See also: FileHandler.files, FileHandler.readImage, tags2name, tagformat  
         
         fn = obj.fileName(varargin{:});
         fns = tagformat2files(fn); 
         tags = name2tags(fn, fns);
      end
      
      
      
      function tr = fileTagRange(obj, varargin) 
         %
         % fn = fileTags(tagspec)
         %
         % description:
         %     returns the list of all tags not specified in tagspec
         %
         % input: 
         %     tagspec  (optional) tag specification as struct or parameter list
         %
         % output:
         %     tr       ranges of the tags as struct: name -> [min, max] or {list of stings}
         %
         % See also: FileHandler.files, FileHandler.readImage, tags2name, tagformat  
         
         tags = obj.fileTags(varargin{:});
         tnames = fieldnames(tags);
         tr = tags(1);
         for i = 1:length(tnames)
            if isnumeric([tags.(tnames{i})])
               tr.(tnames{i}) = unique([tags.(tnames{i})]);
            else
               tr.(tnames{i}) = unique({tags.(tnames{i})});
            end
         end
      end

 
      function cmd = ReadImageCommand(obj, varargin)
         %
         %  cmd = ReadImageCommand(obj, varargin)
         %
         % description:
         %     returns the command to open an image specified by tags
         %
         % input: 
         %     tagspec  (optional) tag specification as struct or parameter list
         %
         % output:
         %     cmd the comand string the opens the image using eval(cmd)
         %
         % See also:  FileHandler.readImage, tags2name, tagformat
         
         cmd = obj.ReadImageCommandFormat;
         cmd = strrep(cmd, '<file>', obj.fileName);
         cmd = tags2name(cmd, varargin{:});

      end
      

      function tnames = tagNames(obj)
         %
         % tnames = tagNames(obj)
         %
         % description:
         %     returns a list of all tag names
         %
         % See also:  FileHandler.readImage, tags2name, tagformat
         
         
         cmd = obj.ReadImageCommand([]);
         tnames = tagformat2tagnames(cmd);
      end
      
 
      function img = readImage(obj, varargin)
         %
         % img = readImage(obj, varargin)
         %
         % description:
         %     returns image data as generated by eval(readDataCommand)
         %
         % output:
         %     img     image data
         %
         % See also: FileHandler.readImageCommand, tags2name, tagformat
         
         img = eval(obj.ReadImageCommand(varargin{:}));     
      end

      
      %%% Saving
      function fn = saveData(obj, fname, data, varargin) %#ok<INUSL>
         %
         % fn = saveData(obj, fname, varargin)
         %
         % description:
         %     save data in the result folder creating subdirectories if necessary
         %
         % input:
         %     fname    the file name relative to the ResultDirectory
         %     varargin optional arguments to save
         %
         % output:
         %     fn       file name where results were saved
         %
         % note: 
         %     only saves one data object at a time
         %
         % See also: FileHandler.loadData
         
         fn = fullfile(obj.ResultDirectory, fname);
         
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
         % load data from the result folder
         %
         % input:
         % fname the file name relative to the ResultDirectory
         % varargin arguments to load
         %
         % output:
         % result the results
         %
         % See also: FileHandler.saveData
         
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
         % See also: FileHandler.info
    
         txt = '';
         txt = [txt, 'Directories:\n'];
         txt = [txt, 'BaseDirectory  : ', obj.BaseDirectory, '\n'];
         txt = [txt, 'ImageDirectory : ', obj.ImageDirectory, '\n'];
         txt = [txt, 'ResultDirectory: ', obj.ResultDirectory, '\n'];  
         
         if ~isempty(obj.ReadImageCommandFormat)
            txt = [txt, 'ReadImageCommandFormat: ', obj.ReadImageCommandFormat, '\n'];  
         end
         if ~isempty(obj.ReadImageFileFormat')
            txt = [txt, 'ReadImageFileFormat   : ', obj.ReadImageFileFormat, '\n'];    
            txt = [txt, 'Number of Files       : ', num2str(length(obj.files)), '\n'];
         end
         
         
         
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
         % See also: FileHandler.infoString
         
         fprintf(obj.infoString());
      end

   end % methods
   
end