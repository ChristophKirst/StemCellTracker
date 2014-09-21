classdef Experiment < FileHandler
%
% Experiment class for storing info, filenames and results of an experiment 
%
% See also: FileHandler, Object, Frame, TimeSeries

   properties
      name@char        = '';        % Name of experiment
      date             = '';        % Date of experiment
      description@char = '';        % Description of the Experiment / notes

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
      
      
      function obj = fromFileExpression(obj, fexpr)
         % set up experiment given the file expression
         
      end
      
      function obj = fromReadCommand(obj, rcmd)
         % set up experimant from read command
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
      
      









classdef FileHandler < matlab.mixin.Copyable
   % FileHandler is a class that organizes file acces and saving
   %
   % See also: Experiment
   
   properties
      % file info
      fbasedirectory      = '';       % all directories are w.r.t to this base directory
      fimagedirectory     = '';       % directory name in which image data is stored relative to fbasedirectory
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

             
            obj.initialize();           
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
            tfrmt = obj.ImageFile();
            fns = tagformat2files(tfrmt);
            
            if isempty(fns)
               warning('FileHandler: ReadImageFileFormat %s, does not point to any files!', obj.ReadImageFileFormat);
            else
               fprintf('FileHandler: %g files of format %s found!\n', length(fns), tfrmt);
               fprintf('FileHandler: first file is %s\n', fns{1});
            end
         end

      end
    

      % directories 
      function ddir = ImageDirectory(obj, varargin)
         % description:
         % returns the image data directorty for the experiment
         %
         % output:
         %     ddir    image data directory
         
         ddir = fullfile(obj.BaseDirectory, obj.ImageDirectoryName);
         
         if (nargin > 1) 
            ddir = tags2name(ddir, varargin{:});
         end
      end
      
      function ddir = ResultDirectory(obj, varargin)
         % description:
         % returns the results directorty for the experiment
         %
         % output:
         %     ddir    result directory
         
         ddir = fullfile(obj.BaseDirectory, obj.ResultDirectoryName);
         
         if (nargin > 1) 
            ddir = tags2name(ddir, varargin{:});
         end
      end

      
      
      %%% Image tagging functionality
      
      
      function fn = ImageFile(obj, varargin)
         %
         % fn = ImageFile(tagspec)
         %
         % description:
         %     returns the image filename specified by the tagsspec
         %
         % input: 
         %     tagspec  (optional) tag specification as struct or parameter list
         %
         % output:
         %     fn       filename
         %
         % See also: FileHandler.ReadImageCommand, FileHandler.readImage, tags2name, tagformat
         
         fn = obj.ReadImageFileFormat;
         fn = tags2name(fn, varargin{:});
         fn = fullfile(obj.ImageDirectory, fn);
      end
      
  
      function cmd = ReadImageCommand(obj, varargin)
         %
         %  cmd = ReadImageCommand(obj, varargin)
         %
         % description:
         %     returns the command to open an image
         %
         % input: 
         %     tagspec  (optional) tag specification as struct or parameter list
         %
         % output:
         %     cmd the comand string the opens the image using eval(cmd)
         %
         % See also:  FileHandler.readImage, tags2name, tagformat
         
         cmd = obj.ReadImageCommandFormat;
         cmd = strrep(cmd, '<file>', obj.ImageFile);
         cmd = tags2name(cmd, varargin{:});
      end
      
      function tags = TagNames(obj, varargin)
         % tags = TagNames(obj)
         %
         % description:
         %     returns all tag names
         %
         % output:
         %     tags    tag names
 
         tags = tagformat2tagnames(obj.ReadImageCommand);
      end
      
      

 
      
      function fns = findImageFiles(obj, varargin)
         %
         % fns = findImageFiles(tagspec)
         %
         % description:
         %     returns all filenames that match the file format 
         %
         % input: 
         %     tagspec  (optional) tag specification as struct or parameter list
         %
         % output:
         %     fns       filenames
         %
         % See also: FileHandler.readImageCommand, FileHandler.readImage, tags2name, tagformat
         
         fns = tagformat2files(obj.ImageFile(varargin{:})); 
      end
      

      function tags = findImageTags(obj, varargin) 
         %
         % fn = findTags(tagspec)
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
         % note:
         %     tags in the ReadCommand are not inferred
         %
         % See also: FileHandler.files, FileHandler.readImage, tags2name, tagformat  
         
         fn = obj.ImageFile(varargin{:});
         fns = tagformat2files(fn); 
         tags = name2tags(fn, fns);
      end
      
      
      
      function tr = findImageTagRange(obj, varargin) 
         %
         % tr = findImageTagRange(tagspec)
         %
         % description:
         %     returns the ranges of all image tags not specified in tagspec
         %
         % input: 
         %     tagspec  (optional) tag specification as struct or parameter list
         %
         % output:
         %     tr       ranges of the tags as struct: name -> [min, max] or {list of stings}
         %
         % See also: FileHandler.files, FileHandler.readImage, tags2name, tagformat  
         
         tags = obj.findImageTags(varargin{:});
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
      
      
      function fn = ResultFile(obj, fname)
         %
         % fn = ResultFile(fname)
         %
         % description:
         %     returns the absolute result filenameof fname
         %
         % input: 
         %     fname  (optional) relative result file
         %
         % output:
         %     fn       filename
         %
         % See also: FileHandler.ReadImageCommand, FileHandler.readImage, tags2name, tagformat
         
         fn = fullfile(obj.ResultDirectory, fname);
      end
      

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
         %      load data from the result folder
         %
         % input:
         %      fname the file name relative to the ResultDirectory
         %      varargin arguments to load
         %
         % output:
         %      result the results
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
         %       returns text with information about the experiment
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
            txt = [txt, 'Number of Files       : ', num2str(length(obj.findImageFiles)), '\n'];
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
         %      prints text with information about the experiment
         %
         % See also: FileHandler.infoString
         
         fprintf(obj.infoString());
      end

   end % methods
   
end



   

