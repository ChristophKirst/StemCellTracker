 classdef ImageSourceFiles < ImageSource
   %
   % ImageSourceFiles class represents images accessed with the tag expression framework, individual files are read with imreadBF
   %
   
   properties       
      ifileexpression  = '';      % tagged file / directory
      ifiletagrange    = struct;  % tag range for file tags
   end

   methods
      function obj = ImageSourceFiles(varargin) % constructor
         %
         % ImageSourceFiles()
         % ImageSourceFiles(taggedfilename)
         % ImageSourceFiles(filenames)
         % ImageSourceFiles(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageSourceFiles') %% copy constructor
               obj = copy(varargin{1});
            elseif ischar(varargin{1})
               obj.fromFileExpression(varargin{1});
            elseif iscellstr(varargin{1})
               obj.fromFileList(varargin{1});
            else
               error('%s: invalid constructor input, expects char at position %g',class(obj), 1);
            end
         else
            obj.fromParameter(varargin);
         end
      end
      
      
      function obj = initializeRawDataFromFile(obj, filename)
         % sets raw data specs using the info in file
         info = imreadBFInfo(filename);
         
         % singleton dimensions are ignored / if more than one series -> warning
         if info.rawCellSize > 1
            warning('%s: initializeRawDataFromFile: more than a single series per file, only using first !', class(obj))
         end

         % copy the raw entries
         cp = {'irawdatasize', 'irawdataformat', 'irawdataclass'};
         for i = 1:length(cp)
            obj.(cp{i}) = info.(cp{i});
         end
      end
      
      function obj = initializeRawCellFromFileExpression(obj, fileExpr, varargin)
         % initializes the cell format and size from the file expr and additional tag ranges

         tagnames = tagExpressionToTagNames(fileExpr);
         if any(cellfun(@length, tagnames) > 1)
            error('%s: initializeRawCellFromFileExpression: tag names need to be single character!', class(obj));
         end

         obj.ifileexpression      = fileExpr;
         obj.ifiletagrange        = tagRangeFromTagExpression(fileExpr, varargin{:});
         %obj.irange               = imfrmtRangeToIndexRange(obj.ifiletagrange, obj.ifiletagrange);
         obj.initializeRangeKeyFromFileTagRange;
         
         obj.irawcellformat       = cell2mat(tagnames);
         obj.irawcellsize         = tagRangeSize(obj.ifiletagrange);

      end


      function obj = initializeRawCellFromFileList(obj, fname, varargin)
         %
         % obj = initializeRawCellFromFileList(obj, fname, varargin)
         %
         % description:
         %    infer properties from a list of image files
   
         [fileExpr, tagnames, tags] = tagExpression(fname, varargin{:}, 'tagnames', num2cell('SUVWABC'));

         if any(cellfun(@length, tagnames) > 1)
            error('%s: initializeRawCellFromFileList: tag names need to be single character!');
         end

         obj.ifileexpression      = fileExpr;
         obj.ifiletagrange        = tagRangeFromTags(tags);
         %obj.irange               = imfrmtRangeToIndexRange(obj.ifiletagrange, obj.ifiletagrange);
         obj.initializeRangeKeyFromFileTagRange;
         
         obj.irawcellformat       = cell2mat(tagnames);
         obj.irawcellsize         = tagRangeSize(obj.ifiletagrange);

      end
               
      function obj = initializeRangeKeyFromFileTagRange(obj)
         
         obj.irangekey = obj.ifiletagrange;
         
         if ~isempty(obj.irangekey)          
            fnames = fieldnames(obj.irangekey);
            rmnames = {};
            for i = 1:length(fnames)
               v = obj.irangekey.(fnames{i});
               if ischar(v)
                  v = {v};
                  obj.irangekey.(fnames{i}) = v;
               end
               if ~iscellstr(v)
                  if iscell(v)
                     v = cell2mat(v);
                     if isequal(sort(v), 1:length(v))
                        rmnames = [rmnames, fnames{i}]; %#ok<AGROW>
                     end
                  end
               end
            end
            
            obj.irangekey = rmfield(obj.ifiletagrange, rmnames);
         end
      end
      


      function obj = fromFileExpression(obj, fileExpr, varargin)
         obj.initializeRawCellFromFileExpression(fileExpr, varargin);
         firstfile = obj.rawFileName(1);
         obj.initializeRawDataFromFile(firstfile);
         
         obj.initializeDataAndCellFromRaw;
      end

      function obj = fromFileList(obj, fileList, varargin)
         obj.initializeRawCellFromFileList(fileList, varargin);
         firstfile = obj.rawFileName(1);
         obj.initializeRawDataFromFile(firstfile);

         obj.initializeDataAndCellFromRaw;
      end


 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% files
      
      function range = fileTagRange(obj, varargin)
         if nargin == 1
            range = obj.ifiletagrange; 
         else
            range = obj.rawRange(varargin{:}); % returns index range
            range = imfrmtRangeFromIndexRange(obj.ifiletagrange, range); % convert to filetag range
            range = imfrmtReformatRange(obj.rawCellDataSize, obj.rawCellDataFormat, obj.rawCellDataFormat, range);
            range = imfrmtParseRangeLast(obj.ifiletagrange, range); % complement all missing ranges
          end
      end

      function tagRange = rawFileTagRangeFromVarargin(obj, varargin)
         % parses raw tag range
         
         if nargin > 1 && isnumeric(varargin{1})
            if nargin > 2
               id = sub2ind(obj.rawCellSize, varargin{:});
            else
               id = varargin{1};
            end
            tagRange = tagFromIndex(obj.ifiletagrange, id);
         else
            tagRange = imfrmtParseRangeLast(obj.ifiletagrange, varargin);
         end
         
         %tagRange = imfrmtRangeFromIndexRange(obj.ifiletagrange, tagRange);
      end

      
      function fsi = fileTagSize(obj, varargin)
         range = obj.fileTagRange(varargin{:});
         fsi = tagRangeSize(range);
      end

      function n = nFiles(obj, varargin)
           n = prod(obj.fileTagSize(varargin{:}));
      end


      function fe = fileExpression(obj)
         fe = obj.ifileexpression;
      end

      function fn = fileName(obj, varargin)
         % file expression for the specified ranges

         % determine range
         range = obj.fileTagRange(varargin{:});
         tags  = tagRangeToTags(range);

         fn = tagExpressionToString(obj.ifileexpression, tags);
         %fn = reshape(fn, imfrmtRangeSize(obj.rawCellSize, obj.rawCellFormat, rawRange);
      end
      
      function [fn, tags] = rawFileName(obj, varargin)
         % file expression for the specified raw tag ranges
             
         % determine range
         range = obj.rawFileTagRangeFromVarargin(varargin{:});
         tags = tagRangeToTags(range);
         
         fn = tagExpressionToString(obj.ifileexpression, tags);
      end
      
  
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      %%% data 

      function d = getRawData(obj, varargin)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw image data in the format given by obj.rawformat
         %     data specs use raw formats and sizes

         % range
         range = imfrmtRangeToIndexRange(obj.irangekey, imfrmtParseRangeLast(varargin{:}));
         
         % check if all cell dims are singeltons
         cid = obj.rawCellIndex(range);
         if length(cid) > 1
            error('%s: getRawData: cell dimensions are not specified to be singeltons!', class(obj));
         end
        
         % get raw file names
         [fn, tag] = obj.rawFileName(range);
         
         tag = imfrmtRemoveRange(tag, obj.rawCellFormat);

         d = imreadBF(fn, tag, 'squeeze', false, 'metadata', false);
         d = imfrmtReformat(d, 'XYZCT', obj.rawDataFormat);   
      end
 
      function d = getRawCellData(obj, varargin)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw cell data in the format given by obj.rawformat
         %     data specs use raw formats and sizes
         
         % range         
         indexRange = imfrmtParseRangeLast(varargin{:});
         range = imfrmtRangeFromIndexRange(obj.irangekey, indexRange);

         % get raw file names
         [fn, tags] = obj.rawFileName(range);
         if ischar(fn)
            fn = {fn};
         end

         % allocate raw data mem
         d = cell(obj.rawCellSize(range));
         
         tags = imfrmtRemoveRange(tags, obj.rawCellFormat);
         for i = 1:length(fn)
            d{i} = imreadBF(fn{i}, tags(i), 'squeeze', false, 'metadata', false);
            %size(d{i})
            d{i} = imfrmtReformat(d{i}, 'XYZCT', obj.rawDataFormat);
            %size(d{i})
         end
      end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% info
      
      function info = infoString(obj)
         info = infoString@ImageSource(obj, 'Files');
         info = [info, '\nexpression: ', obj.fileExpression];
         info = [info, '\nfiles     : ', var2char(obj.nFiles)];
      end

   end
   
   
end