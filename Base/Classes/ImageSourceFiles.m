 classdef ImageSourceFiles < ImageSource
   %
   % ImageSourceFiles class represents images accessed with the tag expression framework, individual files are read with imreadBF
   %
   % note: ranges are assumed to be index ranges / range keys
   %       to access via file tags use XXXFromTagRange routines
   
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
         obj.initializeKeyFromFileTagRange;
         
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
         obj.initializeKeyFromFileTagRange;
         
         obj.irawcellformat       = cell2mat(tagnames);
         obj.irawcellsize         = tagRangeSize(obj.ifiletagrange);

      end
               
      function obj = initializeKeyFromFileTagRange(obj)
         
         obj.ikey = struct;
     
         ftr = obj.ifiletagrange;
         if ~isempty(ftr)          
            fnames = fieldnames(ftr);
%           rmnames = {};
            for i = 1:length(fnames)
               v = ftr.(fnames{i});
               if ischar(v)
                  v = {v};
                  obj.irangekey.(fnames{i}) = v;
               elseif iscellstr(v)
                  obj.irangekey.(fnames{i}) = v;
               else
                  if ~iscell(v)
                     v = num2cell(v);
                  end
                  if ~isequal(sort(cell2mat(v)), 1:length(v)) % in this way accessing file fields via 'xxx' is possible
                     obj.ikey.(fnames{i}) = cellfunc(@num2str, v);
                  end
               end
            end
         end
      end
      


      function obj = fromFileExpression(obj, fileExpr, varargin)
         obj.initializeRawCellFromFileExpression(fileExpr, varargin);
         firstfile = obj.fileNameFromRawRange(1);
         obj.initializeRawDataFromFile(firstfile);
         
         obj.setCellDataSizeAndFormatFromRaw;
      end

      function obj = fromFileList(obj, fileList, varargin)
         obj.initializeRawCellFromFileList(fileList, varargin);
         firstfile = obj.fileNameFromRawRange(1);
         obj.initializeRawDataFromFile(firstfile);

         obj.setCellDataSizeAndFormatFromRaw;
      end


 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% files
      
      function range = fileTagRange(obj, varargin)
         % returns file tags using index range in varargin
         
         range = obj.rawRange(varargin{:});
         range = imfrmtRangeFromIndexRange(obj.ifiletagrange, range); % convert to filetag range
         range = imfrmtRangeFromVarargin(obj.ifiletagrange, range); % complement all missing ranges
      end
      
      function range = fileTagRangeFromRawRange(obj, varargin)
         % varargin is raw range
         range = obj.rawRangeFromRawVarargin(varargin{:});
         range = imfrmtRangeFromIndexRange(obj.ifiletagrange, range); % convert to filetag range
         range = imfrmtRangeFromVarargin(obj.ifiletagrange, range); % complement all missing ranges
      end
      
      function range = fileTagRangeFromTagRange(obj, varargin)
         % varargin is file rag range 
         range = imfrmtRangeFromVarargin(obj.ifiletagrange, varargin); % complement all missing ranges
      end

      
      function range = rangeFromFileTagRange(obj, varargin)
         range = obj.rawRangeFromFileTagRange(varargin{:});
         range = obj.rangeFromRawRange(range);
      end

      function range = rawRangeFromFileTagRange(obj, varargin)
         range = obj.fileTagRangeFromTagRange(varargin{:});
         range = imfrmtRangeToIndexRange(obj.ifiletagrange, range);
      end

      
      function fsi = fileTagSize(obj, varargin)
         range = obj.fileTagRange(varargin{:});
         fsi = imfrmtRangeSize(obj.rawCellSize, obj.rawCellFormat, range);
      end

      function n = nFiles(obj, varargin)
           n = prod(obj.fileTagSize(varargin{:}));
      end

      function fe = fileExpression(obj, varargin)
         fe = obj.ifileexpression;
         fe = tagExpressionToString(fe, varargin);
      end

      function fn = fileName(obj, varargin)
         % file expression for the specified ranges

         % determine range
         range = obj.fileTagRange(varargin{:});
         tags  = tagRangeToTags(range);

         fn = tagExpressionToString(obj.ifileexpression, tags);
         %fn = reshape(fn, imfrmtRangeSize(obj.rawCellSize, obj.rawCellFormat, rawRange);
      end
      
      function [fn, tags] = fileNameFromRawRange(obj, varargin)  
         % determine range
         range = obj.fileTagRangeFromRawRange(varargin{:});
         tags = tagRangeToTags(range);
         
         fn = tagExpressionToString(obj.ifileexpression, tags);
      end
      
      function [fn, tags] = fileNameFromTagRange(obj, varargin)  
         % determine range
         range = obj.fileTagRangeFromTagRange(varargin{:});
         tags = tagRangeToTags(range);
         
         fn = tagExpressionToString(obj.ifileexpression, tags);
      end
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
      %%% data 

      function d = getRawData(obj, range)
         %
         % d = getRawData(obj, range)
         %
         % description:
         %     obtain raw image data in the format given by obj.rawformat
         %     range is index range in raw formats and sizes  
       
         % check if all cell dims are singeltons
         cid = obj.rawCellIndex(range);
         if length(cid) > 1
            error('%s: getRawData: cell dimensions are not specified to be singeltons!', class(obj));
         end
         
         if obj.irawcache
            if isempty(obj.irawcelldata)
               obj.irawceldata = cell(obj.rawCellSize);
            end
            
            if isempty(obj.irawcelldata{cid}) 
               % read full data from file

               % raw file name
               fn = obj.fileNameFromRawRange(range);

               d = imreadBF(fn, 'S', 1, 'squeeze', false, 'metadata', false); % there should be no series structure -> read first series
               d = imfrmtReformat(d, 'XYZCT', obj.rawDataFormat);   
               
               % cach it
               obj.irawcelldata{cid} = d;
            else
               d = obj.irawcelldata{cid};
            end
               
            % reduce to range
            if ~isemptystruct(range)
               d = imfrmtDataSubset(d, obj.rawDataFormat, range);
            end
         else
            % read directly

            % raw file name
            [fn, tag] = obj.fileNameFromRawRange(range);
            
            tag = imfrmtRemoveRange(tag, obj.rawCellFormat);
            
            d = imreadBF(fn, 'S', 1, tag, 'squeeze', false, 'metadata', false); % S = 1 no series structure 1
            d = imfrmtReformat(d, 'XYZCT', obj.rawDataFormat);
         end
      end
 
      function d = getRawCellData(obj, range)
         %
         % d = getRawData(obj, varargin)
         %
         % description:
         %     obtain raw cell data in the format given by obj.rawformat
         %     data specs use raw formats and sizes
              
         if obj.irawcache            
            if isempty(obj.irawcelldata)
               obj.irawcelldata = cell(obj.rawCellSize);
            end

            cid = obj.rawCellIndex(range);

            cidload = cellfun(@isempty, obj.irawcelldata(cid));
            cidload = cid(cidload);

            if ~isempty(cidload)
               
               % get raw file names
               [fn, tags] = obj.fileNameFromRawRange(cidload);
               if ischar(fn)
                  fn = {fn};
               end
               tags = imfrmtRemoveRange(tags, obj.rawCellFormat);

               % read missing data from file
               d = cell(1, length(cidload));
               for i = 1:length(fn)
                  d{i} = imreadBF(fn{i}, tags(i), 'squeeze', false, 'metadata', false);
                  d{i} = imfrmtReformat(d{i}, 'XYZCT', obj.rawDataFormat);
               end
 
               % cache
               obj.irawcelldata(cidload) = d;
            end
            
            % return requested cells
            d = obj.irawcelldata(cid);
            
         else
            
            % get raw file names
            [fn, tags] = obj.fileNameFromRawRange(range);
            if ischar(fn)
               fn = {fn};
            end
            tags = imfrmtRemoveRange(tags, obj.rawCellFormat);

            % allocate raw cell data
            d = cell(obj.rawCellSize(range));
            
            % load
            for i = 1:length(fn)
               d{i} = imreadBF(fn{i}, tags(i), 'squeeze', false, 'metadata', false);
               %size(d{i})
               d{i} = imfrmtReformat(d{i}, 'XYZCT', obj.rawDataFormat);
               %size(d{i})
            end
         end
      end
      
      
      function d = dataFromTagRange(obj, varargin)
         range = obj.fileTagRangeToIndexRange(varargin);
         d = obj.data(range);
      end
      
      
      function c = cellFromTagRange(obj, varargin)
         range = obj.fileTagRangeToIndexRange(varargin);
         c = obj.cell(range);
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