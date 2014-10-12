 classdef ImageSourceBF < ImageSource
   %
   % ImageSourceBF class provdes access to images in bio-format
   % 

   properties 
      ireader        = [];  % instance of loci.formats.ChannelSeparator
      ifilename      = '';  % bio-format file name
   end

   methods
      function obj = ImageSourceBF(varargin) % constructor
         %
         % ImageSourceBF()
         % ImageSourceBF(filename)
         % ImageSourceBF(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageSourceBF') %% copy constructor
               obj = copy(varargin{1});
            elseif ischar(varargin{1})
               obj.fromFile(varargin{1});
            else
               error('%s: invalid constructor input, expects char at position %g',class(obj), 1);
            end
         else
            obj.fromParmeter(varargin{:});
         end
         
         obj.icache = false;  % acces to images is almost as fast, also files are typically very big
         
         %bfinitialize;
      end
      
      function delete(obj)
         if ~isempty(obj.ireader)
            %fprintf('%s: closing file: %s\n', class(obj), obj.ifilename);
            obj.ireader.close();
         end
      end 
      
      function obj = fromParameter(obj, varargin)
         obj = classFromParameter(obj, 'i', varargin{:});
      end
      
      function obj = fromFile(obj, filename, varargin)
         %
         % obj = fromFile(obj, filename, varargin)
         %
         % description:
         %      initialize source from a given file name

         if ~exist(filename, 'file')
            error('%s: fromFile: file %s does not exists!', class(obj), filename);
         end

         % initialize the reader
         obj.ifilename       = filename;
         
         % read form specific series if specified
         rawRange = obj.reshapeRangeToRawRange(varargin);
         if isfield(rawRange, 'S')
            s = {'S', rawRange.S};
         elseif isfield(rawRange, 's')
            s = {'s', rawRange.s};
         else
            s = {};
         end

         [info, obj.ireader] = imreadBFInfo(filename, s{:});
         obj.fromImageInfo(info);
      end
      
      function obj = initializeFromSeries(obj, s)
         obj.initializeFromFile(obj.ifilename, 'S', s);
      end
      
      function obj = fromImageInfo(obj, info)
         fromImageInfo@ImageInfo(obj, info);
      end

      function obj = initializeRangeFromRaw(obj, varargin)
         obj.irange = obj.rawRange(varargin{:});
      end

      function obj = initializeFromTiling(obj, varargin)
         % automatic tiling
         [ts, tf] = obj.tileSizeAndFormat(varargin{:});
         obj.setReshape('S', 'UV', ts);
         obj.setCellFormat(tf);
      end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % bf functions
      
      function fn = currentFilename(obj)
         fn = char(obj.ireader.getCurrentFile());
      end
      
      function tp = type(obj)
         tp = char(obj.ireader.getFormat());
      end
      
      function fn = filename(obj)
         fn = obj.ifilename;
      end
      
      function ir = reader(obj)
         ir = obj.ireader;
      end
      
      function pos = stagePositions(obj, varargin)
         rgs = obj.reshapeRangeToRawRange(imfrmtParseRangeLast(obj.irange, varargin));
         pos = imreadBFStagePositions(obj.ireader, rgs);
         pos = imfrmtReshapeCellData(pos, 'XY', 'ZCTS', obj.dataFormat, obj.cellFormat, obj.reshapeFrom, obj.reshapeTo, obj.reshapeSize);
      end
      
      function pos = timeStamps(obj, varargin)
         pos = imreadBFTimeStamps(obj.ireader, varargin{:});
      end
      
      function pos = exposureTimes(obj, varargin)
         pos = imreadBFExposureTimes(obj.ireader, varargin{:});
      end
  
      function [tsi, tpos] = tileSize(obj,varargin)
         [tsi,tpos] = imreadBFTileSize(obj.ireader, varargin{:});
      end
      
      function [tsi, tfrmt] = tileSizeAndFormat(obj,varargin)
         [tsi, tfrmt] = imreadBFTileSizeAndFormat(obj.ireader, varargin{:});
      end

      function ts = waveLengths(obj,varargin)
         ts = imreadBFWaveLengths(obj.ireader, varargin{:});
      end

      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % ranges
      
      function tr = rawRange(obj, varargin)
         tr = imreadBFRange(obj.ireader, varargin{:});
      end


      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data access

      function d = getRawData(obj, varargin)
         % range
         range = imfrmtParseRangeLast(obj.irange, varargin{:});
         
         % check if all cell dims are singeltons
         cid = obj.rawCellIndex(range);
         if length(cid) > 1
            error('%s: getRawData: cell dimensions are not specified to be singeltons!', class(obj));
         end
        
         d = imreadBF(obj.ireader, range, 'squeeze', false, 'metadata', false);
         d = imfrmtReformat(d, 'XYZCT', obj.rawDataFormat);  
      end
      
      function d = getRawCellData(obj, varargin) 
         % range
         range = imfrmtParseRangeLast(obj.irange, varargin{:});
         
         d = imreadBF(obj.ireader, range, 'squeeze', false, 'metadata', false); 
         
         if ~iscell(d)
            d = {d};
         end
         d = imfrmtReformatCellData(d, 'XYZCT','S', obj.rawDataFormat, obj.rawCellFormat);
         
%          disp BF2
%          size(d)
      end
      

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % plotting 
      
      function plot(obj, varargin)
         imgs = obj.cell(varargin{:});
         implot(imgs{1}, varargin{:});
      end
      
      function plotCells(obj, varargin)
         imgs = obj.cell(varargin{:});
         if nargin> 1 && isnumeric(varargin{1})
            varargin = varargin(2:end);
         end
         implottiling(imgs, varargin{:});
      end

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % info
      
      function info = infoString(obj)
         info = infoString@ImageSource(obj, 'BF');
         info = [info, '\nchannel:   '  var2char(obj.channelName)];
         info = [info, '\nfilename:  ', var2char(obj.filename)];
         info = [info, '\ntype:      ', var2char(obj.type)];
      end

   end
   
   
end