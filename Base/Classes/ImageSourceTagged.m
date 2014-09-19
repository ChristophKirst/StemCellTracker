 classdef ImageSourceTagged < ImageSource
   %
   % ImageSourceTagged class represents images accessed with the tag expression framework
   % 

   properties 
      ireadcommand   = '';  % tagged command to load the data
      iinfocommand   = '';  % tagged command to get ImageInfo of each tagged individual data / image (e.g. 'imread_bf_info(<file>)') should always retunr consistent image format keeping singelton dimensions !
      
      ifilename      = '';  % tagged filename in case tag <file> appears in command or infocommand
      ibasedirectory = '';  % base directory of the files (possibly tagged)
      
      itagranges     = [];  % struct specifying the tag ranges and names (multiplicative tags assumed)
      itagformat     = {};  % specifies how tags are associated with certain dimensions, '' = cell coordinates -> ordered by order in tagranges
      itaginternal   = [];  % specifies the tags that are internally obtained via the read command (1) or need to be looped over externally (0 = default)
      
      %ireadformatordering  = ''; % format ordering in which the read routine returns the data, '' = pqlct
      %isingletagformat = []; % format the read routine returns data when all internal tags are singeltons
   end

   methods
      function obj = ImageSourceTagged(varargin) % constructor
         %
         % ImageSourceTagged()
         % ImageSourceTagged(taggedfilename)
         % ImageSourceTagged(filenames)
         % ImageSourceTagged(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'ImageSourceTagged') %% copy constructor
               obj = copy(varargin{1});
            elseif ischar(varargin{1})
               obj.fromFileFormat(varargin{1});
            elseif iscellstr(varargin{1})
               obj.fromFiles(varargin{1});
            else
               error('%s: invalid constructor input, expects char at position %g',class(obj), 1);
            end
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
         end
         
         obj.icache = false; % for tagging simple caching is usaually not a good idea
      end
      
      % infer properties from a list of image files
      function obj = fromFiles(obj, fname, varargin)
         param = parseParameter(varargin);
         [texpr, ~, tags]   = tagexpr(fname, param);
         obj.ifilename      = texpr;
         
         obj.itagranges     = tags2tagranges(tags);
         obj.itagformat     = obj.tagformatFromTagRanges(tags);
         obj.itaginternal   = zeros(1, length(obj.itagformat)); % by default no internal parameter
         
         obj.ireadcommand   = getParameter(param, 'command', 'imread_bf(''<file>'')');
         obj.iinfocommand   = getParameter(param, 'infocommand', 'imread_bf_info(''<file>'')');
 
         obj.iinfo          = obj.getInfo();
      end
      
      % infer properties from a file tagformat
      function obj = fromFileFormat(obj, filename, varargin)
         param = parseParameter(varargin);
         
         obj.ifilename      = filename;
         
         obj.itagranges     = tags2tagranges(tagexpr2tags(filename, [], param), 'check', true);
         obj.itagformat     = obj.tagformatFromTagRanges();
         obj.itaginternal   = zeros(1, length(obj.itagformat));
         
         obj.ireadcommand   = getParameter(param, 'command', 'imread_bf(''<file>'')');
         obj.iinfocommand   = getParameter(param, 'infocommand', 'imread_bf_info(''<file>'')');
         
         obj.iinfo          = obj.getInfo();
      end
      
      
      % infer settings form read command and filename -> readcommand tags are considered internal data tags !
      function obj = fromReadcommandAndFilename(obj, varargin)
         param = parseParameter(varargin);
         
         filename = obj.ifilename;
         readcmd  = obj.ireadcommand;
         
         if isempty(obj.iinfocommand)
            if ~isempty(strfind(readcmd, 'imread_bf'))
               obj.iinfocommand = strrep(readcmd, 'imread_bf', 'imread_bf_info');
            else
               obj.iinfocommand = getParameter(param, 'infocommand', 'imread_bf_info(''<file>'')');
            end
         end

         intdtags = setdiff(tagexpr2tagnames(readcmd), {'file'});
         
         % the range should be predifined -> ways to infer those but for speed etc we keep it simple
         tgrs = obj.itagranges;
         
         if isempty(tgrs)
            error('%s: cannot infer internal data ranges, specify range of internal data tag rages in itagranges!', class(obj));
         elseif ~isempty(setdiff(intdtags, fieldnames(tgrs)))
            error('%s: cannot infer internal data ranges, specify range of internal data tags %s used in read command!', class(obj), var2char(setdiff(intdtags, fieldnames(tgrs))));
         end

         obj.itagranges     = tags2tagranges(tagexpr2tags(filename, [], param), 'check', true);
         obj.itagranges     = parseParameter(obj.itagranges, tgrs);
 
         obj.itagformat     = obj.tagformatFromTagRanges(); 
         obj.itaginternal   = ismember(fieldnames(obj.itagranges), intdtags)';

         obj.iinfo          = obj.getInfo();
         
      end
      

      % infer tag dimensions form the tags
      function tfrmt = tagformatFromTagRanges(obj)
         tnames = obj.tagnames;
         
         shortnames = num2cell('pqlct');
         longnames  = {'x', 'y', 'z', 'channel', 'time'};
         
         cellnames = 'uvwrs';
         ci = 1;
         
         ntnames = length(tnames); 
         tfrmt = cell(1, ntnames);
         for i = 1:ntnames
            ids = or(ismember(shortnames, tnames{i}), ismember(longnames, tnames{i}));
            if any(ids)
               tfrmt{i} = shortnames{find(ids, 1, 'first')};
            else
               %tfrmt{i} = '';
               tfrmt{i} = cellnames(ci);
               ci = ci + 1;
            end
         end
      end
      
      
      function initializeTagformat(obj)
         cellnames = 'uvwrs';
         tfrmt = obj.itagformat;
         cellnames = setdiff(cellnames, cell2mat(tfrmt));
         
         k = 1;
         for i = 1:length(tfrmt)
            if isempty(tfrmt{i})
               tfrmt{i} = cellnames(k);
               k = k + 1;
            end  
         end
         
         obj.itagformat = tfrmt;
      end
      
      
      function initializeInfo(obj, varargin)
         obj.iinfo = obj.getInfo(varargin{:});
      end

      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% tag routines
      
      
      function obj = setTagRange(obj, name, range)
      %
      % obj = setTagRange(obj, name, range)
      %

         if ~ismember(name, obj.tagnames)
            error('ImagesourceTagged: setTagRange: %s is not a tag!', name);
         end

         if ~iscell(range)
            range = num2cell(range);
         end
         obj.itagranges.(name) = range;
       
         obj.iinfo = obj.getInfo();
      end
      

      % parse optional additonal tag input and reduce tag ranges accordingly
      function tagrs = parseTagRanges(obj, varargin)
         if nargin < 2
            tagrs = obj.itagranges;
         else
            tagrs = parseParameter(varargin);
            
            % sort tags and remove non-tags
            tagrs = rmfield(tagrs, setdiff(fieldnames(tagrs), fieldnames(obj.itagranges)));
            tagrs = parseParameter(obj.itagranges, tagrs);
         end
      end
      
      % return tag values from index id in 1:ntags
      function tvs = ind2tagvalues(obj, id, varargin)
         tagrs = obj.parseTagRanges(varargin{:});
         tvs   = ind2tagvalues(tagrs, id);
      end
      
      function t = ind2tag(obj, id, varargin)         
         tagrs = obj.parseTagRanges(varargin{:});
         t     = ind2tag(tagrs, id);
      end

      function tnames = tagnames(obj)
         if isempty(obj.itagranges)
            tnames = {};
         else
            tnames = fieldnames(obj.itagranges);
         end
      end  
      
      function s = tagrangesize(obj, varargin)
         tagrs = obj.parseTagRanges(varargin{:});
         s     = tagrangesize(tagrs);
      end
      
      function n = ntags(obj, varargin)
         tagrs = obj.parseTagRanges(varargin{:});
         n = prod(tagrangesize(tagrs));
      end
      
      
      
% Todo: handle cell size/datasize changes !!      
%       function obj = addTags(obj, tags)
%          obj.itagranges = parseParameter(obj.itagranges, tags);
%       end
% 
%       function obj = removeTags(obj, tagnames)
%          obj.itagranges = rmfield(obj.itagranges, tagnames);
%       end
% 
%       function obj = setTags(obj, tags)
%          obj.itagranges = tags;
%       end
%       
%       function tags = getTags(obj)
%          tags = obj.itagranges;
%       end
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% commands

      % returns command format by replacing specified tags      
      function cmd = command(obj, varargin)
         cmd = strrep(obj.ireadcommand, '<file>', obj.filename);
         cmd = tagexpr2string(cmd, varargin{:});
      end

      % returns command with tag index i given tag sepcs
      function cmd = ind2command(obj, i, varargin)
         cmd = obj.command(obj.ind2tag(i, varargin{:}));
      end
       
      
      function cmd = infocommand(obj, varargin)
         cmd = strrep(obj.iinfocommand, '<file>', obj.filename);
         cmd = tagexpr2string(cmd, varargin{:});
      end
     
      function cmd = ind2infocommand(obj, i, varargin)
         cmd = obj.infocommand(obj.ind2tag(i, varargin{:}));
      end
      
      
      function fn = filename(obj, varargin)
         fn = tagexpr2string(fullfile(obj.ibasedirectory, obj.ifilename), varargin{:});
      end    
      
      function fn = ind2filename(obj, i, varargin)
         fn = obj.filename(obj.ind2tag(i, varargin{:}));
      end
      

      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% setup and getInfo
      
      function info = getDataInfo(obj)
         % single image info
         tag = obj.ind2tag(1);         % first tag specs -> internal tags are singeltons and assumedto be removed by the infocommand in returned info.iformat and isize!
         cmd = obj.infocommand(tag);
         fprintf('ImageSourceTagged: evaluating command: %s\n', cmd);
         info = eval(cmd);
      end
      

      function info = getInfo(obj, varargin) 
         info = ImageInfo(varargin{:});
         
         % get info of the first data set
         dinfo = obj.getDataInfo(); % info for a single tag, internal tags are singeltons !
         
         % internal raw data format for internal data tags being singeltons
         singletagintfrmt = dinfo.irawformat;
         %objinfo.isingletagformat = singletagintfrmt;

         % full internal data format 
         tint = obj.itaginternal > 0;
         tfrmt = obj.itagformat;
         if any(tint) 
            itfrmt = [tfrmt{tint}];
         else
            itfrmt = '';
         end
         
         if any(ismember(singletagintfrmt, itfrmt))
            error('%s: getInfo: raw internal data and internal data tag formats overlap: %s <> %s', class(obj), var2char(singletagintfrmt), var2char(itfrmt));
         end
         
         fullintfrmt = setdiff('pqlct', setdiff('pqlct', [singletagintfrmt, itfrmt]), 'stable');
          
         
         % full format including external tags contributing to data
         etfrmt = setdiff(cell2mat(tfrmt), itfrmt, 'stable');
     
         fullfrmt = setdiff('pqlct', setdiff('pqlct', [fullintfrmt, etfrmt]), 'stable');
         
         info.idataformat = fullfrmt;
         info.irawformat  = fullfrmt;
          
         % data and cell size
         dsi = dinfo.idatasize;
         tsi = tagrangesize(obj.itagranges);
         
         % ids of single tag data dims in full format
         [stdids, stdpos] = ismember(fullfrmt, singletagintfrmt); stdpos = stdpos(stdids);
         % ids of tags in full format
         [tids, tpos] = ismember(fullfrmt, cell2mat(tfrmt)); tpos = tpos(tids);
         
         isize = ones(1, length(fullfrmt));
         isize(stdids) = dsi(stdpos);
         isize(tids)   = tsi(tpos);
         
         info.idatasize = isize;
         info.irawsize  = isize;
         
         info.pqlctsizeFromFormatAndSize();
         
         % cell format and size
         tfrmt = cell2mat(tfrmt);
         info.icellformat = setdiff(tfrmt, 'pqlct', 'stable');
              
         % sorted by occurence in tag expression 
         info.icellsize = tsi(setdiff(1:length(tsi), tpos)); % remainig tag sizes 
         
         % copy other info
         info.iseries    = dinfo.iseries;
         info.inimages   = dinfo.inimages;
         info.idataclass = dinfo.idataclass;
         info.imetadata  = dinfo.imetadata;
         info.iscale     = dinfo.iscale;
         info.iunit      = dinfo.iunit;
         info.icolor     = dinfo.icolor;
      end
      
      

        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data access
      
      % assumption on the read command is that it returns raw data ordered accroding to pqlct
      % with singelton dimensions removed
      
      % ids in tagranges that map to raw dims (both internal and external data tags)
      % optional output pos is the position in iformat of the corresponding tag ids and zero if not member of the data  
      function [ids, pos] = rawtagids(obj)
         tfrmt = obj.itagformat;
         [ids, pos] = ismember(tfrmt, num2cell('pqlct'));  % ids are the ids of pqlct in tfrmt, pos(ids) ar the positions of tfrmt in pqlct
         %pos = pos(ids);
         %ids = find(ids);
      end
      
      % return the tag names of the data dimensions (both internal and external data tags)
      function dn = rawtagnames(obj)
         dn = obj.tagnames();
         dn = dn(obj.rawtagids);
      end
      
      function dn = datatagnames(obj)
          dn = obj.rawtagnames;
      end

      
      % return the tag ranges of the data dimensions (both internal and external data tags)
      function [tgrs, ids, pos] = rawtagsranges(obj, varargin)
         tgrs = obj.parseTagRanges(varargin{:});
         [ids, pos] = obj.celltagids;
         fnames = fieldnames(tgrs);
         tgrs = rmfield(tgrs, fnames(ids));
      end

      % ids in tagranges that map to cell dims
      % optional output pos is the position in icellformat of the corresponding tag ids 
      function [ids, pos] = celltagids(obj)
         tfrmt = obj.itagformat;
         ids = ~ismember(tfrmt, num2cell('pqlct'));
         if nargout > 1
            pos = zeros(1,length(ids));
            pos(ids) = 1:sum(ids); % we dont have uvwrs coords implemented yet -> change here !!
         end
      end
      
      % return the tag ranges of the data dimensions
      function cn = celltagnames(obj)
         cn = obj.tagnames();
         cn = cn(obj.celltagids);
      end

      % return the tag ranges of the cell dimension
      function [tgrs, ids, pos] = celltagranges(obj, varargin)
         tgrs = obj.parseTagRanges(varargin{:});
         [ids, pos] = obj.rawtagids;
         fnames = fieldnames(tgrs);
         tgrs = rmfield(tgrs, fnames(ids));
      end
       
      
      % internal / external data tags that map to data dims
      function [ids, pos] = rawtagidsinternal(obj)
         [ids, pos] = obj.rawtagids();
         ids = and(ids, obj.itaginternal);
         pos(~ids) = 0;
      end
      
      function [tgrs, ids, pos] = rawtagrangesinternal(obj, varargin)
         tgrs = obj.parseTagRanges(varargin{:});
         [ids, pos] = obj.rawtagidsinternal();
         fnames = fieldnames(tgrs);
         tgrs = rmfield(tgrs, fnames(~ids));
      end
      
      
      function [ids, pos] = rawtagidsexternal(obj, varargin)
         [ids, pos] = obj.rawtagids(varargin{:});
         ids = and(ids, ~obj.itaginternal);
         pos(~ids) = 0;
      end

      function [tgrs, ids, pos] = rawtagrangesexternal(obj, varargin)        
         tgrs = obj.parseTagRanges(varargin{:});
         [ids, pos] = obj.rawtagidsexternal();
         fnames = fieldnames(tgrs);
         tgrs = rmfield(tgrs, fnames(~ids));
      end
      

      
      %%% data 
      
      
      %getData routine: assumptions: tag framework set
      %                              read command returns pqlct ordering with singeltons removed
      %                              asssumes iinfo is initialized (e.g. via getInfo)
      function d = getData(obj, varargin)
         
         tgrs = obj.parseTagRanges(varargin{:});
         
         % check if all cell tags are singeltons
         ctrgs = obj.celltagranges(tgrs);
         csi   = tagrangesize(ctrgs);   
         if any(csi ~= 1)
            error('%s: data: in order to return data array all cell tags %s need to be specified', class(obj), var2char(obj.celltagnames))
         end
         ctag = ind2tag(ctrgs, 1);

         % obtain internal and external data tags and single cell tags
         [rawtagrsint, ~]      = obj.rawtagrangesinternal(tgrs);
         [rawtagrsext, extids] = obj.rawtagrangesexternal(tgrs);
         %datatgrs    = parseParameter(rawtagrsint, rawtagrsext);

         
         % determine size of the data
         si   = obj.iinfo.irawsize;     % full raw data size
         frmt = obj.iinfo.irawformat; % full raw data format
         if isempty(frmt)
            si   = obj.iinfo.idatasize;
            frmt = obj.iinfo.idataformat;
         end
     
         tfrmt = obj.itagformat;
         
         tsi = tagrangesize(tgrs); 
         [tids, tpos] = ismember(tfrmt, num2cell(frmt)); tpos = tpos(tids);
         si(tpos) = tsi(tids);

         d = zeros(si);  % initialize data using pqlct ordering, no singelton removed

         % constant internal data tags
         consttrgs = parseParameter(rawtagrsint, ctag);
         
         % prepare constant part of data asignment format
         asgn = repmat({':'}, 1, length(frmt));
         
         % external data positions
         [~, epos] = ismember(frmt, [tfrmt{extids}]);
         epos = epos > 0;

         % load data: loop over external tags
         nexttags = prod(tagrangesize(rawtagrsext));
         for t = 1:nexttags
            
            [datatagext, dids] = ind2tag(rawtagrsext, t);
            
            cmd = obj.command(consttrgs, datatagext); 
            
            %fprintf('ImageSourceTagged: evaluating command: %s\n', cmd);
            dr = eval(cmd);

            % replace specified external tags with corresponding ids
            asgn(epos) = num2cell(dids);
            d(asgn{:}) = dr;
         end

         % order according to output format
         d = impqlpermute(d, frmt, obj.info.idataformat);
         
         % remove singeltons
         d = squeeze(d);
         
      end
      
      function c = getCellData(obj, varargin)
         tgrs = obj.parseTagRanges(varargin{:});
         
         % check if all cell tags are singeltons
         ctrgs = obj.celltagranges(tgrs);
         csi   = tagrangesize(ctrgs);
         if length(csi) == 1
            c = cell(csi,1);
         else
            c = cell(csi); 
         end
            
         ncells = prod(csi);

    
         % obtain internal and external data tags and single cell tags
         [rawtagrsint, ~]      = obj.rawtagrangesinternal(tgrs);
         [rawtagrsext, extids] = obj.rawtagrangesexternal(tgrs);

         % determine size of the data
         si   = obj.iinfo.irawsize;     % full raw data size
         frmt = obj.iinfo.irawformat; % full raw data format
         if isempty(frmt)
            si   = obj.iinfo.idatasize;
            frmt = obj.iinfo.idataformat;
         end
     
         tfrmt = obj.itagformat;
         
         tsi = tagrangesize(tgrs); 
         [tids, tpos] = ismember(tfrmt, num2cell(frmt)); tpos = tpos(tids);
         si(tpos) = tsi(tids);

         d = zeros(si);  % initialize data using pqlct ordering, no singelton removed
         
         % prepare constant part of data asignment format
         asgn = repmat({':'}, 1, length(frmt));
         
         % external data positions
         [~, epos] = ismember(frmt, [tfrmt{extids}]);
         epos = epos > 0;
         
         nexttags = prod(tagrangesize(rawtagrsext));
         
         for ci = 1:ncells
            ctag = ind2tag(ctrgs, ci);

            % constant internal data tags
            consttrgs = parseParameter(rawtagrsint, ctag);
         
            % load data: loop over external tags

            for t = 1:nexttags
            
               [datatagext, dids] = ind2tag(rawtagrsext, t);
            
               cmd = obj.command(consttrgs, datatagext); 
            
               %fprintf('ImageSourceTagged: evaluating command: %s\n', cmd);
               dr = eval(cmd);

               % replace specified external tags with corresponding ids
               asgn(epos) = num2cell(dids);
               d(asgn{:}) = dr;
            end

            % order according to output format
            d = impqlpermute(d, frmt, obj.info.idataformat);
         
            % remove singeltons
            d = squeeze(d);
            
            c{ci} = d;
         end
      end
      
      
      function d = data(obj, varargin)
         if nargin > 1 && isnumeric(varargin{1}) % acces by id
            tag = obj.ind2tag(varargin{:});
            d = obj.getData(tag);
         else
            d = obj.getData(varargin{:});
         end
      end
            
         
      function d = celldata(obj, varargin)
         if nargin > 1 && isnumeric(varargin{1}) % acces by id
            tag = obj.ind2tag(varargin{:});
            d = obj.getCellData(tag);
         else
            d = obj.getCellData(varargin{:});
         end
      end

      
      
      % data routine should work out of the box using getData
      
      
      function si = datasize(obj, varargin)
         if nargin < 2
            si = obj.iinfo.idatasize;
         else
            si = obj.iinfo.idatasize;
            tgrs = obj.parseTagRanges(varargin{:});
            tsi = tagrangesize(tgrs);
            [tids, tpos] = ismember(obj.itagformat, num2cell(obj.dataformat));
            tpos = tpos(tids);
            si(tpos) = tsi(tids);
            si(si==1) = [];
         end
      end
      
      
      function df = dataformat(obj, varargin)
         if nargin < 2
            df = obj.iinfo.idataformat;
         else
            df = obj.iinfo.idataformat;
            tgrs = obj.parseTagRanges(varargin{:});
            tsi = tagrangesize(tgrs) < 2;
            [tids, tpos] = ismember(obj.itagformat, num2cell(df));
            tids = and(tsi, tids > 0);
            df(tpos(tids)) = [];
         end
      end
      
      
      function cs = cellsize(obj, varargin)
         if nargin < 2
            cs = obj.iinfo.icellsize;
         else
            ctrgs = obj.celltagranges(varargin{:});
            cs   = tagrangesize(ctrgs);
         end
      end

      function cs = cellformat(varargin)
         cs = '';
      end
      
      function plot(obj, varargin)
         imgs = obj.celldata(varargin{:});
         if nargin> 1 && isnumeric(varargin{1})
            varargin = varargin(2:end);
         end
         implottiling(imgs, varargin{:});
      end
     

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      function info = infoString(obj)
         info = infoString@ImageSource(obj);
         info = [info, '\nfilename:    ', obj.filename];
         info = [info, '\nreadcommand: ', obj.ireadcommand];
         info = [info, '\ninfocommand: ', obj.iinfocommand]; 

         info = [info, '\nntags:       ', var2char(obj.ntags)];
         info = [info, '\ntagnames:    ', var2char(obj.tagnames)];
         info = [info, '\ntagrangesize:', var2char(obj.tagrangesize)];
         info = [info, '\ntagformat:   ', var2char(obj.itagformat)];
         info = [info, '\ntaginternal: ', var2char(obj.itaginternal)];
      end

   end
   
   
end