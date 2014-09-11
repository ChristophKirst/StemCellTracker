 classdef ImageSourceTagged < ImageSource
   %
   % ImageSourceTagged class represents images accessed with the tag expression framework
   % 

   properties 
      icommand       = '';  % tagged command to load the data
      iinfocommand   = '';  % tagged command to get ImageInfo of each tagged individual data / image (e.g. 'imread_bf_info(<file>)') should always retunr consistent image format keeping singelton dimensions !
      ifilename      = '';  % tagged filename in case tag <file> appears in command or infocommand
      
      itagranges     = [];  % struct specifying the tag ranges and names (multiplicative tags assumed)
      itagformat     = {};  % specifies how tags are associated with certain dimensions, '' = cell coordinates -> ordered by order in tagranges
      itaginternal   = [];  % specifies the tags that are internally obtained via the read command (1) or need to be looped over externally (0 = default)
      
      %icommandformatordering  = ''; % format ordering in which the read routine returns the data, '' = pqlct
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
      end
      
      % infer properties from a list of image files
      function obj = fromFiles(obj, fname, varargin)
         param = parseParameter(varargin);
         [texpr, ~, tags]   = tagexpr(fname, param);
         obj.ifilename      = texpr;
         
         obj.itagranges     = tags2tagranges(tags);
         obj.itagformat     = obj.tagformatFromTags(tags);
         obj.itaginternal   = zeros(1, length(obj.itagformat)); % by default no internal parameter
         
         obj.icommand       = getParameter(param, 'command', 'imread_bf(<file>)');
         obj.iinfocommand   = getParameter(param, 'infocommand', 'imread_bf_info(<file>)');
 
         obj.XX % some sort of init
      end
      
      % infer properties from a file tagformat
      function obj = fromFileFormat(obj, filename, varargin)
         param = parseParameter(varargin);
         
         obj.ifilename      = filename;
         
         obj.itags          = tags2tagranges(tagexpr2tags(filename, param), 'check', true);
         obj.itagformat     = obj.tagformatFromTags(tags);
         obj.itaginternal   = zeros(1, length(obj.itagformat));
         
         obj.icommand       = getParameter(param, 'command', 'imread_bf(<file>)');
         obj.iinfocommand   = getParameter(param, 'infocommand', 'imread_bf_info(<file>)');
         
         obj.XX
      end
         
      % infer tag dimensions form the tags
      function tfrmt = tagformatFromTags(obj)
         tnames = obj.tagnames(obj);
         
         shortnames = num2cell('pqlct');
         longnames  = {'x', 'y', 'z', 'channel', 'time'};
         
         ntnames = length(tnames); 
         tfrmt = cell(1, ntnames);
         for i = 1:ntnames
            ids = or(ismember(shortnames, tnames{i}), ismember(longnames, tnames{i}));
            if any(ids)
               tfrmt{i} = shortnames{find(ids, 1, 'first')};
            end
         end
      end 
      
      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% tag routines
  
      % parse optional additonal tag input and reduce tag ranges accordingly
      function tagrs = parseTagRanges(obj, varargin)
         if nargin < 2
            tagrs = obj.itagranges;
         else
            tagrs = parseParameter(varargin);
            
            % sort tags and remove non-tags
            tagrs = rmfield(tagrs, setdiff(fieldnames(tagrs), filednames(obj.itagranges)));
            tagrs = catstruct(obj.itagranges, tagrs);
         end
      end
      
      % return tag values from index that lies in 1:ntags
      function tvs = ind2tagvalues(obj, id, varargin)
         tagrs = obj.parseTagRanges(varargin{:});
         tvs   = ind2tagvalues(tagrs, id);
      end
      
      function t = ind2tag(obj, id, varargin)         
         tagrs = obj.parseTagRanges(varargin{:});
         t     = ind2tag(tagrs, id);
      end

      function tnames = tagnames(obj)
         tnames = fieldnames(obj.itagranges);
      end  
      
      function s = tagrangesize(obj, varargin)
         tagrs = obj.parseTagRanges(varargin{:});
         s     = tagrangesize(tagrs);
      end
      
      function n = ntags(obj)
         n = prod(obj.tagrangesize);
      end
      
      
      
% Todo: handle cell size/datasize changes !!      
%       function obj = addTags(obj, tags)
%          obj.itagranges = catstruct(obj.itagranges, tags);
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
         cmd = strrep(obj.icommand, '<file>', obj.ifilename);
         cmd = tagexpr2string(cmd, varargin);
      end

      % returns command with tag index i given tag sepcs
      function cmd = ind2command(obj, i, varargin)
         cmd = obj.command(obj.ind2tag(i, varargin{:}));
      end
       
      
      function cmd = infocommand(obj, varargin)
         cmd = strrep(obj.iinfocommand, '<file>', obj.ifilename);
         cmd = tagexpr2string(cmd, varargin);
      end
     
      function cmd = ind2infocommand(obj, i, varargin)
         cmd = obj.infocommand(obj.ind2tag(i, varargin{:}));
      end
      
      
      function fn = filename(obj, varargin)
         fn = tagexpr2string(obj.ifilename, obj.parseTagRanges(varargin));
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
         info = eval(cmd);
      end
      

      function [info, objinfo] = getInfo(obj, varargin) 
         info = ImageInfo(varargin{:});
         
         % get info of the first data set
         dinfo = getDataInfo(); % info for a single tag, internal tags are singeltons !
         
         % internal daata format for internal data tags being singeltons
         singletagintfrmt = dinfo.iformat;
         objinfo.isingletagformat = singletagintfrmt;


         % full internal data format 
         tint = obj.taginternal;
         tfrmt = obj.tagformat;
         tfrmt(cellfun(@length, tfrmt)==0) = {'u'}; % replace empty slots with dummy -> extenion: with uvwrs 
         itfrmt = cell2mat(tfrmt(tint));
         
         if any(ismember(singletagintfrmt, itfrmt))
            error('%s: getInfo: raw internal data and internal data tag formats overlap: %s <> %s', class(obj), var2char(singletagintfrmt), var2char(itfrmt));
         end
         
         fullintfrmt = setdiff('pqlct', setdiff('pqlct', [singletagintfrmt, itfrmt]));
          
         
         % full format including external tags contributing to data
         etfrmt = setdiff(tfrmt, itfrmt);
         fullfrmt = setdiff('pqlct', setdiff('pqlct', [fullintfrmt, etfrmt]));
         
         info.iformat = fullfrtm;
         info.irawformat = fullfrmt;
          
         % data and cell size
         dsi = dinfo.isize;
         tsi = tagrangesize(obj.itagranges);
         
         % ids of single tag data dims in full format
         [stdids, stdpos] = ismember(fullfrmt, singletagintfrmt); stdpos = stdpos(stdids);
         % ids of tags in full format
         [tids, tpos] = ismember(fullfrmt, cell2mat(tfrmt)); tpos = tpos(tids);
         
         isize = ones(1, length(fullformat));
         isize(stdids) = dsi(stdpos);
         isize(tids)   = tsi(tpos);
         
         obj.isize = isize;
         
         info.pqlctsizeFromFormatAndSize();
         
         % cell format and size
         obj.icellformat = ''; % sorted by occurence in tag expression 
         obj.icellsize = tsi(setdiff(1:length(tsi), tpos)); % remainig tag sizes 
      end
      
      

        
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data access
      
      
      % ids in tagranges that map to data dims (both internal and external data tags)
      % optional output pos is the position in iformat of the corresponding tag ids 
      function [ids, pos] = datatagids(obj)
         tfrmt = obj.itagformat;
         [ids, pos] = ismember(tfrmt, num2cell('pqlct'));
         pos = pos(ids);
         %ids = find(ids);
      end
      
      % return the tag names of the data dimensions (both internal and external data tags)
      function dn = datatagnames(obj)
         dn = obj.fieldnames;
         dn = dn(obj.datatagids);
      end
      
      
      % return the tag ranges of the data dimensions (both internal and external data tags)
      function tgrs = datatagsranges(obj, varargin)
         tgrs = obj.parseTagRanges(varargin{:});
         ids = obj.celltagids;
         fnames = fieldnames(tgrs);
         tgrs = rmfield(tgrs, fnames(ids));
      end

      % ids in tagranges that map to cell dims
      function [ids, pos] = celltagids(obj)
         tfrmt = obj.itagformat;
         ids = ~ismember(tfrmt, num2cell('pqlct'));
         %if nargout > 1
         %   pos = 1:sum(ids); % we dont have uvwrs coords implemented yet -> change here !!
         %end
      end
      
      % return the tag ranges of the data dimensions
      function cn = celltagnames(obj)
         cn = obj.fieldnames();
         cn = cn(obj.celltagids);
      end

      % return the tag ranges of the cell dimension
      function tgrs = celltagranges(obj, varargin)
         tgrs = obj.parseTagRanges(varargin{:});
         ids = obj.datatagids;
         fnames = fieldnames(tgrs);
         tgrs = rmfield(tgrs, fnames(ids));
      end
       
      
      % ids of internal / external data tags that map to data dims
      function ids = datatagidsinternal(obj, varargin)
         ids = obj.datatagids(varargin{:});
         ids = and(ids, obj.taginternal);
      end
      
      function [tgrs, ids] = datatagrangesinternal(obj, varargin)
         tgrs = obj.parseTagRanges(varargin{:});
         ids = obj.datatagidsinternal(tgrs);
         fnames = fieldnames(tgrs);
         tgrs = rmfield(tgrs, fnames(~ids));
      end
      
      
      function ids = datatagidsexternal(obj, varargin)
         ids = obj.datatagids(varargin{:});
         ids = and(ids, ~obj.taginternal);
      end

      function [tgrs, ids] = datatagrangesexternal(obj, varargin)        
         tgrs = obj.parseTagRanges(varargin{:});
         ids = obj.datatagidsexternal(tgrs);
         fnames = fieldnames(tgrs);
         tgrs = rmfield(tgrs, fnames(~ids));
      end
      

      %%% data 
      
      function d = data(obj, varargin)
         tgrs = obj.parseTagRanges(varargin{:});
         
         % check if all cell tags are singeltons
         ctrgs = obj.celltagranges(tgrs);
         csi = tagrangesize(ctrgs);   
         if any(csi ~= 1)
            error('%s: data: in order to return data array all cell tags %s need to be specified', class(obj), var2char(obj.celltagnames))
         end

         % obtain internal and external data tags and single cell tags
         [datatrgsint, intids] = obj.datatagrangesinternal(tgrs);
         [datatrgsext, extids] = obj.datatagrangesexternal(tgrs);
         %datatgrs    = catstruct(datatrgsint, datatrgsext);

         
         % determine size of the data
         si = obj.size;     % full data size
         frmt = obj.format; % full data format
     
         tsi = tagrangesize(tgrs); 
         [tids, tpos] = ismember(obj.itagformat, frmt); tpos = tpos(tids);
         si(tpos) = tsi(tids);

         d = zeros(si);

         % constant internal data tags
         consttrgs = parseParameter(datatrgsint, ctrgs);
         
         % prepare constant part of data signment format
         asgn = repmat({':'}, 1, length(frmt));
         
         % internal tag sizes
         
         for n = 
         
         % ids that get replaced by indices given the external tag id
         asgnids =        
         
         

         % load data: loop over external tags
         nexttags = prod(tagrangesize(datatrgsext));
         for t = 1:nexttags
            
            datatagext  = ind2tag(datatrgsext, t);
            
            cmd = obj.commad(consttrgs, datatagext); 
            dr = eval(cmd);


         % replace specified external tags with corresponding ids
         
         
         
         % remove singelton internal dimensions
         

         
         
         
     
         tsi = tagrangesize(tgrs); 
         [tids, tpos] = ismember(obj.itagformat, frmt); tpos = tpos(tids);
         si(tpos) = tsi(tids);
         
         
         
         
         
         
         % remove singeltons and for full tag ranges put {:}
         singletagfrmt = obj.iinfo.irawformat;
         
         
 
      
   
         % determine cell size
         
         
         % determine data size
         
         
         % load data: loop over tags and asign to correct data
         
         
         
         
         
         
                    asgn = obj.dataAsignment(datatrgsint, extdatatag);
            d(asgn{:}) = dr;
         end
   
      end
      
      
         
      end
      
      
      
      
      
      
      
      

      % format of image data

      % data part of the image format
      % cell part of the tag format
      
%       function texpr = formattag(obj)
%          % tag part of the image format
%          texpr = cell2mat(obj.itagdims);
%          texpr = setdiff(texpr, 'uvwrs');
%       end
%       
%       function dfrmt = formatdata(obj)
%          dfrmt = setdiff(obj.format, obj.formattag);
%       end
      
%       function cfrmt = formatcell(obj)
%          cfrmt = 'uvwrs';
%          cfrmt = cfrmt(1:(length(obj.itagdims) - length(obj.formattag)));
%       end
      
      


      % positions
%       % pos (and member ship) of tags in image format
%       % pos of data in image format
% 
%       function p = tagpos(obj)  % tags not part of image data have pos 0 !
%          p = ismember(obj.tagdims, num2cell(obj.format));
%       end
%       
%       function p = datapos(obj)
%          p = obj.tagpos;
%          p = p(
%       end

      
      
      
            % indices: 

%       function ids = tagids(obj)
%          % indices of tags dims in image format
%          ids = find(ismember(obj.format, obj.formattag));
%       end
%       
%       function ids = dataids(obj)
%          % indices of data dims in image format
%          ids = find(~ismember(obj.format, obj.formattag));
%       end
% 
      % sizes:
      % size of tag dimensions
      % size of tag dimensions in image
      % size of data dimensions in image
      % size of cell dimensions
      % number of tags in total


%       
%       function s = tagsize(obj)
%          s = obj.size;
%          s = s(obj.tagids);
%       end
 
%       function s = datasize(obj)
%          s = obj.size(ids);
%          s = s(obj.dataids);
%       end
%       
%       function s = sizecell(obj)
%          s = obj.tagrangesize(ids);
%          s = s(obj.celltagids);
%       end


      
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %%% data access
      

      

  
%       % size of the cell array that will be returned given optional tag specs 
%       function cs = cellsize(obj, varargin)
%          tagrs = obj.parseTagRanges(varargin{:});
%          cs = tagsize(tagrs);
%          ids = cellfun(@length, setdiff(obj.tagdims, {'u', 'v', 'w', 'r', 's'})) == 0;
%          cs = cs(ids);
%          cs(cs == 1) = [];
%       end
%       
%       % size of the data array that will be returned given optional tag specs 
%       function ds = datasize(obj, varargin)
%          
%          %          s = obj.size(ids);
%          %          s = s(obj.dataids);
%          
%          tagrs = obj.parseTagRanges(varargin{:});
% 
%          ds = tagsize(tagrs);
%          ids = ~(cellfun(@length, setdiff(obj.tagdims, {'u', 'v', 'w', 'r', 's'})) == 0);
%          cs = cs(ids);
%          cs(cs == 1) = [];
%       end
      
      


 
      function s = size(obj, varargin)
         if nargin == 1
            s = size@ImageSource(obj);
         else
            trgs = parseTagRanges(varargi{:}n);
            
            tsize = obj.tagsize(tag);

            
         end
      end
      

      
      
      
      
      
      
      
      
      
      
      
      




      
      
      
      function info = infoString(obj)
         info = infoString@ImageSource(obj);
         info = [info, '\nfilename:   ', obj.ifilename];
         info = [info, '\ncommand:    ', obj.icommand];
         info = [info, '\ninfocommand:', obj.iinfocommand]; 
         info = [info, '\ndataformat: ', obj.idataformat];
         info = [info, '\ntagexpr:  ', obj.itagexpr];
         info = [info, '\nntags:      ', var2char(obj.ntags)];
         info = [info, '\ntagsize:    ', var2char(obj.tagrangesize)];
      end

   end
   
   
end