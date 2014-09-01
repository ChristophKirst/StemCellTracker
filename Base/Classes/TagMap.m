classdef TagMap < matlab.mixin.Copyable
   %
   % TagMap class handles tag and label conversions
   %
   % description:
   %    tags are substrings ('<tagname>') that represent numbers or strings in a tagformat string
   %    labels are shortcuts to acces a configuration of tags
   %
   % example:
   %    tagformat: imread('./image_c<channel>_z<slice>.tif') 
   %    then: replaceTags(tagformat, 'channel', 3, 'slice', 4) gives imread('./image_c3_z4.tif')
   %
   %     if lables has an entry 'DIC' -> {'channel', 0, 'slice', 15)
   %          replaceLabel(tagformat, 'DIC') gives  imread('./image_c0_z15.tif')
   %
   % See also: ImaqeSource, DataSource
   
   properties 
      tags   = containers.Map(); % map tagname -> default value
      label  = containers.Map(); % translation between names and tags, e.g. for channels 'dapi' -> {'channel', 1}
      default = '*';             % default substitution value if not in tags
   end

   methods
      function obj = TagMap(varargin) % constructor
         %
         % TagMap()
         % ImageSourceFile(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 0
            return
         elseif nargin == 1
            if isa(varargin{1}, 'TagMap') %% copy constructor
               obj = copy(varargin{1});
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
      
      
      function repl = replaceDefault(obj, tagformat)
         %
         % repl = replaceDefault(obj, varargin)
         %
         % description:
         %    replaces all tags with values given in varargin         
         
         tagnames = tagformat2tagnames(repl);
         def = {};
         for t = 1:length(tagnames)
            def = [def, tagnames{t}, obj.default]; %#ok<AGROW>
         end
         repl = tags2name(tagformat, def{:}); 
      end
      
      
      function repl = replaceTags(obj, tagformat, varargin)
         %
         % repl = replaceTags(obj, varargin)
         %
         % description:
         %    replaces all tags with values given in varargin
         
         repl = tags2name(tagformat, varargin{:});
      end
      
      function repl = replaceTagsDefault(obj, tagformat, varargin)
         %
         % repl = replaceTags(obj, varargin)
         %
         % description:
         %    replaces all tags with values given in varargin and default value otherwise
         
         repl = tags2name(tagformat, varargin{:});
      end
      
      
      
      
      function repl = replaceLabel(obj, varargin)
         % repl = replaceTags(obj, varargin)
         %
         % description:
         %
         %
         
         
      end
      
      function repl = replace(obj, varargin)
         %
         %
         
      end
      
      
      
   end
   
end