classdef Segment < handle
%
% Segment class handling outline of objects etc
%
% See also: Object, Trajectory, Frame

   properties 
      shape = [];           % the compressed segmentation data
   end
   
   methods
      function obj = Segment(varargin)  % simple constructor
         %
         % Segment()
         % Segment(...,fieldname, fieldvalue,...)
         %
         if nargin == 1 && isa(varargin{1}, 'Segment') %% copy constructor
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
         end
      end
      
      
      function newobj = copy(obj)
         %
         % c = copy(obj)
         %
         % description:
         %     deep copy of the object
         
         nobjs = length(obj);
         newobj(nobjs) = Segment();
         props = properties(newobj);
         for k = 1:nobjs
            for i = 1:length(props)
               newobj(k).(props{i}) = obj(k).(props{i});
            end
         end
      end
      
      function store(obj, varargin)
         %
         % store(mask)
         % store(pixellist, imgsize)
         %
         % description:
         %     compress the label in mask and store them in shape property of this class
         
            obj.shape = bwcompress(varargin{:});
      end
    
      function [idx, isize] = pixelIdxList(obj)
         %
         % [idx, isize] = pixelIdxList(obj);
         %
         % description:
         %     return list of pixel indices that
         
         nobj = length(obj);
         idx = cell(1, nobj);
         for i = 1:nobj
            [idx{i}, isize] = bwuncompress(obj(i).shape, [], 'list');
         end
      end
         
      function m = mask(obj)
         [idx, isize] = obj.pixelIdxList();
         m = false(isize);
         m([idx{:}]) = true;
      end
      
   end % methods
   
end % classdef
      
      