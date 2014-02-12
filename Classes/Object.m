classdef Object < handle
%
% Object class storing info for the objects identified in images and to be tracked / quantified 
%
% See also: Trajectory, Picture, Frame
   properties
      r = [0; 0];    % spatial position
      time = 0;      % time
      
      volume = [];    % area or volume ([] for none)
      intensity = []; % intensity ([] for none)
      id = [];        % id ([] forn none)
      
      type = [];      % data that identifies the object type ([] for none)
      
      %source = [];   % reference to Image class to which this object belongs ([] for none)
   end
   
   %properties (Dependent)
   %   length         % length of data entries obtained when using toArray
   %end
   
   methods
      function obj = Object(data, d, varargin)
         %
         % Object(id, time, r, size, intensity, type)
         % Object(data, dim)
         %
         % input: id, r, time,...   class field entries 
         %        data              column of the from [id time r ...]
         %        dim               spatial dimension
         %
         if nargin > 0
            if nargin ==1 && isa(data, 'Object')  % copy constructor 

               obj.id = data.id;
               obj.time = data.time;
               obj.r = data.r;
               obj.volume = data.volume;
               obj.intensity = data.intensity;
               obj.type = data.type;

            elseif nargin == 2 && isnumeric(d) && length(data) >= 2 + d % Object(data, dim)
               
               if size(data,1) == 1
                  data = data';
               end
               
               obj.id = data(1);
               obj.time = data(2);
               obj.r = data(3:2+d);
               
               n = length(data);
               if n > 2 + d
                  obj.volume = data(3+d);
               end   
               if n > 3 + d
                  obj.intensity = data(4+d);
               end              
               if n > 4 + d
                  obj.type = data(5+d:end);
               end
               
            elseif isnumeric(data) % Object(id, time, r, size, intensity, type)
               
               obj.id = data;
               
               if nargin > 1
                  obj.time = d;
               end
               if nargin > 2
                  obj.r = varargin{1};
               end
               if nargin > 3
                  obj.volume = varargin{2};
               end
               if nargin > 4
                  obj.intensity = varargin{3};
               end
               if nargin > 5
                  obj.type = [ varargin{4:end} ];
               end
            end
         end
      end
                 
      function c = copy(obj)
      %
      % c = copy(obj)
      %
      % deep copy of the object
      %   
         n = length(obj);
         for i=n:-1:1
            c(i) = Object(obj(i));
         end
         
      end
      
      function d = dim(obj)
         if ~isempty(obj)
            d = size(obj(1).r, 1);
         else
            obj
            error('sadsad');
         end
      end
      
      
      %function l = get.length(obj)
      %      l = 1 + obj(1).dim + ~isempty(obj(1).volume) + ~isempty(obj(1).intensity) + ~isempty(obj(1).id) + size(obj(1).type,1);
      %end
      
      function data = toArray(obj)
      %
      % data = toArray(obj)
      %
      % convert data to column vector
      %    
         if length(obj) == 1
            data = [obj.time; obj.r; obj.volume;  obj.intensity; obj.id; obj.type];
         else
            
            d = [];
            if ~isempty(obj(1).volume) 
               d = [obj.volume];
            end
            if ~isempty(obj(1).intensity)
               if ~isempty(d)
                  d = vertcat(d, [obj.intensity]);
               else
                  d = [obj.intensity];
               end
            end

            if ~isempty(obj(1).id)
               if ~isempty(d)
                  d = vertcat(d, [obj.id]);
               else
                  d = [obj.id];
               end
            end 

            if ~isempty(obj(1).type)
               if ~isempty(d)
                  d = vertcat(d, [obj.type]);
               else
                  d = [obj.id];
               end
            end  

            if ~isempty(d)    
               data = vertcat( [obj.time],  [obj.r], d);
            else
               data = vertcat( [obj.time],  [obj.r]);
            end
         end       
      end % toArray
      
 
      function obj = transformCoordinates(obj, R, T, C)
      %
      % obj = transformCoordinates(obj, R, T, C)
      %
      % applies rotation R, scaling C  and translation T to coordinates r
      %  
         n = length(obj);
         
         if n==1
            obj.r = C * R * obj.r + T;
         else
            trnsf = C * R * [ obj.r ] + repmat(T,1,n);
            for i=1:n
               obj(i).r = trnsf(:,i);
            end
         end
         
      end
     
      

   end % methods
end % classdef
      
      