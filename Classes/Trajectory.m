classdef Trajectory < handle
%
% Trajectory class for storing and reshaping trajectory data
%
% See also: Frame, Object

   properties       

      objects  = [];   % objects that constitute the trajectory
      frameids = [];   % id of the frames to which the objects belong
      objids   = [];   % ids of object in each time frame
      
      %timeseries = [];  % Reference to time series class
      %experiment = [];  % Reference to experiment class
   end 
   
   properties (Dependent)
      dim             % spatial dimension
      duration        % length of the trajectory
   end
   
   methods
      
      function obj = Trajectory(varargin)
         %
         %
         % Trajectory()
         % Trajectory(trajectory)
         % Trajectory(...,fieldname, fieldvalue,...)
         %
         
         if nargin == 1 && isa(varargin{1}, 'Trajectory') %% copy constructor
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
      % f = copy(obj)
      %
      % description:
      %    deep copy
      %
         nobjs = length(obj);
         newobj(nobjs) = Frame();
         for k = 1:nobjs
            newobj(k).objects    = obj(k).objects.copy;
            newobj(k).frameids   = obj(k).frameids;
            newobj(k).objids     = obj(k).objids;     
         end 
      end
      
 
      function d = get.dim(obj)
            d = obj(1).objects(1).dim;
      end
      
      function l = get.duration(obj)
         if length(obj) == 1
            l = length(obj.objects);
         else
            l = arrayfun(@length, obj.objects);
         end
      end
      
      
      function t = time(obj)
      %
      % t = time(obj)
      %
      % output:
      %   t    time points of the objects in the trajectory
      %
         if length(obj) > 1 % for array of trajectories
            t = cellfun(@(x) [ x.time ], { obj.objects }, 'UniformOutput', false);
         else               % single image
            t = [ obj.objects.time ];
         end   
      end
      
      
      function xyz = r(obj)
      %
      % xyz = r(obj)
      %
      % output:
      %   xyz    coordinates of the objects in the trajectory
      %
         if length(obj) > 1 % for array of trajectories
            xyz = cellfun(@(x) [ x.r ], { obj.objects }, 'UniformOutput', false);
         else               % single trajectory
            xyz = [ obj.objects.r ];
         end   
      end

      function vol = volume(obj)
      %
      % vol = volume(obj)
      %
      % output:
      %   vol    volumes of the objects in the trajectory
      %
         if length(obj) > 1 % for array of trajectories
            vol = cellfun(@(x) [ x.vol ], { obj.objects }, 'UniformOutput', false);
         else               % single trajectory
            vol = [ obj.objects.volume ];
         end   
      end
      
      function i = intensity(obj)
      %
      % i = intensity(obj)
      %
      % output:
      %   i    intensities of the objects in the trajectory
      %
         if length(obj) > 1 % for array of trajectories
            i = cellfun(@(x) [ x.intensity ], { obj.objects }, 'UniformOutput', false);
         else               % single image
            i = [ obj.objects.intensity ];
         end   
      end

      function t = type(obj)
      %
      % t = type(obj)
      %
      % output:
      % t    type data of the objects in the trajectory
      %
         if length(obj) > 1 % for array of trajectories
            t = cellfun(@(x) [ x.type ], { obj.objects }, 'UniformOutput', false);
         else               % single trajectory
            t = [ obj.objects.type ];
         end   
      end

      function i = id(obj)
      %
      % i = id(obj)
      %
      % output:
      %   i    coordinates of the objects in the trajectory
      %
         if length(obj) > 1 % for array of trajectories
            i = cellfun(@(x) [ x.id ], { obj.objects }, 'UniformOutput', false);
         else               % single trajectory
            i = [ obj.objects.id ];
         end   
      end   
      
      function [trajpos, tpos, objs] = timeSlice(obj, t)
      %
      % [trajpos, tpos, objs] = timeSlice(obj, t)
      %
      % description:
      %   finds all trajectories ids tids and objects objs passing through time t
      %
      % input:
      %   t     time 
      %
      % output:
      %
      %   trajpos  positions in trajectory array obj of trajectories that contain an object at time t
      %   tpos     positions in each trajectory in obj(trajpos) of objects at time t
      %   objs     objects in the trajectories at time t
      %

         if length(obj) == 1 % single trajectory
            tpos = find(obj.objects.time == t);
            %tpos = find(obj.frameids == t);
            if isempty(tpos)
               trajpos = tops;
            else
               trajpos = 1;
            end   
            if nargout > 2
               objs = obj.objects(tpos);
            end   
            
         else % array of trajectories 
            tpos = cellfun(@(x) find(x==t,1), obj.time, 'UniformOutput', false);
            %tpos = cellfun(@(x) find(x==t,1), {obj.frameids}, 'UniformOutput', false);
            trajpos = find(cellfun(@length, tpos)); 
            if nargout > 1              
               tpos = cell2mat(tpos);
            end
            
            if nargout > 2
               n = length(trajpos);
               trajobjs = { obj(trajpos).objects };
               for i = n:-1:1
                  objs(i) = trajobjs{i}(tpos(i));
               end
            end

         end
      end
      
      
      function [trajpos, tpos, objs] = frameSlice(obj, t)
      %
      % [trajpos, tpos, objs] = frameSlice(obj, t)
      %
      % description:
      %   finds all trajectories ids tids and objects objs passing through frame t
      %
      % input:
      %   t     time 
      %
      % output:
      %
      %   trajpos  positions in trajectory array obj of trajectories that contain an object at time t
      %   tpos     positions in each trajectory in obj(trajpos) of objects at time t
      %   objs     objects in the trajectories at time t
      %

         if length(obj) == 1 % single trajectory
            %tpos = find(obj.objects.time == t);
            tpos = find(obj.frameids == t);
            if isempty(tpos)
               trajpos = tops;
            else
               trajpos = 1;
            end   
            if nargout > 2
               objs = obj.objects(tpos);
            end   
            
         else % array of trajectories 
            %tpos = cellfun(@(x) find(x==t,1), obj.time, 'UniformOutput', false);
            tpos = cellfun(@(x) find(x==t,1), {obj.frameids}, 'UniformOutput', false);
            trajpos = find(cellfun(@length, tpos)); 
            if nargout > 1              
               tpos = cell2mat(tpos);
            end
            
            if nargout > 2
               n = length(trajpos);
               trajobjs = { obj(trajpos).objects };
               for i = n:-1:1
                  objs(i) = trajobjs{i}(tpos(i));
               end
            end

         end
      end
      
      
      
      
      function [trajpos, tpos, objs, preobjs, postobjs] = timeSlicePrePost(obj, t, nullobj)
      %
      % [trajpos, tpos, objs, preobjs, postobjs] = timeSlicePrePost(obj, t, dt, nullobj)
      %
      % description:
      %   finds all trajectories tids and objects objs passing through frame t
      %   together with nearest neightbout objects
      %
      % input:
      %   t        time 
      %   nullobj  object return if trajectory ends or starts ([])
      %
      % output:
      %
      %   trajpos  positions of the trajectories in array obj that contain an object at time t
      %   tpos     positions in each trajectory in obj(trajpos) of the object at time t
      %   objs     objects in the trajectories at time t
      %   preobjs  objects at previous time point
      %   postobj  objects at next time point
      %
      % See also: timeSlice
      

         if nargin < 3
            nullobj = [];
         end

         if isempty(nullobj)
            nullobj = Object();
            nullobj.id = -1;
         end
         
      
         if length(obj) == 1
            tpos = find(obj.objects.time == t);
            %tpos = find(obj.frameids == t);
            if isempty(tpos)
               trajpos = tpos;
            else
               trajpos = 1;
            end
            
            if nargout > 2
               objs = obj.objects(tpos);
            end

            if nargout > 3
               if ~isempty(tpos) && tpos > 1
                  preobjs = obj.objects(tpos-1);
               else
                  preobjs = nullobj;
               end
            end

            if nargout > 3
               if ~isempty(tpos) && tpos < length(obj.objects)
                  postobjs = obj.objects(tpos+1);
               else
                  postobjs = nullobj;
               end
            end   
            
           
         else % array of trajectories
            
            tpos = cellfun(@(x) find(x==t,1), obj.time, 'UniformOutput', false);               
            %tpos = cellfun(@(x) find(x==t,1), {obj.frameids}, 'UniformOutput', false);
            trajpos = find(cellfun(@length, tpos));
            
            if nargout > 1
               tpos = cell2mat(tpos);
            end
            
            
            if nargout > 2
               n = length(trajpos);
               trajobjs = { obj(trajpos).objects };
               for i = n:-1:1
                  objs(i) = trajobjs{i}(tpos(i));
               end
            end

            if nargout > 3
               preobjs = repmat(nullobj, 1,n);
               indx = find(tpos > 1);
               for i = indx
                  preobjs(i) =  trajobjs{i}(tpos(i)-1);
               end
            end
            
            if nargout > 4
               postobjs = repmat(nullobj, 1,n);
               indx = find( tpos < cellfun(@length, trajobjs) );
               for i = indx
                  postobjs(i) =  trajobjs{i}(tpos(i)+1);
               end            
            end
         end
  
      end
      
      function [trajpos, tpos, objs, preobjs, postobjs] = frameSlicePrePost(obj, t, nullobj)
      %
      % [trajpos, tpos, objs, preobjs, postobjs] = frameSlicePrePost(obj, t, dt, nullobj)
      %
      % description:
      %   finds all trajectories tids and objects objs passing through frame t
      %   together with nearest neightbout objects
      %
      % input:
      %   t        time 
      %   nullobj  object return if trajectory ends or starts ([])
      %
      % output:
      %
      %   trajpos  positions of the trajectories in array obj that contain an object at time t
      %   tpos     positions in each trajectory in obj(trajpos) of the object at time t
      %   objs     objects in the trajectories at time t
      %   preobjs  objects at previous time point
      %   postobj  objects at next time point
      %
      % See also: timeSlice
      

         if nargin < 3
            nullobj = [];
         end

         if isempty(nullobj)
            nullobj = Object();
            nullobj.id = -1;
         end
         
      
         if length(obj) == 1
            %tpos = find(obj.objects.time == t);
            tpos = find(obj.frameids == t);
            if isempty(tpos)
               trajpos = tpos;
            else
               trajpos = 1;
            end
            
            if nargout > 2
               objs = obj.objects(tpos);
            end

            if nargout > 3
               if ~isempty(tpos) && tpos > 1
                  preobjs = obj.objects(tpos-1);
               else
                  preobjs = nullobj;
               end
            end

            if nargout > 3
               if ~isempty(tpos) && tpos < length(obj.objects)
                  postobjs = obj.objects(tpos+1);
               else
                  postobjs = nullobj;
               end
            end   
            
           
         else % array of trajectories
            
            %tpos = cellfun(@(x) find(x==t,1), obj.time, 'UniformOutput', false);               
            tpos = cellfun(@(x) find(x==t,1), {obj.frameids}, 'UniformOutput', false);
            trajpos = find(cellfun(@length, tpos));
            
            if nargout > 1
               tpos = cell2mat(tpos);
            end
            
            
            if nargout > 2
               n = length(trajpos);
               trajobjs = { obj(trajpos).objects };
               for i = n:-1:1
                  objs(i) = trajobjs{i}(tpos(i));
               end
            end

            if nargout > 3
               preobjs = repmat(nullobj, 1,n);
               indx = find(tpos > 1);
               for i = indx
                  preobjs(i) =  trajobjs{i}(tpos(i)-1);
               end
            end
            
            if nargout > 4
               postobjs = repmat(nullobj, 1,n);
               indx = find( tpos < cellfun(@length, trajobjs) );
               for i = indx
                  postobjs(i) =  trajobjs{i}(tpos(i)+1);
               end            
            end
         end
  
      end
      
      
      
      function data = toArray(obj)
      %
      % data = toArray(obj)
      %
      % returns data array of objects in trajectory
        
         data = obj.objects.toArray;
      end
      
      
      function s = startObjects(obj)  
      %
      % s = startObjects(obj)  
      %
      % returns starting objects of the trajectories
      %
         if length(obj) == 1
            s = obj.objects(1);
         else
            s = cell2mat(cellfun(@(x) x(1), obj.objects));
         end
      end

            
      function e = endObjects(obj)  
      %
      % e = endObjects(obj)  
      %
      % returns starting objects of the trajectories
      %
         if length(obj) == 1
            e = obj.objects(end);
         else
            e = cell2mat(cellfun(@(x) x(end), obj.objects));
            
            %e = cell2mat(cellfun(@(x,y) [x(end), y(end)]', obj.times, obj.ids,'UniformOutput', false))';
         end
      end
 
      function stats = statistics(obj)
      %
      % stats = statistics(obj, data)
      %
      % returns some statistics of the trajectories obj
      %
      
         stats = statisticsTrajectory(obj);
         
      end
      
      
      function imglab = labeledImage(obj)
         if length(obj) ==1 
            objs = obj.objects;
            shape = [objs.segment];
            [idx, isize] = shape.pixelIdxList();
            imglab = zeros(isize);
            for i = 1:length(idx)
               imglab(idx{i}) = obj.frameids(i);
            end
         else
            warning('%s: no methods to generate a labeled image form multiple trajectories', class(obj));
            imglab = [];
         end
      end
      
 
   end % methods
end % classdef
      
      
      
      
      