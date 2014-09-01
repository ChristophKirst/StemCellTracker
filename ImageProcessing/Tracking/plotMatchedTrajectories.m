function plotMatchedTrajectories(frames, traj)
%
% plotMatchedTrajectories(frames, traj)
%
% description:
%    plot data points and trajectories from array of Frame and Trajectory classes
%
% input:
%    frames   array of Frame classes
%    traj     array of Trajectory classes
% 

% plot the data points

if ~isa(frames, 'Frame')
   error('plotMatchedTrajectories: expects array of Frame classes')
end


% figure

nframes = length(frames);
dim = frames.dim;
cols = hsv2rgb([linspace(0,1-1/nframes,nframes)' ones(nframes,2)]);

xyz = frames.r;

hold on
grid on

for t = 1: nframes
   s = frames(t).volume;
   s = s / max(s) .* 200;

   if dim == 2
      scatter(xyz{t}(1,:),xyz{t}(2,:), s, cols(t,:));
   else
      scatter3(xyz{t}(1,:),xyz{t}(2,:),xyz{t}(3,:), s, cols(t,:));
   end
end

clear xyz

% trajectories

if nargin > 1
   
   if ~isa(traj, 'Trajectory')
      error('plotMatchedTrajectories: expects array of Trajectory classes')
   end
   
   ntraj = length(traj);
   cols = hsv2rgb([linspace(0,1-1/ntraj,ntraj)' ones(ntraj,2)]);
   
   for i=1:ntraj

      xyzt = traj(i).r;
      
      if dim == 2
         line(xyzt(1,:), xyzt(2,:), 'Color', cols(i,:))
      else
         line(xyzt(1,:), xyzt(2,:), xyzt(3,:), 'Color', cols(i,:))
      end
      
   end
end


% some legends 
xlabel('x'); ylabel('y'); zlabel('z');
title('matched trajectories between object locations')

hold off


end

