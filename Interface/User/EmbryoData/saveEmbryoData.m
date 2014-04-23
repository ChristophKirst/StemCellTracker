function saveEmbryoData(dirname, frames, trajs, param)
%
% saveEmbryoData(dirname, frames, trajs)
%
% description:
%    Adds trajecktory data to the input files that generated data
%    and saves them under same name in the directory dirname.
%
% input:
%    dirname     directory to write output to
%    frames      array of Frame classes
%    trajs       array of Trajectory classes
%    param       .print.save  print info on saving files
%
% See also: Frame, Trajectory, loadEmbryoData

% label for id

id_label = 'Cell ID';


if nargin < 4
   param = [];
end

print_save = getParameter(param, {'print', 'save'}, 1);


% save data

if ~exist(dirname, 'dir')
   mkdir(dirname)
end

nframes = length(frames);

for t = 1:nframes
         
   % get full data from source file
   fd = importdata(frames(t).filename, ',', 1);
   fd.textdata{1} = [fd.textdata{1} 'Trajectory ID, Prev Cell ID, Next Cell ID, '];


   % create tracking data entries
   [trajpos, ~, objs, preobjs, postobjs] = trajs.frameSlicePrePost(t);
   
   
   % trajpos are the trajectory indices in traj with object at time t
   % tpos    are the positions of the objects in each trajectory in traj(trajpos)
   % objs    are the objects at time t
   % preobjs are the objects in previous time step
   % postobj are the objects in following time step
   
   % above ids based on ordering in data, restore real ids data.toIdsI()
   
   oids = [objs.id];
   preids = [preobjs.id];
   postids = [postobjs.id];

   % find index for id in data file
   labels = strsplit(fd.textdata{1}, ',');
   [~, id_pos] = ismember(id_label,labels);
   if id_pos == 0
      error('saveEmbryoData: Error: cannot find id label in file:\n%s\n', frames(t).filename)
   end
   
   fids = fd.data(:,id_pos); % all ids in file
   nfids = length(fids);   

   % find the object ids in file ids
   [~, ids] = ismember(oids, fids);
   if length(ids) ~= length(oids)
      error('saveEmbryoData.m: Error: index mismatch !\n');
   end
   % ids(1) is index in file of first object in frames(t).objects
   

   % assing trajectories to the indices
   ftids = - ones(1,nfids);
   ftids(ids) = trajpos;
   
   fpreids = - ones(1,nfids);
   fpreids(ids) = preids;
     
   fpostids = - ones(1,nfids);
   fpostids(ids) = postids;
   
   fd.data = [fd.data, ftids', fpreids', fpostids'];

   % filename
   [~, name, ext] = fileparts(frames(t).filename);
   filename = [dirname filesep name ext];

   
   % write data to file
   writeCSV(filename, fd.data, fd.textdata{1});

end

if print_save
   fprintf('saveEmbryoData: saved frames to %d files in %s\n', length(frames), dirname)
end

end