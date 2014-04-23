function saveEmbryoTrajectoryData(filename, frames, trajs, param)
%
% saveEmbryoTrajectoryData(dirname, frames, trajs)
%
% description:
%    Saves all information into asingle csv-file named filename
%    Each row containing the information in the source files pluss a
%    trajectory id by which the rows are sorted. 
%
% input:
%    filename    filename to write output to
%    frames      array of Frame classes
%    trajs       array of Trajectory classes
%    param       .print.save  print info on saving files
%
% See also: Frame, Trajectory, loadEmbryoData, saveEmbryoData

% label for id

id_label = 'Cell ID';

if nargin < 4
   param = [];
end

print_save = getParameter(param, {'print', 'save'}, 2);


% load data

nframes = length(frames);
ntrajs  = length(trajs);

for t = nframes:-1:1
   filedata{t} = importdata(frames(t).filename, ',', 1);
  
   if t < nframes
      if (size(filedata{t}.data,2) ~= dsize)
          error('saveEmbryoTrajectoryData: Error: inconsistent data size in file:%s\n', frames(t).filename)
      end
   else
      dsize = size(filedata{t}.data,2);
   end
   
   % find index for id in data file
   labels = strsplit(filedata{t}.textdata{1}, ',');
   [~, id_pos] = ismember(id_label,labels);
   if id_pos == 0
      error('saveEmbryoTrajectoryData: Error: cannot find id label in file:\n%s\n', frames(t).filename)
   end

   fids{t} = filedata{t}.data(:,id_pos); % all ids in file
   %nfids{t} = length(fids{t});   

   % find the object ids in file ids
   objids{t} = frames(t).id;
 
   [~, ids] = ismember(objids{t}, fids{t});
   if sum(ids > 0) ~= length(objids{t})
      error('saveEmbryoTrajectoryData.m: Error: index mismatch!\nCannot identify all objects in %s.', frames(t).filename);
   end
   objpos{t} = ids;   
   % objpos(1) is row position in file of first object in frames(t).objects  
end

% create large data table of the form
% Trajektory ID, Frame, Cell Id, ...
% sorted by trajectory ID and time frame

nobjects = sum( cellfun(@length, objids) );
data = zeros(nobjects, dsize + 2);
pos = 1;

for tid = 1:ntrajs

   timeids = trajs(tid).timeids;
   oids = trajs(tid).objids;
   
   for tl = 1:length(timeids)
      t = timeids(tl);
      data(pos,:) = [tid t filedata{t}.data(objpos{t}(oids(tl)), :)];
      pos = pos +1;
   end
end
   
 
% write data to file
writeCSV(filename, data, ['Trajectory ID, Frame, ' filedata{1}.textdata{1}]);

if print_save
   fprintf('saveEmbryoTrajectoryData: saved trajectory data to %s\n', filename)
end

end