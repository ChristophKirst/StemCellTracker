function param = setParameter()
%
% param = setParameter()
%
% description:
%     generates struture setting the paramters to standard values
%

%% loading / saving
  
param.load.min = 10;            % skip data file if number of objects (cells|nuclei) is <=load.min (2)
param.load.change = 0.90;       % skip a loading file if fraction of objectnumber change is >=load.change

param.min_time = [];          % start matching from this time onward
param.max_time = [];          % do not match any data > this time
param.max_frames = [];        % stop matching after this number of frames

param.filter = [];            % apply a filter to each frame ([] = none)

param.save = 1;               % save the results

%% matching

% method

param.optimize = true;       % optimize match using an additional optimal coordinate transformation and second matching

% cost matrix

param.cost.creation  = [];   % costs for creating new objects ([])
param.cost.deletion = [];    % costs for deleting objects ([])

param.cutoff.cost = [];      % overall cost cutoff ([] = automatic)
param.cutoff.dist = [];      % cutoff for spatial distance ([] = automatic)
param.cutoff.time = -1;      % cutoff for temporal distances (-1 = ignore time differences)

param.weight.dist = 1.0;          % weight for distances in space (1.0)
param.weight.time.forward = 0.0;  % weight for positive time distances (0.0)
param.weight.time.backward = 0.0; % weight for negative time distances (0.0)
param.weight.volume = 0.0;        % weight for distances in volume (0.0)
param.weight.intensity = 0.0;     % weight for distances in intensity (0.0)
param.weight.type = Inf;          % weight for different types of objects (Inf)



%% printing

param.print.load = 1;           % print info about loading the data files
param.print.save = 1;           % print info about saving the data files

param.print.match.frames  = 1;  % print info about matching the frames
param.print.match.objects = 1;  % print results on matching the objects

param.print.estimates = 1;      % print automatic determined estimates for cutoffs



%% figures

param.figure.match = 0;          % generate figures for each match
param.figure.trajectories = 1;   % generate figure displaying the trajectories
param.figure.stats = 1;          % generate figure on statistics for the trajectories

