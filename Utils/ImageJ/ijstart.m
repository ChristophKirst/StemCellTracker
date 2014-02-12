function mij = ijstart(ijpath)
%
% ijstart(ijpath)
%
% description:
%    starts ImageJ interface
%
% input:
%    ijpath    path to ImageJ or Fiji
%



%%% Fiji version
% dirn = pwd;
% 
% if nargin < 1
%    ijpath = '/home/ckirst/programs/Fiji/scripts/';
% end
% 
% addpath(ijpath);
% 
% Miji(false); %run startup script from Fiji
% 
% cd(dirn)

% return

%%% version for mij

if nargin < 1
   ijpath = '/usr/share/imagej/';
end

addjavapath(ijpath);

mij = ij.ImageJ([], 2);

end




%% Subfunction 
function addjavapath(directory)
    cpath = javaclasspath('-all');

    fns = dir(fullfile(directory, '*.jar'));
    jpath = {};
    for i = 1:length(fns)
        if isempty(cell2mat(regexp(cpath, [filesep fns(i).name '$'])))
            jpath{end + 1} = fullfile(directory, fns(i).name); %#ok<AGROW>
        end
    end
    if ~isempty(jpath)
        javaaddpath(jpath, '-end');
    end
    
    
    fns = dir(directory);
    for i = 1:length(fns)
      if fns(i).isdir && ~strcmp(fns(i).name, '.') && ~strcmp(fns(i).name, '..') 
           addjavapath(fullfile(directory, fns(i).name));
      end
    end
    
end
