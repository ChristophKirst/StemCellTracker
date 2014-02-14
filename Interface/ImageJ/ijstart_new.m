function mij = ijstart(imagej_path)
%
% mij = ijstart(imagej_path)
%
% description:
%     sets up interface to ImageJ and returns MImageJ interface
%
% input:
%     imagej_path   path to ImageJ or Fiji
%
% See also: ijread, ijget, ijput, ijplot, ijplot3d
   
if nargin < 1
   imagej_path = '/home/ckirst/programs/Fiji';
end

ws = warning('off'); %#ok<WNOFF>
addjavapath(imagej_path);  % add this path for imagej.jar (if not using fiji)
addjavapath(fullfile(imagej_path, 'jars'));
addjavapath(fullfile(imagej_path, 'plugins'));
addjavapath(fileparts(mfilename('fullpath')));  % add this directory to make MImageJ available
warning(ws);

% Set the Fiji directory (and plugins.dir which is not Fiji.app/plugins/)
%java.lang.System.setProperty('ij.dir', fiji_directory);
%java.lang.System.setProperty('plugins.dir', fiji_directory);

mij = MImageJ.start();

end


%% Subfunction 
function addjavapath(directory)
    cpath = javaclasspath('-all');
    fns = dir(fullfile(directory, '*.jar'));
    jpath= cell(0);
    for i = 1:length(fns)
        if isempty(cell2mat(regexp(cpath, [filesep fns(i).name '$'])))
            jpath{end + 1} = fullfile(directory, fns(i).name); %#ok<AGROW>
        end
    end
    if ~isempty(jpath)
        javaaddpath(jpath, '-end');
    end
end
