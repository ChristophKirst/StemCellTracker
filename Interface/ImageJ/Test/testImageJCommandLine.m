%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageJ commandLine tools %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

initialize
ijinitialize

%% create testfolder with images

img = syntheticLabeledImage([256, 256, 10], 20, 50);
img = mat2gray(img);

idir = '~/Desktop/TestIJ';
mkdir(idir)

ofn = fullfile(idir, 'Test<Z,2>.tif');

for z = 1:size(img,3);
   imwrite(img(:,:,z), tagExpressionToString(ofn, 'Z', z));
end


%% Start ImageJ

ip = ijpath

d = dirr([idir, '/*']);

ifdir = fullpath(idir);

ev = sprintf('run("Image Sequence...", "open=%s number=%d starting=1 increment=1 scale=100 file=[] sort");',...
   tagExpressionToString(fullfile(ifdir, 'Test<Z,2>.tif'), 'Z', 1), length(d))

fn = fullfile(ifdir, 'IJOpenMacro.txt')
fid = fopen(fn, 'w');
fprintf(fid, ev);
fclose(fid);

cmd = sprintf('cd %s; java -jar ij.jar -macro "%s" &', ip, fn)
system(cmd)


%%
'run("Stack to Hyperstack...", "order=xyczt(default) channels=3 slices=40 frames=5 display=Color");', ...
