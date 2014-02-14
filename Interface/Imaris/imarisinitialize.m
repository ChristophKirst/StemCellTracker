function success = imarisinitialize(imarispath)
%
% success = imarisinitialize(imarispath)
%
% description:
%    initializes Imaris enviroment by making ImarisLib.jar from Imaris XT
%    visible to matlab
%
% input:
%    imarispath  (optional) path to Imaris ([] = automatic detection)
%
% output:
%    success     initialization of ImarisLib.jar succesfull
%

success = false;

% Imaris only runs on Windows and Mac OS X
if ~(ismac() || ispc())
    error('imarisinitialize: Imaris runs only on Windows or Mac OS X.');
end

if nargin < 1 || isempty(imarispath)
   imarispath = {};
else
   imarispath = {imarispath};
end


% Check if Imaris java libs are already installed 
jcp = javaclasspath('-all');
jpc = strfind(jcp, 'ImarisLib.jar');
jpc = ~cellfun(@isempty, jpc);

if any(jpc)  % ImarisLib.jar already installed
   fprintf('imarisinitialize: ImarisLib.jar already installed.'); 
   success = true;
   return;
end

% Search for ImarisLib.jar
ipath = imarispath(imarispath{:});

if isempty(ipath)
   fprintf('imarisinitialize: could not find ImarisLib.jar.'); 
   return;
else
   %add ImarisLib.jar to java libraries
   libjar = fullfile(ipath, 'XT', 'matlab', 'ImarisLib.jar');
   success = javaaddjar(libjar);
   if success; 
      fprintf('imarisinitialize: installed ImarisLib.jar.'); 
   else
      fprintf('imarisinitialize: could not find ImarisLib.jar.'); 
      return;
   end
end


% On Mac OS X, make sure that the Imaris Frameworks folder is at the
% beginning of the dynamic library path to prevent any conflicts with
% the Qt libraries loaded by MATLAB
if ismac()
    dylPath = getenv('DYLD_LIBRARY_PATH');
    imarisFrameworksPath = [ipath, filesep, 'Contents', filesep, 'Frameworks'];
    indx = strfind(dylPath, imarisFrameworksPath);
    if isempty(indx) || indx ~= 1
        % Prepend
        dylPath = [imarisFrameworksPath, ':', dylPath];
        setenv('DYLD_LIBRARY_PATH', dylPath)
    end
end


end