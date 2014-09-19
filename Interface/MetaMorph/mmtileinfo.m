function [dims, wavenames]=getDimsFromScanFile(direc)
%[dims wavenames]=getDimsFromScanFile(direc)
%-----------------------------------------------
%Function to find dimensions of the tiling and
%the names of the wavelengths from the .scan file outputted by
%metamorph. 
%direc = metamorph output directory containing the .scan file
%   assumes only one .scan file in that directory. 
%dims = 2 component vector containing dimension
%wavenames = cell array containing wavelength names

verbose = 0;

fnames=dir([direc filesep '*scan']);
scanfile=[direc filesep fnames(1).name];
ff=fopen(scanfile);

strtofind='Row';
tline=fgetl(ff);
q=1;
while ischar(tline)
    k=strfind(tline,strtofind);
    if k
        if verbose
            disp(tline);
        end
        k1=strfind(tline,'_');
        k2=strfind(tline,'Col');
        k3=strfind(tline,'"');
        n1=str2double(tline((k+3):(k1-1)));
        n2=str2double(tline((k2+3):(k3(end)-1)));
        
        dims=[n1+1 n2+1];
    end
    m=strfind(tline,'WaveName');
    if m
        m1=strfind(tline,'"');
        wavenames{q}=tline((m1(end-1)+1):(m1(end)-1));
        q=q+1;
    end
    tline=fgetl(ff);
end