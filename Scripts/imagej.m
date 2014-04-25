%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ImageJ Interface and Setup %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% To set up imagej you need 
%  - a running version of imagej with the 3d viewer or fiji on your computer (fiji installs the 3d viewer automatically)
%  - javas j3d installed 
% 


%% Find ImageJ

ijpath

% if this fails specifiy the base directory in which image is installed 
% ijpath(dir)

%% Initialize ImageJ

ijinitialize


%% Test if j3d is installed

addpath('./Utils/External/Fiji')
IsJava3DInstalled


%% Install j3d

InstallJava3D


%% Test 3d visualization

ijplot3d(rand(30,30,10,3))

