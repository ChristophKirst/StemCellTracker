%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test ImageSourceTagged Class %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clear classes
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ImageSourceTagged - base class

is = ImageSourceTagged('./Test/Images/hESCells_Stack/W1F127T0001Z<z,2>C1.tif');
is.initialze();

