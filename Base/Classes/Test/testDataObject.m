%%%%%%%%%%%%%%%%%%%%%%%
%%% Test DataObject %%%
%%%%%%%%%%%%%%%%%%%%%%%

%%  
clc
clear all
clear classes

o(5) = DataObject();

%%
o.dataFields()

%%

s = struct('a', num2cell( 1:5))

%%
o.addData(s)

o.dataFields()

%% channel data

o.addChannelData(1, s)

o.dataFields()


%%

o.channel(1, 'a')