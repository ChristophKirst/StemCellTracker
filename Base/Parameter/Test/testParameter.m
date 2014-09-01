%%%%%%%%%%%%%%%%%%%%%%
%%% Test Parameter %%%
%%%%%%%%%%%%%%%%%%%%%%

clear all
clc


%% varargin2parameter

param = varargin2parameter('test', 6, 'hello', 'world', {'test', 100, 'moin.moin', 'moin'}, [], struct(), {})

param2 =  varargin2parameter(param, 'test', 1000, param)

%% parameter2cell

c = parameter2cell(param)
p = cell2parameter(c)


%% isemptyparameter

isemptyparameter([])
isemptyparameter({})
isemptyparameter(struct())

isemptyparameter(param)

%% parseParameter

p = getParameter([], 'x.y', 10)
p = getParameter(param, 'x.y', 100)


%%

param = parseParameter('doll', 10, 'sadads', 100)
param = parseParameter(param, 'x', 100, 'doll', 1000, 'z.o', 99)
param = parseParameter([], 'x', 10)

