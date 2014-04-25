%%%%%%%%%%%%%%%%%%%%%%
%%% Test Parameter %%%
%%%%%%%%%%%%%%%%%%%%%%

clear all
clc


%%

param = setParameter('x.y', 10, 'z.y', 100, 'z.o', 'sfsdf');
param.z


%%

p = getParameter([], 'x.y', 10)
p = getParameter(param, 'x.y', 100)


%%

param = parseParameter('doll', 10, 'sadads', 100)
param = parseParameter(param, 'x', 100, 'doll', 1000, 'z.o', 99)
param = parseParameter([], 'x', 10)

