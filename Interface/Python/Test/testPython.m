close all
clear all

initialize
compilePython

%%
py('eval', 'import numpy')
py('eval', 'x = numpy.array([1,2,3,4])')
py('eval', 'print x')


%%

py('eval', 'import numpy')
py('eval', 'x = numpy.array([[[1,2],[3,4]], [[7, 8], [9, 10]]])')
py('eval', 'print x')

%%

x = pyget('x')


%%

y = [5, 6, 7, 8, 9, 10];

pyset('y', y')
py('eval', 'print y')





%% large data ??


x = rand(5000);



%%

pyset('f', x)
pyset('g', x)


%%

pyset('f', [])


%%

py('exit')
