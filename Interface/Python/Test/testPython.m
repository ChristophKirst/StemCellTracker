%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Python Interface %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
py('eval', 'x = numpy.array([[[1,2,3],[4,5,6]], [[7,8,9], [10,11,12]]])')
py('eval', 'print x')
py('eval', 'print x.shape')
py('eval', 'print x[1,0,0]')

%%

x = pyget('x')
pyprint('x.shape')
size(x)

%%
py('eval', 'print x[0,0]')
x(:,1,1)


%%

pyprint('x[1,0,0]')
x(1,1,1)


%%

y = [5, 6, 7, 8, 9, 10]';

pyset('y', y)
pyprint('y')


%%

pyset('i', int32(6))
pyprint('i')


%%

y = [5, 6, 7; 8, 9, 10];
size(y)
pyset('y', y)
%pyeval('y = y.reshape(3,2)')
pyprint('y')
y


%%
z = pyget('y')
y

isequal(y,z)


%%
y(1,:)
pyprint('y[:,0]')


%%
size(y)
pyprint('y.shape')



%%

y = rand(2,3,4);
size(y)
pyset('y', y)
%pyeval('y = y.reshape(4,3,2)')
pyprint('y')
y



%%
y(1,3,2)
pyprint('y[1,2,0]')


%%
z= pyget('y')
size(z)
y


%%
%pyeval('z = y.swapaxes(1,0)')
%pyeval('z = y.reshape(2,3, order=''F'')')
%pyprint('z')



%% Stings

pyeval('s = ''abc''');
pyprint('s')



%%
s = pyget('s')

%%
pyset('t', s)
pyprint('t')
pyprint('t.__class__')


%% Cells

c = {1;2};
pyset('c', c)

pyprint('c.shape')

pyprint('c')

%%
c = {1; [1;2]}

pyset('c', c)
pyprint('c')


%%
c = {'hello'; 'world'}

pyset('c', c)
pyprint('c')
pyprint('c.shape')
pyprint('c[0]')


%%

c = { {1; 2}; 'hello'; 'world'; 100};
pyset('c', c)

pyprint('c')

%% Large data / Test for Memory leaks

x = rand(5000);

%%

pyset('f', x)
pyset('g', x)


%%

pyset('f', [])


%%

py('clear')


%%

compilePython


%% Structs

% dont work yet

st.x = 100;
st.y = 200;
st.z = 'hello world'
st.q = 0;


pyset('st', st)


%%
pyprint('st')
