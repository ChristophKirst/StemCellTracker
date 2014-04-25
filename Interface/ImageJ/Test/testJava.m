%%
javaclasspath

%%
clear all
javarmdynamicclasspath

%%
javaclasspath


%% compile
clear all
testname = [pwd, '/Test.java']
classname = [testname(1:end-4), 'class']
delete(classname)

res = system(['javac ' testname]);


%% 
%javaaddpath(classname, '-end')
javaaddpath(fileparts(classname), '-end')
javaclasspath

%%
tt = Test()
%%tt = javaObjectEDT('Test')

xx = repmat([1 2 3 4; 5 6 7 8], [1, 1, 3]);
tt.testArray(xx)

%%
ww = cast(xx, 'uint8');
class(ww)


tt.testArray(ww)



%%
%conclusion: arrays are translated to java 1-1
%size(xx) = [s1, s2, s3]