%%%%%%%%%%%%%%%%%
%%% Test Tags %%%
%%%%%%%%%%%%%%%%%

clc
clear all
close all

initialize


%% tagformat2tagnames

texpr = 'test_T<tag1>-<tag2, 5>-<tag1,3>.tif';
[tnames, tagsplit, taginfo] = tagexpr2tagnames(texpr)

taginfo(1)


%% tagexpr2string
tagexpr2string(texpr, 'tag1', 2, 'tag2',  5)
tagexpr2string(texpr, 'tag1', 2)


%% taginfo2tagexpr

taginfo(1).tag = strrep(taginfo(1).tag, taginfo(1).name, 'test');
taginfo(1).name = 'test';

taginfo2tagexpr(tagsplit, taginfo)


%% reducing tags

texpr = 'test_T<tag1>-<tag2, 5>-<tag1,3><tag3>.tif';

clear tags
tags(1).tag1 = 3;
tags(2).tag1 = 2;
val = {1,2};
[tags.tag2] = val{:};
[tags.tag3] = val{:};
[tags.tag4] = val{:};

tagsnew = tagsreduce(tags)


%%
clc
[tagsnew, texprnew] = tagsreduce(tags, texpr)


%% tagexpr2tags

clc
texpr = 'test_T<tag1>-<tag2, 5>-<tag1,3>_<tag3>.tif';
name  = 'test_T01-00003-001_19.tif';
tagexpr2tags(texpr, name)

%%
name  = {'test_T01-00003-001_19.tif', 'test_T01-00004-001_20.tif'};
tags = tagexpr2tags(texpr, name)

tags(2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tagexpr: infer tags from list names or files

clc
[texpr, tnames, tags] = tagexpr('./Test/Images/hESCells_Tiling/*', 'tagnames', {'field', 'test'})
[tags.field]

%%
clc
[texpr, tnames, tags] = tagexpr('./Test/Images/hESCells_Stack/*', 'tagnames', 'z')

[tags.z]

%%
clc
[texpr, tnames, tags] = tagexpr('./Test/Images/mESCells_Wnt/*', 'tagnames', {'t', 'z'});

texpr
tags

tags(5)


%% more complicated file stucture

dirr('./Test/Images/hESCells_Folder/*/*.tif')


%%
[texpr, tnames, tags] = tagexpr('./Test/Images/hESCells_Folder/*/*.tif', 'reduce', false);

texpr
tags

[tags.tag1]
[tags.tag2]
[tags.tag3]
[tags.tag4]


[texpr, tnames, tags] = tagexpr('./Test/Images/hESCells_Folder/*/*.tif', 'reduce', true, 'tagnames', {'t1', 't2'})

%%

texpr = './Test/Images/hESCells_Folder/f<tag1,1>_t<tag2,1>/f<tag1,1>_t<tag2,2>_z<tag3,2>.tif';

tagexpr2files(texpr)

tagexpr2filename(texpr)

%% tagexpr2files - remove no consistent file names if multiple occurences of tags
clc

texpr = './Test/Images/hESCells_Folder/f<tag1,1>_t<tag2,1>/f<tag1,1>_t<tag2,2>_z<tag3,2>.tif';

dirr(tagexpr2filename(texpr))
tagexpr2files(texpr)


%% tagexpr2tags: automatically find filenames

tags = tagexpr2tags(texpr)

tags(5)

 
 
 
%%

texpr = './Test/Images/hESCells_Folder/f<tag1,1>_t<tag2,1>/f<tag1,1>_t<tag2,2>_z<tag3,2>.tif';
[tnames, tagsplit, taginfo] = tagexpr2tagnames(texpr)

%% - type mismatch error

clc
texpr = './Test/Images/hESCells_Folder/f<tag1,d,5>_t<tag2,1>/f<tag1,s,5>_t<tag2,2>_z<tag3,2>.tif';
 
[tnames, tagsplit, taginfo] = tagexpr2tagnames(texpr)


%% - size mismatch error
clc
texpr = './Test/Images/hESCells_Folder/f<tag1,s,5>_t<tag2,1>/f<tag1,s,3>_t<tag2,2>_z<tag3,2>.tif';
 
[tnames, tagsplit, taginfo] = tagexpr2tagnames(texpr)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% tags2tagranges


tags = tagexpr2tags('./Test/Images/hESCells_Folder/f<tag1>_t<tag2,1>/f<tag1>_t<tag2,2>_z<tag3,2>.tif')

trs = tags2tagranges(tags)

tgs = tagranges2tags(trs)


%% tags not multipicative -> should produce error

tags2tagranges(tags, 'check', true)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test some regular expressions
texpr = 'test_T<tag1>-<tag2, 5>.tif';

tt = regexp(texpr, '<(?<name>\w*)\s*(?<k>(,\s*\d*?|\s*?))\s*>', 'names')
{tt.name}
{tt.k}

%%


tt = regexp(texpr, '(?<repl><.*?,*\s*?\d*?\s*?>)', 'names')
tt.repl


%% check for consistency

clc
texpr = './Test/Images/hESCells_Folder/f<tag1,1>_t<tag2,1>/f<tag1,1>_t<tag2,2>_z<tag3,2>.tif';

fname1 = './Test/Images/hESCells_Folder/f1_t1/f1_t01_z01.tif';
 
re = tagexpr2regexp(texpr)
regexp(fname1, re, 'names')


fname2 = './Test/Images/hESCells_Folder/f2_t1/f3_t02_z02.tif';
regexp(fname2, re, 'names')
 

tagexpr2tags(texpr, fname1)
tagexpr2tags(texpr, fname2)




%% TagRanges

clc
tr.tag1 = {1,2,3,4,5};
tr.tag2 = {'a', 'b', 'c'};

tsi = tagrangesize(tr)

prod(tsi)

tvs = ind2tagvalues(tr, 6)




 
