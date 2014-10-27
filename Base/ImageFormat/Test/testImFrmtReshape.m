 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Test Image Formats    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear all
close all
initialize


%%
clc

data =1:(12*5);
data = reshape(data, [12, 5]);
datar = reshape(data, [3,4,5]);

datarr = imfrmtReshapeOnce(data, 'XZ', 'XYZ', 'X', 'XY', [3,4])

size(datar)
size(datarr)

max(max(max(abs(datar-datarr))))

%% with flipping
clc
datar = reshape(data(1:12)', [3,4]);
datarr = imfrmtReshapeOnce(data(1:12)', 'X', 'XY', 'X', 'Xy', [3,4])

max(max(max(abs(datar-datarr))))


%% sequential

data = rand(3,16,12);

clc
datar = imfrmtReshape(data, 'XYZ', 'yXZCT', {'Y', 'Z'}, {'YC', 'ZT'}, {[4,4], [3,4]});

size(data)
size(datar)

%%

data = rand(8,16);

clc
datar = imfrmtReshape(data, 'XY', 'YX', {}, {}, {});

figure(1); clf
implottiling({data,; datar}, 'link', false)



%% size

dsize = [10,20,30];

dsizer = imfrmtReshapeSize(dsize, 'XYZ', 'XyZC', 'Z', 'ZC', [10,3])


%% data to cell conversion

data = rand(3,4,5);

cdata = imfrmtDataToCellData(data, 2)
size(cdata)

datac = imfrmtCellDataToData(cdata);

size(datac)

%% range size

imfrmtRangeSize([10,20,30], 'XYZ', setParameter('Z', 1:3))

%% index to coords

rgs = imfrmtIndexToCoordinate([10,20,30], 'XYZ', 100:110)

%%

imfrmtCoordinateToIndex([10,20,30], 'XYZ', rgs)


%% 

idx = imfrmtRangeToIndex([10,20,30], 'XYZ', struct('Y', 1:3, 'Z', 4))

numel(idx)


%% reformat range

rge = struct('X', 1:4, 'y', 5:8)

rger = imfrmtReformatRange([20,30], 'XY', 'Y', rge)


%% reshape inverse size

si =[10,20,4,3];

rsi = imfrmtReshapeInverseSize(si, 'XYZT', 'XS', {'S', 'X'}, {'ZT', 'XY'}, {[4,3],[10,20]} )


%% reshape inverse range

tgr = struct('X', 1:5, 'Y', 1:3, 'Z', 1:10, 'C', 1:3)


[tgri,resh] = imfrmtReshapeInverseRange([10, 5, 20, 5], 'XYZC', 'UV', {'U', 'V'}, {'XY', 'ZC'}, {[10,5], [20,5]}, tgr)

var2char(tgri.V')
var2char(resh)


% need to know reshape sizes for the inverse ranges in order to obtain data corresponding to data ranges


%%
clc
tgr = struct('Z', 1)

[tgri,resh] = imfrmtReshapeInverseRange([1024, 1300, 20, 5], 'XYZC', 'UV', {'U', 'V'}, {'XY', 'ZC'}, {[1024,1300], [20,5]}, tgr)

tgri.V'
var2char({'reshape', resh})




%% data format

imfrmtReshapeFormat('UVW', 'V', 'XYZ')




%%

clc
tgr = struct('Z', 1);

[tgri,resh] = imfrmtReshapeInverseCellDataRange([512,512, 15], 'XYZ', 'XY', 12, 'U', 'VU', 'V', 'Z', 15, tgr)

var2char({'reshape', resh})



%% Reshape Cell Data
clc

cd = repmat({[1,2;3,4]}, 3,4);

dfrmt = 'XY'; cfrmt = 'CS';
rsf = {'S'};
rst = {'UV'};
rss = {[2,2]};

outdfrmt = 'XYC';
outcfrmt = 'UV';


outcd = imfrmtReshapeCellData(cd, dfrmt, cfrmt, outdfrmt, outcfrmt, rsf, rst, rss)


%%


clc

cd = repmat({[1,2;3,4]}, 4);
for i = 1:4
   cd{i} = cd{i} * i;
end

dfrmt = 'XY'; cfrmt = 'S';
rsf = {'S'};
rst = {'UV'};
rss = {[2,2]};

outdfrmt = 'Yx';
outcfrmt = 'Uv';


outcd = imfrmtReshapeCellData(cd, dfrmt, cfrmt, outdfrmt, outcfrmt, rsf, rst, rss)


var2char(outcd)






