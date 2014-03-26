% TODO !!!
%
%  Compute Relative Intensitites for Surfaces
%
%  Copyright Christoph Kirst
%
%
%  Installation:
%
%  - Copy this file into the XTensions folder in the Imaris installation directory
%  - You will find this function in the Image Processing menu
%
%    <CustomTools>
%      <Menu>
%        <Item name="Add Relative Intensity Statistics" icon="Matlab" tooltip="Compute relative mean Intensities.">
%          <Command>MatlabXT::XTRelativeIntensities(%i)</Command>
%        </Item>
%      </Menu>
%      <SurpassTab>
%        <SurpassComponent name="bpSurfaces">
%           <Item name="Add Relative Intensity Statistics" icon="Matlab" tooltip="Compute relative mean Intensities.">
%               <Command>MatlabXT::XTRelativeIntensities(%i)</Command>
%           </Item>
%        </SurpassComponent>
%      </SurpassTab>
%    </CustomTools>
% 
%
%  Description:
%   
%   Computes relative mean Intensities.
%
%
function XTRelativeIntensities(aImarisApplicationID)

% connect to Imaris Ice interface
if ~isa(aImarisApplicationID, 'Imaris.IApplicationPrxHelper')
  javaaddpath ImarisLib.jar
  vImarisLib = ImarisLib;
  if ischar(aImarisApplicationID)
    aImarisApplicationID = round(str2double(aImarisApplicationID));
  end
  vImarisApplication = vImarisLib.GetApplication(aImarisApplicationID);
else
  vImarisApplication = aImarisApplicationID;
end

% Get Surpass Surfaces Object
vSurpassComponent = vImarisApplication.GetSurpassSelection;    

%Get the selected component
if vImarisApplication.GetFactory.IsSurfaces(vSurpassComponent) 
  %Get Surface Component
  vImarisObject     = vImarisApplication.GetFactory.ToSurfaces(vSurpassComponent);
%elseif vImarisApplication.GetFactory.IsCells(vSurpassComponent)
  %Get Cell Component
%  vImarisObject       = vImarisApplication.GetFactory.ToCells(vSurpassComponent);
else
  msgbox('Your selection is not a valid surfaces object');
  return;
end

%Check if there is a selection
if isempty(vImarisObject)
  msgbox('A surfaces object must be selected');
  return;
end


%Get the selected objects IDs
vIndicesOfSelectedObjects = vImarisObject.GetSelectedIds;

%Check if something is selected
if isempty(vIndicesOfSelectedObjects)
    warndlg('Please select one or more reference objects for the similarity computation');
    return
end

%Get the Name of the Similarity Statistics
vPrompt      = 'Enter the name of the similarity statistics (Similarity prefix will be added)';
vDlgTitle    = 'Statistics Name';
vNumOfLines  = 1;
vDefaultName = {''};
vSimilarityStatisticsName = char(inputdlg(vPrompt,vDlgTitle,vNumOfLines,vDefaultName));
if isempty(vSimilarityStatisticsName), return, end

%Remove Similarity Statistics
vImarisObject.RemoveStatistics(['Similarity' ' ' vSimilarityStatisticsName]);
 
%Get all Statistics
vProgressBar = waitbar(0,'Getting Statistics ...');
try
    vAllStatistics = vImarisObject.GetStatistics;
    vNames       = cell(vAllStatistics.mNames);
    vValues      = vAllStatistics.mValues;
%     vUnits       = cell(vAllStatistics.mUnits); % not used
    vFactors     = cell(vAllStatistics.mFactors);
    vFactorNames = cellstr(char(vAllStatistics.mFactorNames));
    vIds         = vAllStatistics.mIds;
catch er
    close(vProgressBar);
    error('Error in getting statistics');
end
    
%Get categories of the selected objects
vCategoryList     = {};
for vIndex = 1:size(vIndicesOfSelectedObjects, 1)
    j = find(vIds == vIndicesOfSelectedObjects(vIndex), 1);
    if (~isempty(j))
        if (~ismember(vFactors{1, j}, vCategoryList))
            vCategoryList = [vCategoryList(:, :); vFactors{1, j}];
        end
    end
end


if numel(vCategoryList) > 1
    %Show Category Dialog
    [vCategorySelection,vCategoryOk] = listdlg('ListString', vCategoryList, ... 
                                               'SelectionMode','Single',...
                                               'ListSize',[400 200],...
                                               'Name','Category - Selection Box',...
                                               'PromptString',{'Please select the category you want to compare:'});
    if vCategoryOk < 1
        close(vProgressBar);
        return
    end
else
    vCategorySelection = 1;
end

%Get the selected category
vCategory = vCategoryList{vCategorySelection};

%Get selected objects of the selected category
vObjectsSelectedByCategory = [];
for vIndex = 1:size(vIndicesOfSelectedObjects, 1)
    j = find(vIds == vIndicesOfSelectedObjects(vIndex), 1);
    if (~isempty(j))
        if (strcmp(vFactors{1, j},vCategory) == 1)
            vObjectsSelectedByCategory = [vObjectsSelectedByCategory(:, :); vIndicesOfSelectedObjects(vIndex)];
        end
    end
end

vIndicesOfSelectedObjects = vObjectsSelectedByCategory;

%Extract only selected object(s) related statistics
vRelevantStatistics = {};
for vIndex = 1:size(vNames, 1)
    if ismember(vIds(vIndex), vIndicesOfSelectedObjects)
        if isempty(strfind(vNames{vIndex}, 'Similarity'))
            vRelevantStatistics = [vRelevantStatistics(:, :); vNames(vIndex)];
        end
    end
end

vStatistics = unique(vRelevantStatistics);

%Show statistics selection dialog
[vSelection,vOk] = listdlg('ListString',   vStatistics, ... 
                           'SelectionMode','Multiple',...
                           'ListSize',     [400 200],...
                           'Name',         'Compute Similarity - Selection Box',...
                           'PromptString', {'Please select the Statistics for the Similarity Computation:'});
if vOk < 1
    close(vProgressBar);
    return
end

%Get statistics selected by the user (that should be used for the similarity computation)
vSimilarityStatistics = vStatistics(vSelection);
for vIndex = 1:size(vSelection, 2)
    vSimilarityStatistics = {vSimilarityStatistics{:}, vStatistics{vSelection(vIndex)}};
end

vSelectedObjectsMatrix       = [];
vObjectsMatrix               = [];

%Create matrix containing statistics values (every column contains all the values of the specific statistics)
try
    waitbar(0.2, vProgressBar, 'Creating Statistics Matrix ...');
    for j = 1:size(vSimilarityStatistics,2) 
        vFoundStatistic = strmatch(vSimilarityStatistics(j), vNames);  
        vIndexOfSelectedObjects = 1;
        vIndexOfAllObjects = 1;
        for k = min(vFoundStatistic) : max(vFoundStatistic)
            vObjectsMatrix(vIndexOfAllObjects, j) = vValues(k); 
            vIndexOfAllObjects = vIndexOfAllObjects + 1;
            if (ismember(vIds(k), vIndicesOfSelectedObjects))
                vSelectedObjectsMatrix(vIndexOfSelectedObjects, j) = vValues(k);
                vIndexOfSelectedObjects = vIndexOfSelectedObjects+1;
            end
        end
    end
catch er
    close(vProgressBar);
    error('Error in creating statistics matrix');
end

%Scale matrices values into the range [0, 1]
for vIndex = 1:size(vObjectsMatrix, 2)
    vMin = min(vObjectsMatrix(:, vIndex));
    vMax = max(vObjectsMatrix(:, vIndex));
    vObjectsMatrix(:, vIndex) = (vObjectsMatrix(:, vIndex) - vMin)/(vMax - vMin);
    vSelectedObjectsMatrix(:, vIndex) = (vSelectedObjectsMatrix(:, vIndex) - vMin)/(vMax - vMin);
end


try
    %Compute the centroid C of the selected objects
    waitbar(0.5, vProgressBar, 'Computing Centroid ...');
    %[IDX,C] = kmeans(vSelectedObjectsMatrix,1);
    C = 1:size(vSelectedObjectsMatrix,2);
    for vIndex = 1:size(vSelectedObjectsMatrix,2)
        C(vIndex) = mean(vSelectedObjectsMatrix(:,vIndex), 1);
    end

    vTransposeOfAllObjectsMatrix = transpose(vObjectsMatrix);

    %Compute the similarity matrix
    waitbar(0.7, vProgressBar, 'Computing Similarity Matrix ...');
    vCentroidSimilarity  = 1:size(vObjectsMatrix, 1);
    for vIndex = 1:size(vObjectsMatrix, 1)
        vCentroidSimilarity(vIndex)  = norm(transpose(C)-vTransposeOfAllObjectsMatrix(:,vIndex));
    end
catch er
    close(vProgressBar);
    error('Error in computing the similarity matrix');
end


%Add new statistics: Similarity
try
    waitbar(0.8, vProgressBar, 'Adding Similarity Statistics ...');
    
    %Remove factor 'Collection'
    vCollectionIndices = strmatch('Collection', vFactorNames);
    vFactorNames(vCollectionIndices) = [];
    vFactors(vCollectionIndices, :) = [];
    
    vFoundStatistic = strmatch(vSimilarityStatistics(1), vNames);
    vSize = numel(vFoundStatistic);
    aSimilarityNames   = cell(vSize, 1); 
    aSimilarityValues  = 1:vSize;
    aSimilarityUnits   = cell(vSize, 1);
    aSimilarityFactors = cell(size(vFactors, 1), vSize);
    aSimilarityIds     = 1:vSize;
    for j = 1:vSize
       k = vFoundStatistic(j);
       aSimilarityNames{j}      = ['Similarity', ' ', vSimilarityStatisticsName]; 
       aSimilarityValues(j)     = vCentroidSimilarity(j); 
       aSimilarityUnits{j}      = 'um'; 
       aSimilarityFactors(:, j) = vFactors(:, k); 
       aSimilarityIds(j)        = vIds(k);
    end
    vImarisObject.AddStatistics(aSimilarityNames, aSimilarityValues, ...
      aSimilarityUnits, aSimilarityFactors, vFactorNames, aSimilarityIds);
catch er
    close(vProgressBar);
    error('Error in adding similarity statistics');
end

waitbar(1.0, vProgressBar, 'Finishing ...');
vHelpMessage = msgbox(['"Similarity' ' ' vSimilarityStatisticsName '"' ' ' 'statistics value has been added. It can be used in Filter, Statistics and Color coding tabs.'], ...
                      'Information', 'help')
close(vProgressBar);



