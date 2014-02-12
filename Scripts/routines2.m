    
%% segmentation Distance B

[labels_out,d]=IdentifySecPropagateSubfunction(PrelimPrimaryLabelMatrixImage,OrigImage,ThresholdedOrigImage,1.0);
            labels_out(d>DistanceToDilate) = 0;
            labels_out((PrelimPrimaryLabelMatrixImage > 0)) = PrelimPrimaryLabelMatrixImage((PrelimPrimaryLabelMatrixImage > 0));
            RelabeledDilatedPrelimSecObjectImage = labels_out;
            
            
                    %%% Removes objects that are not in the edited EditedPrimaryLabelMatrixImage.
        Map = sparse(1:numel(PrelimPrimaryLabelMatrixImage), PrelimPrimaryLabelMatrixImage(:)+1, EditedPrimaryLabelMatrixImage(:));
        LookUpColumn = full(max(Map,[], 1));
        LookUpColumn(1)=0;
        FinalLabelMatrixImage = LookUpColumn(RelabeledDilatedPrelimSecObjectImage+1);
        
        
        %% propagation
        
        
        
        %%% STEP 2: Starting from the identified primary objects, the secondary
        %%% objects are identified using the propagate function, written by Thouis
        %%% R. Jones. Calls the function
        %%% "IdentifySecPropagateSubfunction.mexmac" (or whichever version is
        %%% appropriate for the computer platform being used), which consists of C
        %%% code that has been compiled to run quickly within Matlab.

        % 2007-Jul-16 Kyungnam: If you want to get additional outputs, then
        % add more output arguments as follows:
        %%% [PropagatedImage, dist, diff_count, pop_count] = IdentifySecPropagateSubfunction(PrelimPrimaryLabelMatrixImage,OrigImage,ThresholdedOrigImage,RegularizationFactor);
        PropagatedImage = IdentifySecPropagateSubfunction(PrelimPrimaryLabelMatrixImage,OrigImage,ThresholdedOrigImage,RegularizationFactor);

        %%% STEP 3: We used the PrelimPrimaryLabelMatrixImage as the
        %%% source for primary objects, but that label-matrix is built
        %%% before small/large objects and objects touching the
        %%% boundary are removed.  We need to filter the label matrix
        %%% from propagate to make the labels match, and remove any
        %%% secondary objects that correspnd to size- or
        %%% boundary-filtered primaries.
        %%%
        %%% Map preliminary labels to edited labels based on maximum
        %%% overlap from prelim to edited.  We can probably assume
        %%% that no borders are adjusted during editing (i.e., any
        %%% changes from Prelim to Edited only involves removing
        %%% entire objects), but this is safer.
        %%%
        %%% (add one so that zeros are remapped correctly.)
        PrelimToEditedHist = sparse(EditedPrimaryLabelMatrixImage(:) + 1, PrelimPrimaryLabelMatrixImage(:) + 1, 1);
        [ignore, PrelimToEditedRemap] = sort(PrelimToEditedHist, 1);
        PrelimToEditedRemap = PrelimToEditedRemap(end, :) - 1;
        %%% make sure zeros map to zeros (note the off-by-one for the
        %%% index because Matlab doesn't do 0-indexing).
        PrelimToEditedRemap(1) = 0;
        EditedLabelMatrixImage = PrelimToEditedRemap(PropagatedImage + 1);

        %%% STEP 4:
        %%%
        %%% Fill holes (any contiguous, all-0 regions that are
        %%% surrounded by a single value).
        FinalLabelMatrixImage = CPfill_hole
        
        
        
        %% intensity intensity
        
        
        BlurredImage = conv2(OrigImage,f,'same') ./ conv2(ones(size(OrigImage)),f,'same');
        Objects = BlurredImage > Threshold;
        
        if strcmp(FillHolesOption,'Yes')
           Objects = imfill(double(Objects),'holes');                            % Fill holes
        end

        
        %%% Smooth images for maxima suppression
        if strcmpi(SizeOfSmoothingFilter,'Automatic')
           SizeOfSmoothingFilter=2.35*MinDiameter/3.5;
        end
        BlurredImage = CPsmooth(OrigImage,'Gaussian Filter',SizeOfSmoothingFilter,0);
        
        if strcmpi(MaximaSuppressionSize,'Automatic')
           MaximaSuppressionSize = round(MinDiameter/1.5);
        end
        MaximaMask = getnhood(strel('disk', MaximaSuppressionSize));
        
        %%% Initialize MaximaImage
        MaximaImage = ResizedBlurredImage;
        %%% Save only local maxima
        MaximaImage(ResizedBlurredImage < ...
           ordfilt2(ResizedBlurredImage,sum(MaximaMask(:)),MaximaMask)) = 0;
        
        
        if GetThreshold
           %%% Remove dim maxima
           %%% TODO: THIS IS THE MEAN THRESHOLD, SHOULDN'T IT
           %%% BE THE ORIG THRESHOLD?
           
           MaximaImage = MaximaImage > Threshold;
        end
        %%% Shrink to points (needed because of the resizing)
        if all(MaximaImage(:))
           MaximaImage = zeros(size(MaximaImage));
        else
           MaximaImage = bwmorph(MaximaImage,'shrink',inf);
        end
        
        
        if strcmp(WatershedTransformImageType,'Intensity')
           %%% Overlays the objects markers (maxima) on the inverted original image so
           %%% there are black dots on top of each dark object on a white background.
           Overlaid = imimposemin(1 - OrigImage,MaximaImage);
        end
        
        %%% Calculate the watershed transform and cut objects along the boundaries
        WatershedBoundaries = watershed(Overlaid) > 0;
        Objects = Objects.*WatershedBoundaries;
        
        %%% Label the objects
        Objects = bwlabel(Objects);
        
        %%% Remove objects with no marker in them (this happens occasionally)
        %%% This is a very fast way to get pixel indexes for the objects
        tmp = regionprops(Objects,'PixelIdxList');
        for k = 1:length(tmp)
           %%% If there is no maxima in these pixels, exclude object
           if sum(MaximaImage(tmp(k).PixelIdxList)) == 0
              Objects(tmp(k).PixelIdxList) = 0;
           end
        end
        
        %%% Merge small objects
        
        if strcmp(MergeChoice,'Yes')
           NumberOfObjectsBeforeMerge = max(Objects(:));
           Objects = MergeObjects(Objects,OrigImage,[MinDiameter MaxDiameter]);
           NumberOfObjectsAfterMerge = max(Objects(:));
           NumberOfMergedObjects = NumberOfObjectsBeforeMerge-NumberOfObjectsAfterMerge;
        end
        
        %%% Get diameters of objects and calculate the interval
        %%% that contains 90% of the objects
        tmp = regionprops(Objects,'EquivDiameter');
        Diameters = [0;cat(1,tmp.EquivDiameter)];
        SortedDiameters = sort(Diameters);
        NbrInTails = max(round(0.05*length(Diameters)),1);
        Lower90Limit = SortedDiameters(NbrInTails);
        Upper90Limit = SortedDiameters(end-NbrInTails+1);
        
        %%% Locate objects with diameter outside the specified range
        tmp = Objects;
        if strcmp(ExcludeSize,'Yes')
           %%% Create image with object intensity equal to the diameter
           DiameterMap = Diameters(Objects+1);
           %%% Remove objects that are too small
           Objects(DiameterMap < MinDiameter) = 0;
           %%% Will be stored to the handles structure
           SmallRemovedLabelMatrixImage = Objects;
           %%% Remove objects that are too big
           Objects(DiameterMap > MaxDiameter) = 0;
        else
           %%% Will be stored to the handles structure even if it's unedited.
           SmallRemovedLabelMatrixImage = Objects;
        end
        %%% Store objects that fall outside diameter range for display
        DiameterExcludedObjects = tmp - Objects;
        
        %%% Remove objects along the border of the image (depends on user input)
        tmp = Objects;
        if strcmp(ExcludeBorderObjects,'Yes')
           PrevObjects = Objects;
           Objects = CPclearborder(Objects);
           
           %%% CODE TO REMOVE BORDERS FROM ELLIPSE CROPPED OBJECTS
           if sum(PrevObjects(:)) == sum(Objects(:))
              try %#ok Ignore MLint
                 CropMask = CPretrieveimage(handles,['CropMask',ImageName],ModuleName);
                 CropBorders = bwperim(CropMask);
                 BorderTable = sortrows(unique([CropBorders(:) Objects(:)],'rows'),1);
                 for z = 1:size(BorderTable,1)
                    if BorderTable(z,1) ~= 0 && BorderTable(z,2) ~= 0
                       Objects(Objects == BorderTable(z,2)) = 0;
                    end
                 end
              end
           end
        end
        %%% Store objects that touch the border for display
        BorderObjects = tmp - Objects;
        
        %%% Relabel the objects
        [Objects,NumOfObjects] = bwlabel(Objects > 0);
        FinalLabelMatrixImage = Objects;
        
                  
                  