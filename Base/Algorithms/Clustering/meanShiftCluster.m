function [clpos, clids, cls] = meanShiftCluster(pts, ksize)
%
% [clpos,clids,cls] = meanShiftCluster(pts, ksize)
%
% description:
%     menashift clustering on pts using using a flat kernel of size ksize
%
% input:
%     pts    array of points as row vectors
%     ksize  width of flat kernel
%
% output:
%     clpos  array of final clustre positions as row vectors
%     clids  cluster ids of the data points
%     cls    ids of data points in the clusters as cell array
%
% note: 
%     based on K. Funkunaga and L.D. Hosteler, "The Estimation of the Gradient of a  Density Function, with Applications in Pattern Recognition"

np        = size(pts,2);
nc        = 0;
k2        = ksize^2;
ipids     = 1:np;
stopThresh= 1e-3*ksize;                           
clpos     = [];                            

beenhere     = zeros(1,np,'uint8');       
ninitp       = np;                          %number of points to posibaly use as initilization points
clusterVotes = zeros(1,np,'uint16');        %used to resolve conflicts on cluster membership


while ninitp

    tempInd      = ceil( (ninitp-1e-6)*rand);   % pick a random seed point
    stInd        = ipids(tempInd);              % use this point as start of mean
    mMean        = pts(:,stInd);                % intilize mean to this points location
    mMembers     = [];                          % points that will get added to this cluster                          
    newClusterVotes = zeros(1,np,'uint16');        % used to resolve conflicts on cluster membership

    
    while 1     % loop untill convergence
        

        sqDistToAll = sum((repmat(mMean,1,np) - pts).^2, 1);   %dist squared from mean to all points still active
        inInds      = find(sqDistToAll < k2);                  %points within bandWidth
        newClusterVotes(inInds) = newClusterVotes(inInds)+1;         %add a vote for all the in points belonging to this cluster
        
        
        oldMean   = mMean;                                  %save the old mean
        mMean     = mean(pts(:,inInds),2);                  %compute the new mean
        mMembers  = [mMembers, inInds];                     %#ok<AGROW> %add any point within bandWidth to the cluster
        beenhere(mMembers) = 1;                             %mark that these points have been visited


        % if mean doesn't move much stop this cluster
        if norm(mMean-oldMean) < stopThresh
            
            %check for merge posibilities
            mergeWith = 0;
            for cN = 1:nc
                distToOther = norm(mMean-clpos(:,cN));     %distance from posible new clust max to old clust max
                if distToOther < ksize/2                    %if its within bandwidth/2 merge new and old
                    mergeWith = cN;
                    break;
                end
            end
            
            
            if mergeWith > 0    
               % something to merge
                clpos(:,mergeWith)       = 0.5*(mMean+clpos(:,mergeWith));      %#ok<AGROW> %record the max as the mean of the two merged (biased twoards new ones)

                clusterVotes(mergeWith,:)= clusterVotes(mergeWith,:) + newClusterVotes;    %add these votes to the merged cluster
            else
               %its a new cluster
                nc                = nc+1;        %increment clusters
                clpos(:,nc)       = mMean;      %#ok<AGROW> %record the mean  

                clusterVotes(nc,:)= newClusterVotes;
            end

            break;
        end

    end

    ipids      = find(beenhere == 0);          % we can initialize with any of the points not yet visited
    ninitp     = length(ipids);                % number of active points in set

end


clusterVotes
[~,clids] = max(clusterVotes,[],1);             % a point belongs to the cluster with the most votes

% cluster cells
if nargout > 2
    cls = cell(nc,1);
    for cN = 1:nc
        mMembers = find(clids == cN);
        cls{cN} = mMembers;
    end
end