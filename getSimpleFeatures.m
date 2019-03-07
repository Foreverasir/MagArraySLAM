function [segmentCell, magSteps, magStepsIndex, accSteps, accStepsIndex] = getSimpleFeatures(OX_NUM, mag, varMag, varAcc, proThreshold, anchorThreshold, stepThreshold)
% Generate the Segment Infomation and the simple feature matrix.
% Input:
%   mag: the magnatic data.
%   varMag: the variance of the magnetic signal.
%   varAcc: the variance of the motion signal.
%   proThreshold: the threshold of prominence of Peaks.
%   anchorThreshold: the threshold of mag-step variance
%   stepThreshold: the threshold of motion-step variance or other
%
% Output:
%   segmentCell: A cell structure which contains the head & tail index for
%   each segment and the anchors index in each segment.
%       segmentCell(i,1): the head and tail index of the Mag-Array
%       detected --- [head, tail]
%       segmentCell(i,2): the step over the array info --- [head,tail,motion amount covered, Motion step bound, index;...]
%       segmentCell(i,3): the feature matrix of each segment.
%           Matrix: peaks, locations, width, prominence, ratio pro/width,
%           magStepIndex, accStepIndex, label

segmentCell = {};

% Mag Step detection 
j = 1;
magSteps = zeros(size(varMag));
magStepsIndex = zeros(ceil(size(varMag,1)/10), 6);
tempFlag = 0;
for i = 1:size(varMag,1)
    if varMag(i) > anchorThreshold
        magSteps(i) = 100;
        if tempFlag == 0
            tempFlag = 1;
            magStepsIndex(j,1) = i;
        end
    else
        if tempFlag == 1
            tempFlag = 0;
            magStepsIndex(j,2) = i;
            j = j+1;
        end
    end
end
magStepsIndex = magStepsIndex(1:j-1,:);
% Curation: Merge the close two 'step'
tempThreshold = 10;
if size(magStepsIndex,1)>1
    for i = 1:size(magStepsIndex,1)-1
        if magStepsIndex(i+1,1)-magStepsIndex(i,2) < tempThreshold
            magSteps(magStepsIndex(i,2):magStepsIndex(i+1,1)) = 100;
        end 
    end
end
j = 1;
tempFlag = 0;
for i = 1:size(magSteps,1)
    if magSteps(i) > 0
        if tempFlag == 0
            tempFlag = 1;
            magStepsIndex(j,1) = i;
        end
    else
        if tempFlag == 1
            tempFlag = 0;
            magStepsIndex(j,2) = i;
            j = j+1;
        end
    end
end
magStepsIndex = magStepsIndex(1:j-1,:);

% Step detection by motion.
j = 1;
accSteps = zeros(size(varAcc,1),1);
accStepsIndex = zeros(ceil(size(varAcc,1)/10),2);
tempFlag = 0;
varAccSum = sum(varAcc,2);
for i = 1:size(varAcc,1)
    if varAccSum(i) > stepThreshold
        accSteps(i) = 10;
         if tempFlag == 0
            tempFlag = 1;
            accStepsIndex(j,1) = i;
        end
    else
        if tempFlag == 1
            tempFlag = 0;
            accStepsIndex(j,2) = i;
            j = j+1;
        end
    end
end
accStepsIndex= accStepsIndex(1:j-1,:);

% Curation the range of step.
for i = 1:size(accStepsIndex,1)-1
    temp = accStepsIndex(i,2)-accStepsIndex(i,1);
    if temp > 24
        temp = log2(temp);
    end
    if accStepsIndex(i,2)+temp > accStepsIndex(i+1,1)
        accSteps(accStepsIndex(i,2):accStepsIndex(i+1,1))=10;
        continue;
    end
    temp = accStepsIndex(i+1,2)-accStepsIndex(i+1,1);
    if temp > 24
        temp = log2(temp);
    end
    if  accStepsIndex(i+1,1) - temp < accStepsIndex(i,2)
        accSteps(accStepsIndex(i,2):accStepsIndex(i+1,1))=10;
    end
end
tempFlag = 0;
j = 1;
accHelpMatrix = [];
for i = 1:size(accSteps,1)
    if accSteps(i) > 0
        if tempFlag == 0
            tempFlag = 1;
            accStepsIndex(j,1) = i;
        end
    else
        if tempFlag == 1
            tempFlag = 0;
            accStepsIndex(j,2) = i;
            temp = zeros(1,size(accSteps,1));
            temp(accStepsIndex(j,1):accStepsIndex(j,2))=1;
            accHelpMatrix=[accHelpMatrix;temp];
            j = j+1;
        end
    end
end
accStepsIndex= accStepsIndex(1:j-1,:);

%% Check each segment and remove the slope. Meanwhile set the feature matrix.
[highpeaksY, highlocsY, highWidthY, highProminenceY] = findpeaks(mag,'MinPeakProminence',proThreshold,'Annotate','extents');
[lowpeaksY, lowlocsY, lowWidthY, lowProminenceY]= findpeaks(-mag,'MinPeakProminence',proThreshold,'Annotate','extents');
ratioForHighPeaksY = highProminenceY ./ highWidthY;
ratioForLowPeaksY = lowProminenceY ./ lowWidthY;
OXMatrix = sortrows([[highpeaksY, highlocsY, highWidthY, highProminenceY, ratioForHighPeaksY, zeros(size(highlocsY)), zeros(size(highlocsY)), zeros(size(highlocsY))]; [lowpeaksY, lowlocsY, lowWidthY, lowProminenceY, ratioForLowPeaksY, zeros(size(lowlocsY)), zeros(size(lowlocsY)), ones(size(lowlocsY))]],2);

% OXHelp: which step contains the peak.
OXHelp = zeros(size(OXMatrix,1),2);
for i = 1:size(OXMatrix,1)
    for j = 1:size(magStepsIndex,1)
        if magStepsIndex(j,1) < OXMatrix(i,2) && OXMatrix(i,2) < magStepsIndex(j,2)
            OXHelp(i,1) = j;
            break;
        end
    end
    for j = 1:size(accStepsIndex,1)
        if accStepsIndex(j,1) < OXMatrix(i,2) && OXMatrix(i,2) < accStepsIndex(j,2)
            OXHelp(i,2) = j;
            break;
        end
    end
end

% Remove the peaks not in the mag-step and move the mag-step which contains no peaks.
helpIndex = [];
for i = 1:size(OXHelp,1)
    if OXHelp(i,1) == 0
        helpIndex = [helpIndex,i];
    end
end
OXMatrix(helpIndex,:) = [];
OXHelp(helpIndex,:) = [];

helpIndex = [];
for i = 1:size(magStepsIndex,1)
    if find(OXHelp(:,1)==i)
        continue;
    else
        helpIndex = [helpIndex,i];
        magSteps(magStepsIndex(i,1):magStepsIndex(i,2)) = 0;
    end
end
magStepsIndex(helpIndex,:) = [];

%% Bind the steps to one segement.
% bind the mag-step to the ox
for i = 1:size(OXMatrix,1)
    for j = 1:size(magStepsIndex,1)
        if magStepsIndex(j,1) < OXMatrix(i,2) && OXMatrix(i,2) < magStepsIndex(j,2)
            OXMatrix(i,6) = j;
            break;
        end
    end
    for j = 1:size(accStepsIndex,1)
        if accStepsIndex(j,1) < OXMatrix(i,2) && OXMatrix(i,2) < accStepsIndex(j,2)
            OXMatrix(i,7) = j;
            break;
        end
    end
end

% caculate the cover ratio
magHelpMatrix = zeros(size(magSteps,1), size(magStepsIndex,1));
for i = 1:size(magStepsIndex,1)
    for j = 1:size(accStepsIndex,1)
        if j == 1 
            if magStepsIndex(i,1) < accStepsIndex(j,2)
                t1 = j;
            end
        elseif magStepsIndex(i,1) > accStepsIndex(j-1,2) && magStepsIndex(i,1) < accStepsIndex(j,2)
            t1 = j;
        end
        
        if j == size(accStepsIndex,1)
            if magStepsIndex(i,2) > accStepsIndex(j,1)
                t2 = j;
            end
        elseif magStepsIndex(i,2) > accStepsIndex(j,1) && magStepsIndex(i,2) < accStepsIndex(j+1,1)
            t2 = j;
        end
    end
    
    magHelpMatrix(magStepsIndex(i,1):magStepsIndex(i,2),i) = 1;
    
    % coincidence amount
    [coincidence, temp] = max(accHelpMatrix([t1,t2],:) * magHelpMatrix(:,i));
        
    if t1 == t2
        magStepsIndex(i,4) = t1;
        magStepsIndex(i,6) = 0; % step-on array flag 
    else

        if temp == 1
            magStepsIndex(i,4) = t1;
        else
            magStepsIndex(i,4) = t2;
        end
        ratio = coincidence / (magStepsIndex(i,2)-magStepsIndex(i,1));
        if ratio < 0.512
            magStepsIndex(i,6) = 1;
        end
    end
    
    magStepsIndex(i,3) = coincidence;
    if magStepsIndex(i,3) < 0
        magStepsIndex(i,3) = 0;
    end
%     magStepsIndex(i,3) = sum(accSteps(magStepsIndex(i,1):magStepsIndex(i,2))./10);
     
end
% bind the acc step to the magstep
% [~,temp] = max(accHelpMatrix * magHelpMatrix);
% magStepsIndex(:,4) = temp';


%% SegmentCell Generate
j = 1;
% Bind magStepIndex to each segment, this method only can be used for two steps which means the LIMIT should be 4.
magStepsIndex(1,5) = 1;
for i= 1:size(magStepsIndex,1)-1
    magStepsIndex(i,5) = i;
    magStepsIndex(i+1,5) = i+1;
    if magStepsIndex(i,4)+1==magStepsIndex(i+1,4)
        segmentCell(j,1) = {[magStepsIndex(i,1),magStepsIndex(i+1,2)]};
        segmentCell(j,2) = {[magStepsIndex(i,:);magStepsIndex(i+1,:)]};
        j = j+1;
    elseif ( i==1 || (magStepsIndex(i,4)-1 ~= magStepsIndex(i-1,4)) ) && (magStepsIndex(i,4)+1 ~= magStepsIndex(i+1,4))
        segmentCell(j,1) = {[magStepsIndex(i,1),magStepsIndex(i,2)]};
        segmentCell(j,2) = {magStepsIndex(i,:)};
        j = j+1;
    end
end
if size(magStepsIndex,1)==1 || ((size(magStepsIndex,1)>1) && (magStepsIndex(end,4)-1 ~= magStepsIndex(end-1,4)))
    segmentCell(j, 1) = {[magStepsIndex(end,1),magStepsIndex(end,2)]};
    segmentCell(j, 2) = {magStepsIndex(end,:)};
end

magSteps = magSteps .* 0.1;
% Set related OX to the segment
for i = 1:size(segmentCell,1)
    oxInfo = [];
    magsi = segmentCell{i,2};
    for j =1:size(magsi,1)
        currentOX = OXMatrix(OXMatrix(:,6)==magsi(j,5),:);
        % Redo the oxInfo according to the magStepIndex.
        [hpk, hl, hw, hp] = findpeaks(mag(magsi(j,1)-proThreshold:magsi(j,2)+proThreshold),'MinPeakProminence',proThreshold,'Annotate','extents');
        [lpk, ll, lw, lp]= findpeaks(-mag(magsi(j,1)-proThreshold:magsi(j,2)+proThreshold),'MinPeakProminence',proThreshold,'Annotate','extents');
        rhp = hp ./ hw;
        rlp = lp ./ lw;
        tOX = sortrows([[hpk, hl+magsi(j,1)-proThreshold-1, hw, hp, rhp, zeros(size(hl)), zeros(size(hl)), zeros(size(hl))]; [lpk, ll+magsi(j,1)-proThreshold-1, lw, lp, rlp, zeros(size(ll)), zeros(size(ll)), ones(size(ll))]],2);
        for k = 1:size(tOX,1)
            t = find(currentOX(:,2) == tOX(k,2));
            if t > 0
                currentOX(t,3:5) = tOX(k,3:5);
            end
        end        
        
        if currentOX(1,5) <= currentOX(end,5)
            if currentOX(1,5) < 0.75 
                currentOX(1,:) = [];
            end
        else
            if currentOX(end,5) < 0.75
                currentOX = currentOX(1:end-1,:);
                currentOX(end,:) = [];
            end
        end
            
        oxInfo = [oxInfo; currentOX];
    end
    segmentCell(i,3) = {oxInfo};
    
    oxInfo = handleUnnecessaryOX(oxInfo, OX_NUM, segmentCell{i,2});
    segmentCell(i,4) = {oxInfo};
end

end

function result = handleUnnecessaryOX(oxInfo, Limit, index)
    if size(oxInfo,1) <= Limit
        result = oxInfo;
    else
        switch size(index,1)
            case 1
                [~,b] = min(oxInfo(:,5));
                if b == 1 || b == size(oxInfo,1)
                    oxInfo(b,:) = [];
                    result = handleUnnecessaryOX(oxInfo,Limit,index);
                else
                    result = oxInfo;
                end
                
               
            case 2
                ox1 = oxInfo(oxInfo(:,6) == index(1,5),:);
                ox2 = oxInfo(oxInfo(:,6) == index(2,5),:);
                [a1,b1] = min(ox1(:,5));
                [a2,b2] = min(ox2(:,5));
                if a1<=a2 && (b1 == 1 || b1 == size(ox1,1))
                    ox1(b1,:) = [];
                    result = handleUnnecessaryOX([ox1;ox2], Limit, index);
                elseif a1 > a2 && (b2 == 1 || b2 == size(ox2,1))
                    ox2(b2,:) = [];
                    result = handleUnnecessaryOX([ox1;ox2], Limit, index);
                else
                    result = oxInfo;
                end
            otherwise
                % TBD
        end
    end
end
