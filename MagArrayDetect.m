% Main file for mag-array experiments
% Encapsulation for functions should be considered this time. 2018.12.22
clear;

OX_NUM_SET = 4;
% Flag defined here.
PigtureFlag = true;
MovingFilterFlag = true;

%% Read the csv file
% Attention we delete the first column in csv file.
rawData =  csvread('data11.3\lab1-9-1.csv', 1, 1);
% Remove the head and tail data
rawData = rawData(10:size(rawData(:, 1))-10, :);
dataSize = size(rawData);
% Get the data to use 
rawTimestamp = rawData(:, 1);
rawAcc = rawData(:, 2:4);
rawLinearAcc = rawData(:, 5:7);
rawGyro = rawData(:, 8:10);
rawGravity = rawData(:, 11:13);
rawMag = rawData(:, 14:16); 
% Projection  
unitGravity = rawGravity;
for i = 1 : dataSize(1)
    unitGravity(i,:) = unitGravity(i,:)./norm(rawGravity(i,:));
end
magY = sum(rawMag .* unitGravity, 2);
projectData = rawMag - magY .* unitGravity;

%% Preprocess
% moving mean filter
windowSize = 5;
magY = movingMeanFilt(magY, windowSize);    
%acc = movingMeanFilt(rawAcc, windowSize);
linearAcc = movingMeanFilt(rawLinearAcc,windowSize);

rawTimestamp = rawTimestamp(windowSize:dataSize(1), :);
dataSize(1) = size(rawTimestamp, 1);

% S-G filter
sgMagY = sgolayfilt(magY,10,21);
sgMagYDiff = diff(sgMagY);
sgMagYDiff2 = movingMeanFilt(sgMagYDiff, windowSize*2);
%sgAcc = sgolayfilt(acc,10,21);
sgLinearAcc = sgolayfilt(linearAcc,10,21);

% Variance 
varSgMagY = sgMagY;
varSgLAcc = sgLinearAcc;
for i = 1:dataSize(1)
    if i >= dataSize(1)-windowSize
        varSgMagY(i) = varSgMagY(i-1);
        varSgLAcc(i,:) = varSgLAcc(i-1,:);
    else
        varSgMagY(i) = var(sgMagY(i:i+windowSize));
        varSgLAcc(i,:) = var(sgLinearAcc(i:i+windowSize,:));
    end
end


%% Get the Mag-Data Segment Infomation and simple feature matrix.
proYLimit = 10;
% segmentCell: A cell structure which contains the head & tail index for
%   each segment and the anchors index in each segment.
%       segmentCell(i,1): the head and tail index of the Mag-Array
%       detected --- [head, tail]
%       segmentCell(i,2): the step over the array info --- [head,tail,motion amount covered, Motion step bound, index;]
%       segmentCell(i,3): the feature matrix of each segment.
%           Matrix: peaks, locations, width, prominence, ratio pro/width, magStepIndex, accStepIndex, label
[segmentCell, magSteps, magStepsIndex, accSteps, accStepsIndex] = getSimpleFeatures(OX_NUM_SET, sgMagY, varSgMagY, varSgLAcc, proYLimit, 20, 1);

%% Caculate the new rectangle-wavelet 2019.2.28
recWave= getRecWave(sgMagY, segmentCell, magSteps);


%% Analyse the step on the array and the velocity influence
[stepVelocity, velocity] = motionAnalyse(sgLinearAcc, accStepsIndex, rawTimestamp, true);
% oxStatus:each cell corresponds to the segment.  
%   {:,1} = motion-step covering ratio, the mean of around velocity
%   {:,2} = ratio for Prominence. n * n+1, the last column means the mean
%   {:,3} = ratio for Width
%   {:,4} = metric matrix
%   {:,5} = twins Index [twins1, twins2, ratio Pro, ratio Width, mag-step Index]

coverRatio = magStepsIndex(:,3) ./ (magStepsIndex(:, 2) - magStepsIndex(:,1) + 1);
oxStatus = getOXStatus(segmentCell, coverRatio, velocity, magStepsIndex);

dtwRes = dtwWave(recWave);

%% Advanced features by motion & Mag Waves
% Set the curation threshold here.
TWINS_PRO_BOUND = 3.25;
TWINS_WIDTH_BOUND = 2.75;

% each line advOX: info about motion[original, curationed], info about mag, seq of motion, seq of mag
for i = 1:size(segmentCell,1)
    oxInfo = segmentCell{i,3};
    oxV = oxStatus{i,1};
    oxProR = oxStatus{i,2};
    oxWidthR = oxStatus{i,3};
    oxTwins = oxStatus{i,5};
    
    % if the v == 0,may curse no sequence, so give the least speed.
    for j = 1:size(oxV,1)
        if oxV(j,2) == 0
            oxV(j,2) = 10000;
            oxV(j,2) = min(oxV(:,2))/2;
        end
    end
    
    % if the step over the array, give award.
    oxV(oxV(:,1) > 0.618,1) = 1;
    advOX{i,1} = [oxV(:,1) .* oxV(:,2) .* oxInfo(:,3), oxV(:,1) .* oxV(:,2) .* oxInfo(:,3)];
    advOX{i,2} = [oxProR(:,end) .* oxWidthR(:, end), oxProR(:,end) .* oxWidthR(:, end)];
    
    advOX{i,3} = {get01Sequence(oxInfo,advOX{i,1}(:,1))};
    advOX{i,4} = {get01Sequence(oxInfo,advOX{i,2}(:,1))};
    
    % Curation by the twins.
    for j = 1:size(oxTwins,1)
        curationLevel = 0;
        if oxTwins(j,3) < TWINS_PRO_BOUND && oxTwins(j,3) > 1/TWINS_PRO_BOUND
            curationLevel = curationLevel + 0.5;
        elseif oxTwins(j,3) <= 1/TWINS_PRO_BOUND
            curationLevel = curationLevel + 0.5 * oxTwins(j,3);
        else
            curationLevel = curationLevel + 0.5 * (1 / oxTwins(j,3));
        end
        if oxTwins(j,4) < TWINS_WIDTH_BOUND && oxTwins(j,4) > 1/TWINS_WIDTH_BOUND
            curationLevel = curationLevel + 0.5;
        elseif oxTwins(j,4) <= 1/TWINS_WIDTH_BOUND
            curationLevel = curationLevel + 0.5 * oxTwins(j,4);
        else
            curationLevel = curationLevel + 0.5 * (1 / oxTwins(j,4));
        end
        oxTwins(j,5) = curationLevel;
        
        delta1 = (advOX{i,1}(oxTwins(j,1),1) - advOX{i,1}(oxTwins(j,2),1)) / 2;
        delta2 = (advOX{i,2}(oxTwins(j,1),1) - advOX{i,2}(oxTwins(j,2),1)) / 2;
        advOX{i,1}(oxTwins(j,1),2) = advOX{i,1}(oxTwins(j,1),2) - delta1 * curationLevel;
        advOX{i,1}(oxTwins(j,2),2) = advOX{i,1}(oxTwins(j,2),2) + delta1 * curationLevel;
        advOX{i,2}(oxTwins(j,1),2) = advOX{i,2}(oxTwins(j,1),2) - delta2 * curationLevel;
        advOX{i,2}(oxTwins(j,2),2) = advOX{i,2}(oxTwins(j,2),2) + delta2 * curationLevel;
    end
    
    % Generate the new 01 sequence
    advOX{i,3} = [advOX{i,3};get01Sequence(oxInfo,advOX{i,1}(:,2))];
    advOX{i,4} = [advOX{i,4};get01Sequence(oxInfo,advOX{i,2}(:,2))];
end

%% Judgement from 01-Sequence
% the function has threshold RATIO_LIMIT = 3.25;
% each line result: motion result, acc result
result = judgeFrom01Sequence(advOX, oxStatus, OX_NUM_SET, magStepsIndex);

%% Show the figure.
if PigtureFlag
%     figure
%     set(gcf,'color','w')
%     subplot(2,1,1)
%     plot(rawTimestamp, velocity,'LineWidth',1.2); 
%     grid on;
%     legend('X','Y','Z');
%     xlabel('Time(s)');
%     ylabel('v(m/s)')
%     title('velocity from S-G Linear Acc');
%     
%     subplot(2,1,2);
%     plot(rawTimestamp, [sgMagY,magSteps],'LineWidth',1.2); 
%     grid on;
%     xlabel('Time(s)');
%     ylabel('Magnetic Field(uT)')
%     title('S-G magY');
    
%     figure
%     set(gcf,'color','w')
%     subplot(2,1,1)
%     plot(rawTimestamp, rawMag(windowSize:end,2),'LineWidth',1.5); 
%     grid on;
%     xlabel('Time(s)');
%     ylabel('Magnetic Field(uT)')
%     title('Raw Mag Data');
    
    figure
    set(gcf,'color','w')
    plot(rawTimestamp, [sgMagY,magSteps],'LineWidth',1.2); 
    grid on;
    xlabel('Time(s)');
    ylabel('Magnetic Field(uT)')
    title('S-G magY');
%     
%     figure
%     set(gcf,'color','w')
%     plot(rawTimestamp(2:end), [sgMagYDiff,magSteps(2:end)],'LineWidth',1.2); 
%     grid on;
%     xlabel('Time(s)');
%     ylabel('Magnetic Field(uT) Diff')
%     title('S-G magY');
%     
%     figure
%     set(gcf,'color','w')
%     plot(rawTimestamp(windowSize*2+1:end), sgMagYDiff2,'LineWidth',1.2); 
%     grid on;
%     xlabel('Time(s)');
%     ylabel('Magnetic Field(uT) Diff')
%     title('S-G magY');
    
%     
%     for i = 1:size(segmentCell,1)
%         figure
%         set(gcf,'color','w')
%         plot(rawTimestamp(segmentCell{i,1}(1)-10:segmentCell{i,1}(2)+10), [sgMagY(segmentCell{i,1}(1)-10:segmentCell{i,1}(2)+10),magSteps(segmentCell{i,1}(1)-10:segmentCell{i,1}(2)+10)]); 
%         grid on;
%         xlabel('Time(s)');
%         ylabel('Magnetic Field(uT)')
%         title('S-G magY');
%     end

end