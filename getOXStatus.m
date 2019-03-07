function [oxStatus] = getOXStatus(segmentCell, coverRatio, velocity, magStepsIndex)
% oxStatus:each cell corresponds to the segment. 
%   {:, 1} = motion-step covering ratio, the mean of around velocity
%   {:, 2} = ratio for Prominence. n * n+1, the last column means the sum
%   {:, 3} = ratio for Width
%   {:, 4} = metric matrix
%   {:, 5} = twins Index1, index2, ratio Pro, ratio Width, mag step Index.

for i = 1:size(segmentCell,1)
    oxInfo = segmentCell{i,3};
    current = zeros(size(oxInfo,1),2);
    current(:,1) = coverRatio(oxInfo(:,6));
    t1 = [];
    t2 = [];
    oxMap = mapminmax(oxInfo(:,2:6)',0,1)';
    oxMap(:,5) = oxMap(:,5)*2;
    metric = squareform(pdist(oxMap));
    metric(metric == 0) = 10000; 

    for j = 1:size(oxInfo,1)
        % Only take the x-v now.
        current(j,2) = mean(velocity((oxInfo(j,2)-2:oxInfo(j,2)+2),1));
        t1 = [t1, oxInfo(:,4) ./ oxInfo(j,4)];
        t2 = [t2, oxInfo(:,3) ./ oxInfo(j,3)];
        temp = 1000;
        for k = 1:size(oxInfo,1)
            if oxInfo(k,end) ~= oxInfo(j,end) && metric(j,k) < temp
                temp = metric(j,k);
                metric(j, size(oxInfo,1)+1) = k;
            end
        end
    end
    t1 = [t1,zeros(size(t1',1),1)];
    t2 = [t2,zeros(size(t2',1),1)];
    t1(:,end) = (sum(t1,2)-1) ./ (size(t1,2)-2);
    t2(:,end) = (sum(t2,2)-1) ./ (size(t2,2)-2);
    
    % Generate the twins index. No check now.
    temp = zeros(size(oxInfo,1),1);
    twins = [];
    if size(metric,1) < size(metric,2) 
    for j = 1:size(oxInfo,1)-1
        if metric(metric(j,end),end) == j && temp(j) == 0 && temp(metric(j,end)) == 0 && oxInfo(j,6)==oxInfo(metric(j,end),6)
            if magStepsIndex(oxInfo(j,6),6) == 0
                twins = [twins;j, metric(j,end), t1(j,metric(j,end)), t2(j,metric(j,end)), oxInfo(j,6)];
                temp(j) = 1;
                temp(metric(j,end)) = 1;
            end
        end
    end
    end
    
    oxStatus(i,1) = {current};
    oxStatus(i,2) = {t1};
    oxStatus(i,3) = {t2};
    oxStatus(i,4) = {metric};
    oxStatus(i,5) = {twins};
end
end

