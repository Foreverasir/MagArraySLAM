function [stepVelocity, v] = motionAnalyse(acc, accStepsIndex, time, opflag)
% TODO: Now we just use the raw linear acc, but we should consider the
% projection and the forward direction.
% Input:
%   acc: the acc data.
%   accStepIndex: the motion step info.
%   time: the timestamp
% Output:
%   stepVelocity: a cell to put different size data

stepVelocity = {};
v = zeros(size(acc));
for i = 1:size(accStepsIndex,1)
    velocity = zeros(accStepsIndex(i,2)-accStepsIndex(i,1)+1,3);
    y = acc(accStepsIndex(i,1):accStepsIndex(i,2),:);
    temp = (y(end,:)+y(1,:))./(accStepsIndex(i,2)-accStepsIndex(i,1)-1);
    y(2:end-1,:) = y(2:end-1,:) - temp;
    y(1,:) = 0.0;
    y(end,:) = 0.0;
    x = time(accStepsIndex(i,1):accStepsIndex(i,2));
    for j=2:size(velocity,1)
        velocity(j, :) = velocity(j-1,:)+(x(j)-x(j-1))*(y(j,:)+y(j-1,:))/2;
    end
    % zero-velocity curation. i*delta/(n-1)
    temp = [1:accStepsIndex(i,2)-accStepsIndex(i,1)]' * (velocity(end,:)./(accStepsIndex(i,2)-accStepsIndex(i,1)));
    velocity(2:end,:) = velocity(2:end,:) - temp;
    if opflag
        velocity(:,1) = -velocity(:,1);
        velocity(:,2) = -velocity(:,2);
    end
    [peaks,locs,width,pro] = findpeaks(velocity(:,1),'Annotate','extents');
    
    v(accStepsIndex(i,1):accStepsIndex(i,2),:) = velocity; 
    
    stepVelocity(i,1) = {velocity(:,1)};
    stepVelocity(i,2) = {[peaks,locs,width,pro]};
end
%     figure
%     subplot(2,1,1)
%     set(gcf,'color','w') 
%     plot(y); 
%     grid on;
%     legend('X','Y','Z');
%     title('Linear Acc within one step');
%     subplot(2,1,2)
%     set(gcf,'color','w') 
%     plot(velocity); 
%     grid on;
%     legend('X','Y','Z');
%     title('Velocity within one step');
end