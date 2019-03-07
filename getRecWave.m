function recWave = getRecWave(mag, segmentCell, magSteps)
% 
recWave = {};

for i = 1:size(segmentCell,1)
    stepInfo = segmentCell{i,2};
    oxInfo = segmentCell{i,3};
    range = stepInfo(:,1:2);
    currentRec = {};
    for j = 1:size(range,1)
        oxCurrent = oxInfo(oxInfo(:,6) == stepInfo(j,5),:);
        cmag = mag(range(j,1):range(j,2));
        
        a = mean(mag(range(j,1)-5:range(j,1)));
        b = mean(mag(range(j,2):range(j,2)+5));
        startmean = zeros(5,1) + a;
        endmean = zeros(5,1) + b;
        
        
        wave = zeros(range(j,2)-range(j,1)+1, 1);
        cutpoint = zeros(size(oxCurrent,1)+1,1);
        cutpoint(1) = 1;
        for k = 2:size(oxCurrent,1)
            cutpoint(k) = ceil((oxCurrent(k,2)+oxCurrent(k-1,2))/2 - range(j,1) + 1);
        end
        cutpoint(end) = range(j,2) - range(j,1) + 1;
        for k = 2:size(cutpoint,1)
            if oxCurrent(k-1,end) == 1
                wave(cutpoint(k-1):cutpoint(k)) = -oxCurrent(k-1,1);
            else
                wave(cutpoint(k-1):cutpoint(k)) = oxCurrent(k-1,1);
            end
        end
        wave = [startmean;wave;endmean];
        wave = wave - (a+b)/2;
        radio = (max(wave)-min(wave))/20;
        wave = wave/radio;
        
        currentRec(j) = {wave};
%         figure
%         subplot(2,1,1)
%         set(gcf,'color','w') 
%         plot(cmag); 
%         grid on;
%         title('cMag');
% %         text(cutpoint',zeros(size(cutpoint',1),1),'|','color','b');
% %         text(oxCurrent(:,2)-range(j,1),zeros(size(oxCurrent,1),1),'|','color','r');
% 
%         subplot(2,1,2)
%         set(gcf,'color','w') 
%         plot(wave);
%         grid on; 
%         title('wave');
        
    end
    recWave(i) = {currentRec};
end
end

