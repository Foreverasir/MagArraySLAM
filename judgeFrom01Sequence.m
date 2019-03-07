function result = judgeFrom01Sequence(advOX, oxStatus, LIMITNUM, magStepsIndex)
% Generate the result.
% Input:
%   advOX: the advanced features cell.
%   oxStatus: you know that.
% Output:
%   result: 01 result.

if size(advOX,1) ~= size(oxStatus,1)
    return;
end

result = {};
magStepCount = 1;
for i = 1:size(advOX,1)
    oxNum = size(oxStatus{i,1},1);
    limitNum = LIMITNUM;
    
    seqAcc = advOX{i,3}{2};
    seqMag = advOX{i,4}{2};
    rA = [];
    rM = [];
    rABackup = [];
    rMBackup = [];
    rAa = '';
    rAb = '';
    
    segA = getStepSeg(seqAcc, '#');
    segM = getStepSeg(seqMag, '#');
    
    for j = 1:size(segA,2)
        [tA,tMetricA,sum_A,index_A,sum0A,index0A,sum1A,index1A] = analyseStepSeq(segA{j});
        [tM,tMetricM,sum_M,index_M,sum0M,index0M,sum1M,index1M] = analyseStepSeq(segM{j});
        
        if size(tA,2) ~= size(tM,2)
            return;
        end     
        
        twins = oxStatus{i,5};
        twin = [];
        if size(twins,1) > 0
            twin = twins(twins(:,5) == magStepCount,:);
        end
        
        flag = magStepsIndex(magStepCount,6);
        
        [rAa, rAb, f1] = OXJudge(segA{j}{1}, tA, tMetricA, sum0A, sum1A, twin, flag);
        [rMa, rMb, f2] = OXJudge(segM{j}{1}, tM, tMetricM, sum0M, sum1M, twin, flag);
        
        % 01 ->0011
        rAa = extand01(rAa,f1,segA,limitNum);
        rMa = extand01(rMa,f2,segM,limitNum);

        % Two results cooprate.
        
        rA = [rA, rAa];
        rABackup = [rABackup, rAb];
        rM = [rM, rMa];
        rMBackup = [rMBackup, rMb];
        
        magStepCount = magStepCount + 1;
    end    
    
    result{i,1} = {rA};
    result{i,2} = {rM};
    result{i,3} = {rABackup};
    result{i,4} = {rMBackup};
end

end

% Some fundmental functions 
function stepSeg = getStepSeg(seq, label)
    a = find(seq == label);
    for j = 1:size(a, 2)
        if j == 1
            stepSeg{j} = {seq(1 : a(1)-1)};
        else
            stepSeg{j} = {seq(a(j-1)+1 : a(j)-1)};
        end
    end 
end

function [t,tMetric,s,a,s0,a0,s1,a1] = analyseStepSeq(stepSeq)
    s = sum(stepSeq{1} == '_');
    a = find(stepSeq{1} == '_');
    s0 = sum(stepSeq{1} == '0');
    s1 = sum(stepSeq{1} == '1');
    a0 = find(stepSeq{1} == '0');
    a1 = find(stepSeq{1} == '1');
    t = [];
    tMetric = [];
    % get the num of 0 or 1
    for k = 1:size(a,2)
        if k == 1
            t = [t, a(k) - 1];
        else
            t = [t, a(k) - a(k-1) - 1];
        end
    end
    for k = 1:size(t,2)
        tMetric = [tMetric;abs(t - t(k))];
    end
end

function [r, rBackup, flag] = OXJudge(seq, t, metric, sum0, sum1, twin, stepOn)
    
    flag = size(t,2);
    
    RATIO_LIMIT = 3.25;
    a = '';
    b = '';
    switch size(t,2)
        case 1
            if sum0>0
                a = '00';
            elseif sum1>0
                a = '11';
            end
            r = a;
            rBackup = a;
            
        case 2
            ratio = sum0 / sum1;
            
            if ratio < RATIO_LIMIT && ratio > 1/RATIO_LIMIT
                if seq(1)=='0'
                    a = '01';
                    b = '01';
                else
                    a = '10';
                    b = '10';
                end
            elseif ratio >= RATIO_LIMIT
                if seq(1)=='0'
                    a = '0001';
                    b = '001';
                else
                    a = '1000';
                    b = '100';
                end
            elseif ratio <= 1/RATIO_LIMIT
                if seq(1)=='0'
                    a = '0111';
                    b = '011';
                else
                    a = '1110';
                    b = '110';
                end     
            end
            
            if size(twin,1) > 0 && twin(1,1)-twin(1,2) == -1
                if seq(1)=='0'
                    r = '01';
                else
                    r = '10';
                end
                rBackup = a;
            else                
                r = a;
                rBackup = b;
            end
            
            if stepOn > 0
                r = '**';
                if seq(1)=='0'
                    rBackup = '01';
                else
                    rBackup = '10';
                end
            end
            
        case 3
            s1 = seq(1);
            s2 = seq(t(1)+2);
            s3 = seq(t(1)+t(2)+3);
            
            % TODO: consider twin
            if size(twin,1) == 1 && twin(1,1)-twin(1,2) == -1
                if twin(1,1) == 1
                    if seq(1)=='0'
                        a = '0100';
                    else
                        a = '1011';
                    end
                    
                    if t(3) > t(2) && t(3) > t(1)
                        b = a;
                    elseif t(3) < t(2) && t(3) < t(1)
                        b = a; % The original inference set as backup
                        a = [s1,s2];
                    else
                        % twin is wrong 
                        b = [s1,s1,s2,s3];
                    end
                    
                elseif twin(1,1) == 2
                    if s2=='0'
                        a = '1101';
                    else
                        a = '0010';
                    end
                    
                    if t(1) > t(2) && t(1) > t(3)
                        b = a;
                    elseif t(1) < t(2) && t(1) < t(3)
                        b = a;
                        a = [s2,s3];
                    else
                        b = [s1,s2,s3,s3];
                    end
                end         
            else                
                if metric(2,1) < metric(2,3)
                    if seq(1)=='0'
                        a = '0100';
                    else
                        a = '1011';
                    end
                    if t(3) > t(2) && t(3) > t(1)
                        b = a;
                    elseif t(3) < t(2) && t(3) < t(1)
                        b = a; % The original inference set as backup
                        a = [s1,s2];
                    else
                        % twin is wrong 
                        b = [s1,s1,s2,s3];
                    end
                    
                elseif metric(2,1) > metric(2,3)
                    if seq(1)=='0'
                        a = '0010';
                    else
                        a = '1101';
                    end
                    if t(1) > t(2) && t(1) > t(3)
                        b = a;
                    elseif t(1) < t(2) && t(1) < t(3)
                        b = a;
                        a = [s2,s3];
                    else
                        b = [s1,s2,s3,s3];
                    end
                end
            end
            
            r = a;
            if stepOn > 0
                r = '***';
            end
            rBackup = b;
            
        case 4
            s1 = seq(1);
            s2 = seq(t(1)+2);
            s3 = seq(t(1)+t(2)+3);
            s4 = seq(t(1)+t(2)+t(3)+4);
            % Strong condition
            if size(twin,1) == 2 && twin(1,1)-twin(1,2) == -1 && twin(1,2) - twin(2,1) == -1 && twin(2,1) - twin(2,2) == -1
                a = [s1,s2,s3,s4];
                b = a;
                if t(2)>t(1)+t(4) && t(3)>t(1)+t(4)
                    b = [s2,s3];
                end
                if t(1)>t(2)+t(3) && t(4)>t(2)+t(3)
                    if t(1) > t(4)
                        b = [s1,s1,s1,s4];
                    elseif t(1) < t(4)
                        b = [s1,s4,s4,s4];
                    else
                        b = [s1,s4];
                    end
                end
            elseif size(twin,1) == 2 && twin(1,1) - twin(2,2) == -3 && twin(1,2) - twin(2,1) == -1
                if t(1) > t(4)
                    a = [s1,s1,s1,s4];  
                    b = [s1,s2,s3,s4];
                elseif t(1) < t(4)
                    a = [s1,s4,s4.s4];
                    b = [s1,s2,s3,s4];
                else
                    a = [s1,s2,s3,s4];
                    b = [s1,s2,s3,s4];
                end
            % one twin to be considered                   
            else
                b = [s1,s2,s3,s4];
                if t(2)>t(1)+t(4) && t(3)>t(1)+t(4)
                    a = [s2,s3];
                elseif t(1)>t(2)+t(3) && t(4)>t(2)+t(3)
                    if t(1) > t(4)
                        a = [s1,s1,s1,s4];
                    elseif t(1) < t(4)
                        a = [s1,s4,s4,s4];
                    else
                        a = [s1,s4];
                    end
                else
                    
                end
            end
                      
            r = a;
            if stepOn > 0
                r = '****';
            end
            rBackup = b;
        otherwise
            r = '#';
            rBackup = '#';
    end
end

function rt = extand01(r,flag,seg,limitNum)
    rt = r;
    if size(seg,2) == 1 && limitNum >= 4
        if flag == 2
            if r == '01'
                rt = '0011';
            elseif r == '10'
                rt = '1100';
            end
        end
    end
end








