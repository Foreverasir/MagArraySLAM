function res = dtwWave(waveCell)
% return the result of dtw
% The template is now defined as follows:
%   OneStepSingle{}
%   OneStepDouble{}
res = {};

stable = zeros(1,5);
short = ones(1,10)*10;
long = ones(1,20)*10;


OneStepSingle(1) = {[stable,long*2,stable]'};
OneStepSingle(2) = {[stable,-long*2,stable]'};
OneStepSingle(3) = {[stable,short,-short,stable]'};
OneStepSingle(4) = {[stable,-short,short,stable]'};

OneStepDouble(1) = {[stable,short,-short,short,-short,stable]'};
OneStepDouble(2) = {[stable,-short,short,-short,short,stable]'};
OneStepDouble(3) = {[stable,long,-long,stable]'};
OneStepDouble(4) = {[stable,-long,long,stable]'};

OneStepDouble(5) = {[stable,short,-short,long,stable]'};
OneStepDouble(6) = {[stable,-short,short,-long,stable]'};
OneStepDouble(7) = {[stable,long,-short,short,stable]'};
OneStepDouble(8) = {[stable,-long,short,-short,stable]'};

OneStepDouble(9) = {[stable,short*0.5,-short*1.5,-long*1.5,stable]'};
OneStepDouble(10) = {[stable,-short*0.5,short*1.5,long*1.5,stable]'};
OneStepDouble(11) = {[stable,long*1.5,short*1.5,-short*0.5,stable]'};
OneStepDouble(12) = {[stable,-long*1.5,-short*1.5,short*0.5,stable]'};

for i = 1:size(waveCell,2)
    currentCell = waveCell{i};
    c = [];
    if size(currentCell,2) > 1
        for j = 1:size(currentCell,2)
            c = [c;dtw(OneStepSingle{1},currentCell{j}),dtw(OneStepSingle{2},currentCell{j}),dtw(OneStepSingle{3},currentCell{j}),dtw(OneStepSingle{4},currentCell{j})];
        end
    else
        c = [dtw(OneStepDouble{1},currentCell{1}),dtw(OneStepDouble{2},currentCell{1}),dtw(OneStepDouble{3},currentCell{1}),dtw(OneStepDouble{4},currentCell{1}), ...
            dtw(OneStepDouble{5},currentCell{1}),dtw(OneStepDouble{6},currentCell{1}),dtw(OneStepDouble{7},currentCell{1}),dtw(OneStepDouble{8},currentCell{1}), ...
            dtw(OneStepDouble{9},currentCell{1}),dtw(OneStepDouble{10},currentCell{1}),dtw(OneStepDouble{11},currentCell{1}),dtw(OneStepDouble{12},currentCell{1})];
    end
    res(i) = {c};
end
end

