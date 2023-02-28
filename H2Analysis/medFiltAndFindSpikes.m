function [medData, spikeVecI] = medFiltAndFindSpikes(data, parSt, plotF)

% function medFiltAndFindSpikes(data, parSt, plotF)
%
% This function was designed for H2 recordings (default parameters are set
% to that). it median filters the voltage trace and removes the spikes from
% the recording. It then finds the peak of each spike (by differences from
% the median filtered trace)
%
% INPUT
% data -    NX1 Voltage trace already converted to mV
% parSt -   (optional) parameter structure that can include any of the
%           following:
%           .smWin -    window size in sample from moving averge (before median
%                       filter, to get rid of small peaks). default 21
%           .medWin -   window size for median filter (in samples). default
%                       351
%           .threshV -  threshold size in mV for the min difference to
%                       detect between the original trace and the medfilt
%                       one. Default 4.5
% plotF -   (optional). logical flag for plotting original data overlaid
%           with spikes and medfilt trace
%
% OUTPUT
% medData -     NX1 median filtered trace
% spikeVec -    NX1 logical vector marking position of spike peak

smWin = 21; 
medWin = 351;
threshV = 4.5;


if nargin > 1
    
    if isfield(parSt, 'smWin')
        smWin = parSt.smWin;
        assert(smWin > 0 & floor(smWin) == ceil(smWin), 'smWin should be a positive integer')    
    end
    if isfield(parSt, 'medWin')
        medWin = parSt.medWin;
        assert(medWin > 0 & floor(medWin) == ceil(medWin), 'medWin should be a positive integer')
    end
    if isfield(parSt, 'threshV')
        threshV = parSt.threshV;
        assert(threshV > 0, 'threshV should be a positive number')
    end
end

spikeVecI = zeros(size(data));
spikeVecV = nan(size(data));
smData = smooth(data, smWin); 
medData = medfilt1(smData, medWin, 'truncate');
diffI = data - medData > threshV;
[relGV, relGI] = SplitVec(diffI, 'equal', 'firstVal', 'bracket');

relBrac = relGI(relGV == 1, :);

for jj=1:size(relBrac,1)
    [maxV, maxI] = max(data(relBrac(jj,1):relBrac(jj,2)));
    spikeVecI(maxI+relBrac(jj,1)-1) = 1; 
    spikeVecV(maxI+relBrac(jj,1)-1) = maxV;
end


if nargin == 3 && plotF == 1
    clf
    hold on 
    plot(data)
    plot(spikeVecV, 'or')
    plot(medData, 'linewidth', 2)
    hold off
end


end
