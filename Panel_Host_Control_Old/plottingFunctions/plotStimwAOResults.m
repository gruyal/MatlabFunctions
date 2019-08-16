function varargout = plotStimwAOResults(pStructAO, uStimNum, fh)

% function plotStimwAOResults(pStructAO, uStimNum)
%
% This function take the protocol structure of an AO experiment and plots
% the results of one stimulus with all the AO combinations. 
%
% INPUT
%
% pStructAO -       Structure generated by runProtocolStructwAO. Should
%                   include the following fields.
%   .stim   -       stimulus structure that contains .data field and
%                   .comInds (read in runPosFunProtocolwAO)
% uStimNum  -       Stimulus number as indexed by comInds. 
% fh        -       figure handle (optiona). If given all axes will be
%                   plotted on this figure. If not a clf will be
%                   used.
%                   This is to allow an easier use of
%                   plotStimStructAOResults
%
% OUTPUT
% axh       -       axes handle to the plotted axes
%
% NOTE function assume AO channel was used as TTL pulse and will look draw
% gray rectangles when AO > 4.5
%
currCh = 3;
msConvFac = 10^3;

assert(isfield(pStructAO, 'stim'), 'protocol structure is missing stim field')
assert(isfield(pStructAO.stim, 'data'), 'Stim is missing data field')
assert(isfield(pStructAO.stim, 'combInds'), 'Stim is missing combInds field')

allStimCombInds = vertcat(pStructAO.stim.combInds);
uniStimInds = unique(allStimCombInds, 'rows');

assert(ismember(uStimNum, uniStimInds(:,1)), 'numStim is not included in combInds(:,1)')

relUStim = uniStimInds(uniStimInds(:,1) == uStimNum, :);
numUniStim = size(relUStim, 1);

allRelStimInds = cell(1, numUniStim);
for ii=1:numUniStim
    allRelStimInds{ii} = find(arrayfun(@(x) isequal(allStimCombInds(x, :), relUStim(ii,:)), 1:length(allStimCombInds))); % could also used a simple for loop
end

posCell = generatePositionCell(0.05, 0.95, 0.05, 0.95, -1, 0.05, numUniStim);
axh = zeros(1, numUniStim);

% Colors
cMap = cbrewer('seq', 'Reds', max(cellfun(@length, allRelStimInds)));

rectCol = [1, 1, 1] * 0.9;

if nargin < 3
    clf
else
    figure(fh)
end

minYY = -60;
maxYY = -60;

for ii=1:numUniStim
    axh(ii) = axes('position', posCell{ii});
    hold on
    for jj=1:length(allRelStimInds{ii})
        
        tempData = pStructAO.stim(allRelStimInds{ii}(jj)).data{1};
        tempPos = pStructAO.stim(allRelStimInds{ii}(jj)).data{2};
        tempAO = pStructAO.stim(allRelStimInds{ii}(jj)).data{3};
        [spAOBrac, spAOVal] = SplitVec(tempAO(:,2), 'equal', 'bracket', 'firstval');
        counter = 0;
        firstPulse = 0;
        aoTimeStp = (tempAO(:,1) - tempData(1,1))/msConvFac;
        for kk=1:length(spAOVal)
            
            if spAOVal(kk) > 4.5
                counter = counter+1;
                if counter == 1
                    firstPulse = aoTimeStp(spAOBrac(kk, 1))-1;
                end
                rectangle('position', [aoTimeStp(spAOBrac(kk,1))-firstPulse, -100, aoTimeStp(spAOBrac(kk,2))-aoTimeStp(spAOBrac(kk,1)), 200], ...
                          'facecolor', rectCol, 'edgecolor', rectCol)
            end
        end
        
        tempLength = size(tempPos, 1);
        currentDat = tempData(:,currCh)*10;
        yVal = quantile(currentDat, 0.1);
        
        %colInd = 2*(jj-1) + 1;
        
        plot((tempData(:,1)-tempData(1,1))/msConvFac-firstPulse, currentDat, 'linewidth', 2, 'color', cMap(jj, :))
        plot(double((tempPos(:,1) - tempData(1,1)))/msConvFac-firstPulse, ones(1, tempLength)* yVal, ...
             'o', 'markerfacecolor', cMap(jj, :), 'markeredgecolor', cMap(jj, :))
        % subtract tempData(1,1) since the clock doesn't start at zero
        % divide by 10^3 to convert to ms (clock is 1MHz) and multioly by
        % 10 to convert to mV
        
        ind5 = find(tempPos(:,2) == 5);
        markInd = ind5:5:tempLength;
        plot(double((tempPos(markInd,1) - tempData(1,1)))/msConvFac-firstPulse, ones(1, length(markInd))* yVal, 'o', 'markerfacecolor', 'k', 'markeredgecolor', 'k')
        
        
        tempYMin = min(currentDat);
        tempYMax = max(currentDat);
        
        if tempYMin < minYY
            minYY = tempYMin;
        end
        if tempYMax > maxYY
            maxYY = tempYMax;
        end
        
    end
    hold off
end

if nargout == 1
    varargout{1} = axh;
elseif nargout == 2
    varargout{1} = axh;
    varargout{2} = [minYY, maxYY];
end



end
