function quickRevPhiTrialPlot(LogData)

% this function plots the frame position, LmR, freq, and 1 and 3 overlaid
% for a quick check of the output. 
% Designed specifically for the revPhi protocol of 12-10-19
% also fits 02 26 20 (rotation protocol)
%
% input
% LogData -     structure generated from TDMS files (G4_TDMS_folder2struct)


timeConv = 1000000; %converts seconds to microseconds (TDMS timestamps are in micros)
firstTP = double(LogData.ADC.Time(1)) / timeConv ;
numReps = 7; 


% adding pattern
patComInd = cellfun(@(x) strcmp(x, 'Set Pattern ID'), LogData.Commands.Name);
patComStopInd = cellfun(@(x) strcmp(x, 'Stop-Display'), LogData.Commands.Name);
funcComInd = cellfun(@(x) strcmp(x, 'Set Pattern Function ID'), LogData.Commands.Name);
patVal = LogData.Commands.Data(patComInd);
funcVal = LogData.Commands.Data(funcComInd);
patNum = cellfun(@(x) hex2dec(x(1)), patVal);
funcNum = cellfun(@(x) hex2dec(x(1)), funcVal); % 2 speeds and third for optomotor
patNum(patNum == 11) =  0; % making open loop pattern zero
% parsedPatNum = reshape(patNum, [], numReps); 

patTime = double(LogData.Commands.Time(patComInd)) /timeConv - firstTP; 
patStopTime = double(LogData.Commands.Time(patComStopInd)) /timeConv - firstTP; 
timeVec = double(LogData.Frames.Time) / timeConv - firstTP; 
timeVecV = double(LogData.ADC.Time(1,:)) / timeConv - firstTP;  % since this time vector is of a different length
timeInds = arrayfun(@(x) find(timeVec > x, 1, 'first'), patTime);
timeStopInd = find(timeVec > patStopTime, 1, 'first');
if isempty(timeStopInd) 
    timeStopInd = length(timeVec);
end
timeStopIndV = find(timeVecV > patStopTime, 1, 'first');
if isempty(timeStopIndV)
    timeStopIndV = length(timeVecV);
end
timeIndsV = arrayfun(@(x) find(timeVecV > x, 1, 'first'), patTime);
patValVec = zeros(length(timeVecV), 3); 
timeIndsV = [timeIndsV, timeStopIndV];

for ii=1:length(timeIndsV)-1
    %remove for original function
    patFuncInd = funcNum(ceil(ii/2));    % since there are twice as many pattern changes (itertrial is open loop)
%     patValVec(timeInds(ii):timeInds(ii+1), patFuncInd) = patNum(ii); 
    patValVec(timeIndsV(ii):timeIndsV(ii+1), patFuncInd) = patNum(ii); 
end

% timeInds = timeInds(1:end-1); % remove for original function
tempTimeInds = reshape(timeInds, [], numReps);
tempTimeIndsV = reshape(timeIndsV(1:end-1), [], numReps);
tempInds = [tempTimeInds(1, :), timeStopInd];
tempIndsV = [tempTimeIndsV(1, :), size(LogData.ADC.Time,2)];

% since chldren are flipped
patCol = flipud(cbrewer('qual', 'Set1', 3));
patWid = fliplr([0.5, 1, 2]); 

for rr=1:numReps
    
    figure('position', [500+100*rr, 20*(rr-1), 1000, 1300])
    
    relInds = tempInds(rr):tempInds(rr+1);
    relIndsV = tempIndsV(rr):tempIndsV(rr+1);
    axh = gobjects(1,4);

    axh(1) = subplot(4,1,1);

    plot(timeVec(relInds), LogData.Frames.Position(relInds))

    axh(2) = subplot(4,1,2);

    yyaxis left
    plot(double(LogData.ADC.Time(1,relIndsV)) / timeConv - firstTP, medfilt1(LogData.ADC.Volts(1, relIndsV), 11))
    line([timeVec(relInds(1)), timeVec(relInds(end))], [0,0], 'color', 'k')
    yyaxis right
    plot(timeVecV(relIndsV), patValVec(relIndsV, :))
    
    for cc=1:size(patValVec, 2) 
        axh(2).Children(cc).LineStyle = '-'; 
        axh(2).Children(cc).Color = patCol(cc,:); 
        axh(2).Children(cc).LineWidth = patWid(cc); 
    end

    axh(3) = subplot(4,1,3);

    yyaxis left
    plot(double(LogData.ADC.Time(1,relIndsV)) / timeConv - firstTP, LogData.ADC.Volts(4, relIndsV))
    yyaxis right
    plot(timeVecV(relIndsV), patValVec(relIndsV,:))
    
    for cc=1:size(patValVec, 2)
        axh(3).Children(cc).LineStyle = '-';
        axh(3).Children(cc).Color = patCol(cc,:); 
        axh(3).Children(cc).LineWidth = patWid(cc); 
    end


    axh(4) = subplot(4,1,4);
    hold on 

    yyaxis left
    plot(double(LogData.Frames.Time(relInds)) / timeConv - firstTP, LogData.Frames.Position(relInds))
    yyaxis right
    plot(double(LogData.ADC.Time(1,relIndsV)) / timeConv - firstTP, LogData.ADC.Volts(1, relIndsV))
    hold off


    for ii=1:4
        axh(ii).XLim = [double(LogData.ADC.Time(1,relIndsV(1))), double(LogData.ADC.Time(1,relIndsV(end)))] / timeConv - firstTP;
    end

    linkaxes(axh, 'x')
    zoom xon


end

end