function quickTestTrialPlot(LogData)

% this function plots the frame position, LmR, freq, and 1 and 3 overlaid
% for a quick check of the output. 
% Designed specifically for the test protocol (fixation and optomotor)
%
% input
% LogData -     structure generated from TDMS files (G4_TDMS_folder2struct)

timeConv = 1000000; %converts seconds to microseconds (TDMS timestamps are in micros)
firstTP = double(LogData.Frames.Time(1)) / timeConv ;

% adding pattern
patComInd = cellfun(@(x) strcmp(x, 'Set Pattern ID'), LogData.Commands.Name);
patVal = LogData.Commands.Data(patComInd);
patNum = cellfun(@(x) str2double(x(1)), patVal);
patNum(patNum == 5) =  0; % making open loop pattern zero

patTime = double(LogData.Commands.Time(patComInd)) /timeConv - firstTP; 
timeVec = double(LogData.Frames.Time) / timeConv - firstTP; 
timeInds = arrayfun(@(x) find(timeVec > x, 1, 'first'), patTime);
patValVec = zeros(size(timeVec)); 

for ii=1:length(timeInds)-1
    patValVec(timeInds(ii):timeInds(ii+1)) = patNum(ii); 
end

figure('position', [100, 20, 1000, 1300])

axh = gobjects(1,4);

axh(1) = subplot(4,1,1);

plot(timeVec, LogData.Frames.Position)

axh(2) = subplot(4,1,2);

yyaxis left
plot(double(LogData.ADC.Time(1,:)) / timeConv - firstTP, LogData.ADC.Volts(1, :))
line([timeVec(1), timeVec(end)], [0,0], 'color', 'k')
yyaxis right
plot(timeVec, patValVec)

axh(3) = subplot(4,1,3);

yyaxis left
plot(double(LogData.ADC.Time(1,:)) / timeConv - firstTP, LogData.ADC.Volts(4, :))
yyaxis right
plot(timeVec, patValVec)


axh(4) = subplot(4,1,4);
hold on 

yyaxis left
plot(double(LogData.Frames.Time) / timeConv - firstTP, LogData.Frames.Position)
yyaxis right
plot(double(LogData.ADC.Time(1,:)) / timeConv - firstTP, LogData.ADC.Volts(1, :))
hold off


for ii=1:4
    axh(ii).XLim = [0, double(LogData.Frames.Time(end)) / timeConv - firstTP];
end

linkaxes(axh, 'x')
zoom xon

end