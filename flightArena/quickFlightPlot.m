function quickFlightPlot(LogData)

% this function plots the frame position, LmR, freq, and 1 and 3 overlaid
% for a quick check of the output
%
% input
% Log -     structure generated from TDMS files (G4_TDMS)folder2struct)

timeConv = 1000000; %converts seconds to microseconds (TDMS timestamps are in micros)
firstTP = double(LogData.Frames.Time(1)) / timeConv ;


clf

axh = gobjects(1,4);

axh(1) = subplot(4,1,1);

plot(double(LogData.Frames.Time) / timeConv - firstTP, LogData.Frames.Position)

axh(2) = subplot(4,1,2);

plot(double(LogData.ADC.Time(1,:)) / timeConv - firstTP, LogData.ADC.Volts(1, :))

axh(3) = subplot(4,1,3);

plot(double(LogData.ADC.Time(1,:)) / timeConv - firstTP, LogData.ADC.Volts(4, :))

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



end