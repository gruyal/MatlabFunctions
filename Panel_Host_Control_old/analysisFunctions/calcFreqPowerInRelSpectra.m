function freqRes = calcFreqPowerInRelSpectra(dataVec, rangeMat, sampF)

% function freqRes = calcFreqPowerInRelSpectra(dataVec, rangeMat)
%
% This function does fft on the data and provides the results in the given
% ranges. 
%
%
% INPUT 
% 
% dataVec -             NX1 vector, with first vector being time and second
%                       measured results. 
%
% rangeMat -            MX2 matrix. Each of the M pairs is a range of freq (in Hz) from which results will be provided
%
% sampF -               (optional) sampling frequency. If not given assumed
%                       20KHz
%
% OUTPUT
%
% freqRes -             max power for each of the ranges. 



assert(all(diff(rangeMat,[],2) > 0), 'ranges provided are inadequate, make sure second index is larger than first')

if nargin < 3
    sampF = 20000; 
end

stimLen = length(dataVec); 

fftRes = fft(dataVec);

P2 = abs(fftRes/stimLen);
P1 = P2(1:floor(stimLen/2)+1);
P1(2:end-1) = 2*P1(2:end-1);

fRes = sampF*(0:stimLen/2)/stimLen;

freqRes = zeros(1, size(rangeMat,1));

for ii=1:length(freqRes)
    tempI = fRes > rangeMat(ii,1) & fRes < rangeMat(ii,2);
    freqRes(ii) = max(P1(tempI));
end








end