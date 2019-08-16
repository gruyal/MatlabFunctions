function [medTime, varargout] = getFrameTimes(pStruct, framePos)

% function [medTime, varargout] = getFrameTimes(pStruct, frameNums)
%
% This function gets the recorded time from each of the frames given in
% framePos. Reports the median of frame time, and as an optional output also reports 0.1 and 0.9 quantiles
%
% INPUT
%
% pStruct -         porotocolStructure as generated by runPosFuncProtocol (with
%                   .stim field for each stimulus and .data field for each stim)
% framePos -        X positions of the pattern presented at a particular stim
%                   (extracted by consolidateData). Can be integer of 1XN integer vector
%
% OUTPUT
% medTime -         median time for each of the positions given in framePos
% qTimes -          (optional). if more than one output is given also reports the
%                   Q0.1 and Q0.9 of each time (given as a table. 
%
% Note: framePos should be unique, if the given position appears more than once function will abort 

assert(isvector(framePos), 'framePos should be either an integer or a vector')
assert(isfield(pStruct, 'stim'), 'pStruct must have a stim field')
assert(isfield(pStruct.stim, 'data'), 'pStruct.stim must have a data field')

numStim = length(pStruct.stim);
times = nan(numStim, length(framePos));

for ii=1:numStim
    
    if ~isempty(pStruct.stim(ii).data)
        tempTimes = pStruct.stim(ii).data{2};
        for jj=1:length(framePos)
            logInd = tempTimes(:,2) == framePos(jj);
            switch sum(logInd)
                case 2
                    fprintf('Aborted in Stim %d\n', ii)
                    error('Position appears more than once in Stim')
                case 1
                    times(ii, jj) = tempTimes(logInd, 1);
                case 0
                    fprintf('Found no framePos %d in stim %d\n', framePos(jj), ii) 
            end
        end
    end
end

medTime = nanmedian(times);

if nargout > 1
    qTimes = quantile(times, [0.1, 0.9]);
    varargout{1} = qTimes;
end







end