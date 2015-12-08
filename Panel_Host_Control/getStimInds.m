function stimInds = getStimInds(pStruct, desiredStim)

% function stimInds = getStimInds(pStruct)
%
% This function takes the protocolStructure generated by any of the
% createXXXProtocol functions and finds the indices for a particular
% stimulus based on the relInds property of the stimulus.
% 2015 12
% Function has been modified to include combined protocols, and therefore 
% to look for combInds instead of relInds. Zeros can no loger be used as wild cards 
% (because if the way they are used in combInds). NaNs are wild cards.
% cards
%
% INPUT
%
% pStruct -         protocolStruct with stim and relInds fields for each stimulus
%                   (as generated by the createXXXProtocol.
% desiredStim -     1X4 vector for the desired relInds values. NaNs can be
%                   used as wild cards. 
%
% OUTPUT
%
% stimInds -        structure with val and Inds fields for each unique combination
%                   of the results

assert(isvector(desiredStim), 'desiredStim should be a 1X4 vector')
assert(length(desiredStim) == 4, 'desiredStim should be a 1X4 vector')

assert(isfield(pStruct, 'stim'), 'protocolStructure is missing stim field')

if isfield(pStruct.stim, 'relInds')
    allInds = vertcat(pStruct.stim.relInds);
elseif isfield(pStruct.stim, 'combInds')
    allInds = vertcat(pStruct.stim.combInds);
else
    error('stim field is missing relInds or combInds field')
end


relvInds = allInds(:, ~isnan(desiredStim)); % chooses only columns where wild cards aren't present
relvStim = desiredStim(~isnan(desiredStim)); 

allRelvStim = allInds(all(bsxfun(@eq, relvInds, relvStim), 2), :);
uRelvStim = unique(allRelvStim, 'rows');
assert(~isempty(uRelvStim), 'No matching inds to desired inds')

for ii=1:size(uRelvStim)
    stimInds(ii).val = uRelvStim(ii, :);
    stimInds(ii).inds = find(all(bsxfun(@eq, allInds, uRelvStim(ii, :)), 2))'; % transposed so that it would be printed more easily on screen
end


end