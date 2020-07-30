function gratTab = normMovBarShiftTable(movBarShiftTable, relPD, relMaxE)

% function gratTab = normMovBarShiftTable(movBarShiftSt, PDFlag)
%
% This function takes the table from movebarshift protocol and adds the
% proper normilized columns (also correct so that trajectories indicies
% would be useful across protocols that both used the components and didnt)
%
%
% INPUT 
% movBarShiftTable -        Table from movBarShift protocol
% relPD -                   flag for whether PD is with the orientation of
%                           the protocol or reversed (labeled 1 or -1). 
%                           Needed since orientation for this protocol runs
%                           only from 0-3 (like singleBar) 
% relMaxE -                 max excitatory peak for the positions (new
%                           zero)
%
% OUTPUT
% gratTab -                 saem as the input table with the added relevant
%                           columns: 
%                           normStartPos, normEndPos,normShiftPos, normPosSeq
%                           all coreccted for maxEPos and PD  
%                           trajLabel - labels so that trajectories can be
%                           compared with protocols that had components in
%                           them 
%                           normTrajInd - trajInd corrected w/ the
%                           appropriate labels
%
% !!! NOTE: function is designed only for full trajectories and not for
% their components - will exclude them from the table


assert(ismember(relPD, [-1,1]), 'relPD should be either -1 or 1')
gratTab = movBarShiftTable(ismember(movBarShiftTable.trajInd, 1:4), :); 

if height(gratTab) < height(movBarShiftTable)
    warning('!!! excluded component trajectories  !!!' ) 
end

trajLabs = {'PN'; 'PP'; 'NP'; 'NN'};
altTajInds = [3;4;1;2];

% dont need to take width into account since always label the leading edge
% (unlike w/ singleBar)

if relPD == 1
    normStartPos = (gratTab.startPos - relMaxE) * relPD;
    normEndPos = (gratTab.endPos - relMaxE) * relPD;
    normPosSeq = cellfun(@(x) (x-relMaxE) * relPD, gratTab.posSeq, 'uniformoutput', 0);
    
    trajLabel = trajLabs(gratTab.trajInd);
    normTrajInd = gratTab.trajInd; 
else
    normStartPos = (gratTab.startPos - relMaxE - gratTab.width + 1) * relPD;
    normEndPos = (gratTab.endPos - relMaxE - gratTab.width + 1) * relPD;
    normPosSeq = cellfun(@(x,y) (x - relMaxE - y + 1 ) * relPD, gratTab.posSeq, arrayfun(@(x) {x}, gratTab.width), 'uniformoutput', 0);
    
    normTrajInd = altTajInds(gratTab.trajInd);
    trajLabel  = trajLabs(normTrajInd);
end
    
% if any(isnan(gratTab.shiftPos))
%     
%    normShiftPos = (gratTab.shiftPos - relMaxE) * relPD; 
%         
% else
%     normShiftPos = gratTab.shiftPos; 
%     relI = gratTab.shiftPos < 999; % in case the table already replaced NaNs with 999
%     normShiftPos(relI) = (gratTab.shiftPos(relI) - relMaxE) * relPD; 
% end

% replace NaN with 999 to allow for row comparisons (for traj comparisons) 
% normShiftPos(isnan(normShiftPos)) = 999;
gratTab.shiftPos(isnan(gratTab.shiftPos)) = 999;


gratTab = [gratTab, table(normStartPos, normEndPos, normPosSeq, trajLabel, normTrajInd)];

end