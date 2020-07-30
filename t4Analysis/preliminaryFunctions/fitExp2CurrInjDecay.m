function fitParams = fitExp2CurrInjDecay(protocolStructAO, verbose)

% this function fits an exponent to the decay phase of a current injection
% protocol. Only fit the strongest injection (last) after responses are
% aligned and the mean is taken.
%
% INPUT
% protocolStructAO  -       protocol structure from a step current injection protocl
% verbose  -                (optional) logical. if true plots the results (default 1). 
%
% OUTPUT
% fitParams -               coefficients and GoF valuse from the fit in a
%                           table format


if nargin < 2
    verbose =1;
end 


fOpt = fitoptions('method', 'nonlinearleastsquare', ...
                  'lower', [-100, -100, 0], ...
                  'upper', [100, 100, 100], ...
                  'startpoint', [-50, 10, 20]);

decayMod = fittype('a+b*(exp(-x/c))', 'coefficients', {'a', 'b', 'c'}, 'options', fOpt);


datToMvConv = 10;
datToMsConv = 10^-3;


allStimInds = vertcat(protocolStructAO.stim.relInds);

totStim = max(allStimInds(:,4));
relIndsSt = getStimInds(protocolStructAO, [1, nan, nan, totStim]);
relInds = relIndsSt.inds;
respCell = cell(1, length(relInds));

for ii=1:length(relInds)

    stimDat = protocolStructAO.stim(relInds(ii)).data{3};
    decayTimeInd = find(stimDat(:,2) >0, 1, 'last');
    decayTimestamp = stimDat(decayTimeInd, 1);
    
    respDat = protocolStructAO.stim(relInds(ii)).data{1}(:, [1,3]);
    relRespDat = respDat(respDat(:,1) - decayTimestamp > 0, :);
    respCell{ii} = relRespDat(:,2) * datToMvConv;
    
end

timeVec = relRespDat(:,1) - relRespDat(1,1);

numSamp = min(cellfun(@length, respCell));

cropResp = cellfun(@(x) x(1:numSamp), respCell, 'uniformoutput', 0);

meanResp = mean([cropResp{:}],2);
timeVec = timeVec(1:numSamp) * datToMsConv;


[expFit, expGoF] = fit(timeVec, meanResp, decayMod);


if verbose
    figure
    plot(timeVec, meanResp)
    hold on
    plot(expFit)
    hold off
end

offset = expFit.a;
scale = expFit.b;
tau = expFit.c;
rSq = expGoF.rsquare;

fitParams = table(offset, scale, tau, rSq);


end




