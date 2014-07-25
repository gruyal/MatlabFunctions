function stepData = injectStepsNewPC(maxCurr, numSteps, reps)

% function data = injectSteps(maxCurr, reps, injectRate)
%
% This function injects a 500ms current steps and records the responses in both
% current and voltage traces. The function assumes a 400pA/V conversion 
% ratio (defined in the amplifier gains)   
%
% INPUT
% maxCurr - maxCurr to be injected in pA
% numSteps - number fo steps between maxCurr and rest
% reps -    (optional), number of times to repeat the injection (defualt 5)
% 
%
% OUTPUT
% stepData -    a NX2 matrix with voltage and current as both channels. 
%           Sampling rate is 10K and N is determined by both maxCurr (due to the preset rate)
%           and the number of repeats. 
% Data is atomatically saved in the working directory. 



% Defining core variables
sampRate = 100;
interStepInt = sampRate/2; % 0.5 sec between steps
direc = pwd;


if nargin < 3
    reps = 5;
end

conversionRatio = 400/1; % 400pA/V
maxVol = maxCurr/conversionRatio * 5; % conversion ratio from a velfunc to voltage

stepLength = sampRate/2; % 500ms
stepVals = linspace(0, maxVol, numSteps+1);
stepVals = stepVals(2:end); %gets rid of zero

valInds = repmat(1:numSteps, 1, reps);
permInds = randperm(length(valInds));
shufInds = valInds(permInds);

oneStep = [ones(1, stepLength), zeros(1, interStepInt)];

sigMat = zeros(length(oneStep), numSteps*reps);

for ii=1:length(shufInds)
    sigMat(:,ii) = oneStep*stepVals(shufInds(ii));
end


anaSig = vertcat(zeros(interStepInt, 1), reshape(sigMat, [], 1), zeros(interStepInt, 1));
make_velocityfunction_image(anaSig);

tempVarNames = who('global');
connChk = sum(cellfun(@(x) strcmp(x, 'tcpHandle'), tempVarNames)); 

if ~logical(connChk)
    error('No tcp connection detected - run "runPanelHost"')
end


Panel_tcp_com('set_mode', [0, 5]);
Panel_tcp_com('set_funcy_freq' , sampRate);
Panel_tcp_com('set_velfunc_id', [2, 1]);

Panel_tcp_com('start_log')
Panel_tcp_com('start')

pause(length(anaSig)/sampRate)
Panel_tcp_com('stop')
Panel_tcp_com('stop_log')


% % Plotting data
% if numSteps > 6
%     subpar = [2, ceil(numSteps/2)];
% else
%     subpar = [1, numSteps];
% end
% 
% figure;
% set(gcf, 'position', [180, 700, 1600, 300])
% axh = zeros(1, numSteps);
% for ii=1:numSteps
%     tempInds = find(anaSig == stepVals(ii));
%     spInds = SplitVec(tempInds, 'consecutive');
%     spInds = cellfun(@(x) [x', x(end) + (1:250)], spInds, 'uniformoutput', 0);
%     axh(ii) = subplot(subpar(1), subpar(2), ii);
%     hold on
%     cellfun(@(x) plot(1:length(x), stepData(x, 2), 'b'), spInds)
%     hold off
%     title(num2str(stepVals(ii)))
% end
% 
% yylim = get(axh(:), 'ylim');
% yylim = vertcat(yylim{:});
% yymax = max(yylim(:,2));
% yymin = min(yylim(:,1));
% 
% set(axh(:), 'ylim', [yymin, yymax])



end
