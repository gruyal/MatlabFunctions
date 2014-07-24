function funcHamd = listCreateProtocolFunctions

% function funcHamd = listCreateProtocolFunctions
%
% This function lists the available createXXXProtocol functions and
% provides a function handle to the desired one

fprintf(' \n Available createXXXProtocol functions: \n')

avFuncCell = {'createGratingProtocol'; ...
              'createGratingFadeProtocol'; ...
              'createEdgeProtocol'; ...
              'createEdgeFadeProtocol'; ...
              'createCenterSurroundProtocol'; ...
              'createMovingLinearObjectProtocol'...
              };

for ii=1:length(avFuncCell)
    fprintf('%d. %s \n', ii, avFuncCell{ii})
end

while 1
    inp = input('Please provide the desired function number:');
    if inp > 0 && inp < (length(avFuncCell)+1)
        break
    end
    fprintf('\n')
end

fprintf('selected %s \n', avFuncCell{inp})

funcHamd = str2func(avFuncCell{inp});



end