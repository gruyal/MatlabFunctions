function condVec = calcCond(tr, td)
    
t=0:0.2:5000;
t0=0;
t1 = 160;

condVec = nan(size(t));



for ii=1:length(t)
    
    if t(ii) >= t0 && t(ii) < t1  
        
        tmp = (td - exp( (-t(ii) + t0)/td )*td + (-1 + exp( (-t(ii) + t0)/tr))*tr) / (td - tr);
        
    elseif t(ii) >= t1
            
        tmp = -(exp( (-t(ii) + t0)/td )*td - exp((-t(ii) + t1)/td)*td - exp((-t(ii) + t0)/tr)*tr +exp((-t(ii) + t1)/tr)*tr)/(td - tr);
        
    end
    
    condVec(ii) = tmp;

end

end