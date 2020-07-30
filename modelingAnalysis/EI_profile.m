function y = EI_profile(x, t, v, tau, lam, amp)

y = zeros(size(x));

dx = x - v*t;
slam = lam + v * tau;
dlam = lam - v * tau;

if(v>0)
    y(dx>0) =  exp(-dx(dx>0)/lam) * lam * tau  / slam;
    y(dx<=0) = lam * tau * ( -2*v*tau * exp(dx(dx<=0)/(v*tau)) +...
        slam * exp(dx(dx<=0)/lam) ) / (slam*dlam);
else
    y(dx<=0) =  exp(dx(dx<=0)/lam) * lam * tau  / dlam;
    y(dx>0) = lam * tau * ( 2*v*tau * exp(dx(dx>0)/(v*tau)) +...
        dlam * exp(-dx(dx>0)/lam) ) / (slam*dlam);    
end

y = y * amp;


end