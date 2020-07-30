function y = I_profile2(x, t, tau, lam, amp)

y = zeros(size(x));

v = 5.75 / t; 
dx = x - v*t;
slam = lam + v * tau;
dlam = lam - v * tau;

y(dx>0) =  exp(-dx(dx>0)/lam) * lam * tau  / slam;

y(dx<=0) = lam * tau * ( -2*v*tau * exp(dx(dx<=0)/(v*tau)) +...
 slam * exp(dx(dx<=0)/lam) ) / (slam*dlam);

y = y * amp;
