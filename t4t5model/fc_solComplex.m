function f = fc_solComplex(t_tot, t_on, t_off, Ti, Tr, Td)

% This function was copied from t5_complex

t0 = t_tot(t_on);
t_tot = t_tot - t0; %shift everything so t0 is 0
t1 = t_tot(t_off);  %so now t1 is effDur

%solve where t0 = 0, so shift everything over
F1 = @(t) (exp(-t./Td).*(Td.*exp(t./Td).*(Ti - Tr) + (Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr^2.*exp(t./Td - t./Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr) + (Td^2.*exp(-t./Td))./(Td.*Ti + Td.*Tr - Ti.*Tr - Td^2);
F2 = @(t) (exp(-t./Td).*(Td^2 - Td^2.*exp(t1./Td)))./(Td.*Ti + Td.*Tr - Ti.*Tr - Td^2) + (exp(-t./Td).*((Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr^2.*exp(t./Td - t./Tr))./(Td - Tr) - (Td.*Ti^2.*exp(t./Td - t./Ti + t1./Ti))./(Td - Ti) + (Td.*Tr^2.*exp(t./Td - t./Tr + t1./Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);
f = [zeros(size(t_tot(1:t_on-1)));F1(t_tot(t_on:t_off-1));F2(t_tot(t_off:end))]';

    if any(isnan(f) | isinf(f))
        %if nans, solve the long way. t0 is 0, and t1 is effDur
            I1 = @(t) 1 - exp(-t/Ti);
            H1 = @(t) (Tr*exp(-t/Tr))/(Ti - Tr) - (Tr - Ti + Ti*exp(-t/Ti))/(Ti - Tr);
            F1 = @(t) (exp(-t/Td)*(Td*exp(t/Td)*(Ti - Tr) + (Td*Ti^2*exp(t/Td - t/Ti))/(Td - Ti) - (Td*Tr^2*exp(t/Td - t/Tr))/(Td - Tr)))/(Td*Ti - Td*Tr) + (Td^2*exp(-t/Td))/(Td*Ti + Td*Tr - Ti*Tr - Td^2);

            A = I1(t1); %find init conds for decay
            B = H1(t1);
            C = F1(t1);

            t_tot = t_tot-t1; %shift frame where t1 is 0
            t1 = 0;
            I2 = @(t) A.*exp(-t./Ti); %evaluate decays from this time (t1 = 0)
            H2 = @(t) exp(-t./Tr).*(B - (A.*Ti)./(Ti - Tr)) + (A.*Ti.*exp(-t./Ti))./(Ti - Tr);
            F2 = @(t) exp(-t./Td).*(C - ((Td.*Tr.*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr) - (A.*Td.*Ti^2)./(Td - Ti))./(Td.*Ti - Td.*Tr)) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr.*exp(t./Td - t./Tr).*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);

        counter = 0;
        t_true = 0;
        while any(isnan(f) | isinf(f)) %if still nans, evaluate each ODE
            counter = counter+1;
            if counter> 100 %no unchecked while loops in cluster
                break
            end            
            tmp_idx = find(isnan(f) | isinf(f),1) -1;  %find where function NaNs

            t_true = t_true + t1; %make t1 now relative to the old t1, so it is the time of integration until NaN
            t1 = t_tot(tmp_idx) - t_true; %t1 is integration time
            A = I2(t1); %find initial conditions at this time
            B = H2(t1);
            C = F2(t1);

            I2 = @(t) A.*exp(-t./Ti); %evaluate decays from this time
            H2 = @(t) exp(-t./Tr).*(B - (A.*Ti)./(Ti - Tr)) + (A.*Ti.*exp(-t./Ti))./(Ti - Tr);
            F2 = @(t) exp(-t./Td).*(C - ((Td.*Tr.*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr) - (A.*Td.*Ti^2)./(Td - Ti))./(Td.*Ti - Td.*Tr)) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr.*exp(t./Td - t./Tr).*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);
            f(tmp_idx+1:end) = F2(t_tot(tmp_idx+1:end) - t_tot(tmp_idx)); %add to output vector

        end %repeat while we have nans
    end
end