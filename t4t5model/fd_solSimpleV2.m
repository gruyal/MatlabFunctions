function f = fd_solSimpleV2(t_tot, t_on, t_off, Ti, Tr, Td, b)

% This function was copied from t5_simpleV2

    t0 = t_tot(t_on);
    t_tot = t_tot - t0;
    t1 = t_tot(t_off);

    A = b;

    t_tot = t_tot - t1; %shift frame again so that t1 is 0
    t1 = 0;

    F2 = @(t) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (A.*Td.*Ti.*Tr.*exp(t./Td - t./Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr) - (A.*Td.*Ti.*exp(-t./Td))./(Td.*Ti + Td.*Tr - Ti.*Tr - Td^2);
    f = [zeros(size(t_tot(1:t_off-1)));F2(t_tot(t_off:end))]';

    if any(isnan(f) | isinf(f))
        %if nans define all functions explicitly and solve piecewise
        A = b;
        B = 0;
        C = 0;
        I2 = @(t) A.*exp(-t./Ti);
        H2 = @(t) exp(-t./Tr).*(B - (A.*Ti)./(Ti - Tr)) + (A.*Ti.*exp(-t./Ti))./(Ti - Tr);
        F2 = @(t) exp(-t./Td).*(C - ((Td.*Tr.*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr) - (A.*Td.*Ti^2)./(Td - Ti))./(Td.*Ti - Td.*Tr)) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr.*exp(t./Td - t./Tr).*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);
        f = [zeros(size(t_tot(1:t_off-1)));F2(t_tot(t_off:end))]';

        counter = 0;
        t_true = t1;
        while any(isnan(f) | isinf(f))
            counter = counter+1;
            if counter> 100 %no unchecked while loops in cluster
                break
            end            
            tmp_idx = find(isnan(f) | isinf(f),1) -1;  %find where function NaNs
            t_true = t_true+t1; %keep running index of last stopping point
            t1 = t_tot(tmp_idx) - t_true; %t1 is integration time, so difference from current stopping point and last stopping point

            A = I2(t1); %evaluate funciton just before nans. these are new init conds to decay
            B = H2(t1);
            C = F2(t1);

            I2 = @(t) A.*exp(-t./Ti); %these functions evaluate from 0
            H2 = @(t) exp(-t./Tr).*(B - (A.*Ti)./(Ti - Tr)) + (A.*Ti.*exp(-t./Ti))./(Ti - Tr);
            F2 = @(t) exp(-t./Td).*(C - ((Td.*Tr.*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr) - (A.*Td.*Ti^2)./(Td - Ti))./(Td.*Ti - Td.*Tr)) - (exp(-t./Td).*((A.*Td.*Ti^2.*exp(t./Td - t./Ti))./(Td - Ti) - (Td.*Tr.*exp(t./Td - t./Tr).*(A.*Ti - B.*Ti + B.*Tr))./(Td - Tr)))./(Td.*Ti - Td.*Tr);

            f(tmp_idx+1:end) = F2(t_tot(tmp_idx+1:end) - t_tot(tmp_idx)); %add to output vector
        end %repeat while we have nans

    end
end
