function f = fc_solSimpleV2(t_tot, t_on, t_off, Tr, Td)

% This function was copied from t5_simpleV2

    t0 = t_tot(t_on);
    t_tot = t_tot - t0; %shift everything so t0 is 0
    t1 = t_tot(t_off);  %so now t1 is effDur

    %solve where t0 = 0, so shift everything over
    F1 = @(t) (Td - Tr + Tr.*exp(-t./Tr))./(Td - Tr) - (Td.*exp(-t./Td))./(Td - Tr);
    F2 = @(t) exp(t1./Td).*exp(-t./Td).*((Td - Tr + Tr.*exp(-t1./Tr))./(Td - Tr) - (Td.*exp(-t1./Td))./(Td - Tr) + (Tr.*exp(t1./Td - t1./Tr).*exp(-t1./Td).*(exp(t1./Tr) - 1))./(Td - Tr)) - (Tr.*exp(t./Td - t./Tr).*exp(-t./Td).*(exp(t1./Tr) - 1))./(Td - Tr);
    f = [zeros(size(t_tot(1:t_on-1)));F1(t_tot(t_on:t_off-1));F2(t_tot(t_off:end))]';

    if any(isnan(f) | isinf(f))
        %if nans, solve the long way. t0 is 0, and t1 is effDur
            H1 = @(t) 1 - exp(-t./Tr);
            %F1 defined above

            B = H1(t1);
            C = F1(t1);

            t_tot = t_tot-t1; %shift frame where t1 is 0
            t1 = 0;
            H2 = @(t) B.*exp(-t./Tr);
            F2 = @(t) exp(-t./Td).*(C + (B.*Tr)./(Td - Tr)) - (B.*Tr.*exp(-t./Tr))./(Td - Tr);

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
            B = H2(t1);
            C = F2(t1);

            H2 = @(t) B.*exp(-t./Tr);
            F2 = @(t) exp(-t./Td).*(C + (B.*Tr)./(Td - Tr)) - (B.*Tr.*exp(-t./Tr))./(Td - Tr);
            f(tmp_idx+1:end) = F2(t_tot(tmp_idx+1:end) - t_tot(tmp_idx)); %add to output vector

        end %repeat while we have nans

    end
end