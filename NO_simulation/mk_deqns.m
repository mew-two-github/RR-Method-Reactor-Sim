function d = mk_deqns(c,T)
    % Species  1:NO, 2:CO, 3:NO.S, 4:N.S, 5:O.S, 6:N2, 7:CO.S, 8:CO2.S, 9:N2O.S
    % 10:N2O, 11:CO2
    k_f = pre_exp_f*exp(-1*E_forward/(R*T))
    k_b = pre_exp_b*exp(-1*E_back/(R*T))
    r = rates(c,k_f,k_b);
    
end
function rk = rates(c,k_f,k_b)
    rk = zeros(8,1);
    % Species  1:NO, 2:CO, 3:NO.S, 4:N.S, 5:O.S, 6:N2, 7:CO.S, 8:CO2.S, 9:N2O.S
    % 10:N2O, 11:CO2
    % equations incomplete but simple differential eqns
    rk(1) = kf(1)*c(1) - k_b(1)*c(3);
    rk(2) = kf(2) - k_b(2);
    rk(3) = kf(3) - k_b(3);
    rk(4) = kf(4) - k_b(4);
    rk(5) = kf(5) - k_b(5);
    rk(6) = kf(6) - k_b(6);
    rk(7) = kf(7) - k_b(7);
    rk(8) = kf(8) - k_b(8);
end