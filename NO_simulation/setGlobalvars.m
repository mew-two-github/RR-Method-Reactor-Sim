function setGlobalvars()
    global n_steps E_forward E_back pre_exp_f pre_exp_b R;
    % Reaction conditions
    n_steps = 8;
    E_forward = [0.6; 10^11; 10^11; 10^11; 10^13; 0.89; 10^11; 10^13];
    E_back = [10^13; 10^13; 0.001; 10^11; 0.001; 10^13; 10^11; 0.005;];
    pre_exp_f = [0;12.6;27;21;5.6;0;23.2;3.6];
    pre_exp_b = [26;36.5;0;3.6;0;32;18.4;0];
    R = 8.314; %SI
end