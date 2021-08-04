% Follows non-linear reaction kinetics
%% Reaction conditions
n_steps = 4;
E_forward = [0; 0; 71.9; 85]*1000;
E_back = [78; 118; 41.6; 41.6]*1000;
pre_exp_f = [6.3*10^5; 3*10^4; 1.9*10^9; 1.3*10^12];
pre_exp_b = [4.5*10^10; 1.1*10^12; 7.6*10^11; 2.9*10^4];
T = 673;%Kelvin
% Pressure in atm
p_H2 = 0.6;
p_C4H8 = 0.2; 
p_C4H10 = 0.2;
R = 8.314; %SI
%% Calculate rate constants using Arrhenius equation
forw_rates = pre_exp_f.*exp(-E_forward/(R*T));
back_rates = pre_exp_b.*exp(-E_back/(R*T));
% Multiplying by pressure terms to get 'omega' as represented in the paper
w_f = forw_rates.*[p_H2;p_C4H8;1;1;];
w_b = back_rates.*[1;1;1;p_C4H10];
w = w_f./w_b; %Basically eqbm constant with the known pressure term 
% for minimising number of terms in later expressions
wi = 1./w;
%% Resistances
Ri = zeros(n_steps,1);
Ri(1) = 1/w_f(1)*(1 + sqrt(prod(wi(2:4))) + w(2) + sqrt((w(2)*w(3))/w(4)))^2;
Ri(2) = 1/w_f(2)*(1 + sqrt(w(1)) + w(2)/prod(w(1:4)) + sqrt(1/w(1))/w(4));
Ri(3) = 1/w_f(3)*(wi(2)*sqrt(wi(1)))*(1 + sqrt(w(1)) + w(2) + sqrt(wi(1))*wi(4))^2;
Ri(4) = 1/w_f(4)*prod(wi(1:3))*(1 + sqrt(w(1)) + w(2) + sqrt(w(1))*w(2)*w(3))^2;
%% Overall rate
rOR =( 1-1/prod(w))/sum(Ri);
