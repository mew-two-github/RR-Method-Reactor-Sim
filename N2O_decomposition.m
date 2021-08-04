% Follows linear reaction kinetics
%% Reaction Conditions
n_steps = 5;
E_forward = [0; 41; 0; 135; 117]*1000;
E_back = [16; 138; 28; 239; 0]*1000;
pre_exp_f = [10^6; 10^13; 10^13; 4.5*10^11; 10^13];
pre_exp_b = [10^13; 10^13; 10^13; 10^13; 10^6];
p_N2O = 100*10^-6; p_N2 = 100*10^-6; p_O2 = 50*10^-6; %in ppm
R = 8.314;
T = 800; %K
%% Calculate rate constants using Arrhenius equation
forw_rates = pre_exp_f.*exp(-E_forward/(R*T));
back_rates = pre_exp_b.*exp(-E_back/(R*T));
% Multiplying by pressure terms to get 'omega' as represented in the paper
w_f = forw_rates.*[p_N2O;1;p_N2O;1;1];
w_b = back_rates.*[1;p_N2;1;p_N2;p_O2];
w = w_f./w_b; %Basically eqbm constant with the known pressure term 
% for minimising number of terms in later expressions
%% Resistances
Ri = zeros(n_steps,1);
Ri(1) = 1./w_f(1)*(1+1/prod(w(2:5))+1/prod(w(3:5))+1/prod(w(4:5))+1/w(5));
Ri(2) = 1/w(1)/w_f(2)*(1+w(1)+1/prod(w(3:5))+1/prod(w(4:5))+1/w(5));
Ri(3) = 1/w_f(3)*1/(prod(w(1:2)))*(1+w(1)+prod(w(1:2))+1/prod(w(4:5))+1./w(5));
Ri(4) = 1/w_f(4)*1/prod(w(1:3))*(1+w(1)+prod(w(1:2))+prod(w(1:3))+1/w(5));
Ri(5) = 1/w_f(5)/prod(w(1:4))*(1+w(1)+prod(w(1:2))+prod(w(1:3))+prod(w(1:4)));
%% Overall rate
rOR =( 1-1/prod(w))/sum(Ri);
