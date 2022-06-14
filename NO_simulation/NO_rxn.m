clear; close all;
%% Set global variables
global n_steps E_forward E_back pre_exp_f pre_exp_b R;
setGlobalvars();
%% Microkinetic model
% Species  1:NO, 2:CO, 3:NO.S, 4:N.S, 5:O.S, 6:N2, 7:CO.S, 8:CO2.S, 9:N2O.S
% 10:N2O, 11:CO2
initial = zeros(10,1);
% Strategy:
% 1. For given feed flow rates find concentrations: ci = Fi/sum(F)*Ctotal.
% Ctotal = Preactor/RTreactor
% 2. dFi/dV = rate of change of concentration of that species
% 3. For a given volume simultaneously solve the above set of odes. i.e. 
% write a function that performs step 1 and 2 and feed it to ode45 imposing
% limits on the volume.
%% RR method
% Strategy:
% Unlike the microkinetic model, only rates that are of concern to us are
% those of the terminal species (initial reactants and final products).
% Essentially follow the same procedure, except pretend that you have only 
% the terminal species whose rates are given by the RR method eqns overall
% rate. some reactants/products might have rOR1 + rOR2 as their rates since
% this is a parallel reaction.
