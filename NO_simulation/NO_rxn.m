clear; close all;
%% Set global variables
global n_steps E_forward E_back pre_exp_f pre_exp_b R;
setGlobalvars();
%% Microkinetic model
% Species  1:NO, 2:CO, 3:NO.S, 4:N.S, 5:O.S, 6:N2, 7:CO.S, 8:CO2.S, 9:N2O.S
% 10:N2O, 11:CO2
initial = zeros(10,1);

