close all;clear;
%% Initial Conditions
pN2O = 100*10^-6; pN2 = 100*10^-6; pO2 = 50*10^-6; 
T0 = 800; %K
V = 25; % Volume of the reactor in m^3
[rOR,R] = func(pN2O,pN2O,pN2,pO2,T0); %Just to check whether the function is working fine
% %% Solving for time
% pFinal = 0.001;
% fun = @(p)(1./func(p,pN2O,pN2,pO2,T0));
% t = -1*integral(fun,pN2O,pFinal)*V/(8.314*T0);
% %% Implement it for different pFinals 
% % so we get rates and resistances at different instants
% n = 30;
% p = linspace(pN2O,pFinal,n);
% times = zeros(n,1); rates = zeros(n,1);Rs = zeros(5,n);
% for i = 1:30
%     times(i) = integral(fun,p(i),pN2O)/(8.314*T0);
%     if times(i) < 0  % Happens when pFinal is less than eqbm pressure
%         times = times(1:i-1);
%         rates = rates(1:i-1);
%         p = p(1:i-1);
%         Rs = Rs(:,1:i-1);
%         break
%     end
%     [rates(i),Rs(:,i)] = func(p(i),pN2O,pN2,pO2,T0);
% end
% plot(times,rates);
% xlabel('time');ylabel('rates');title('r_O_R as a f(t)');
% figure();
% for i = 1 : 5
%     subplot(3,2,i);
%     plot(times,Rs(i,:));
%     xlabel('time');ylabel('R'+string(i)); title('Plot of resistances as a f(time)');
% end
%% Solving for Pexit
t0 = 5*10^-7;
fun = @(t,p)(-1*func(p,pN2O,pN2,pO2,T0)*8.314*T0);
[t,P] = ode45(fun,[0,t0],pN2O);
plot(t,P);
title('P_N_2_O_e_x_i_t vs t at 800K');ylabel('Pressure'); xlabel('Time');
T = 600:25:1000;
Pexit = length(T);
for i = 1:length(T)
    fun = @(t,p)(-1*func(p,pN2O,pN2,pO2,T(i))*8.314*T(i));
    [t,P] = ode45(fun,[0,t0],pN2O);
    Pexit(i) = P(end);
end
figure();
plot(T,Pexit);
title('P_{N2O,exit} vs T'); ylabel('Pressure'); xlabel('Temperature');

%% function to evaluate the Resistances and rOR
function [rOR, Ri] = func(p_N2O,pN2O0,pN20,pO20,T)
    % Reaction Conditions
    n_steps = 5;
    E_forward = [0; 41; 0; 135; 117]*1000;
    E_back = [16; 138; 28; 239; 0]*1000;
    pre_exp_f = [10^6; 10^13; 10^13; 4.5*10^11; 10^13];
    pre_exp_b = [10^13; 10^13; 10^13; 10^13; 10^6];
    %p_N2O = 100*10^-6; p_N2 = 100*10^-6; p_O2 = 50*10^-6; %in ppm
    R = 8.314;
    % Assuming their initial values
    p_O2 = (pN2O0-p_N2O)/2 + pO20;
    p_N2 = (pN2O0-p_N2O) + pN20;
    %T = 800; %K
    % Calculate rate constants using Arrhenius equation
    forw_rates = pre_exp_f.*exp(-E_forward/(R*T));
    back_rates = pre_exp_b.*exp(-E_back/(R*T));
    % Multiplying by pressure terms to get 'omega' as represented in the paper
    w_f = zeros(5,length(p_N2O));
    w_b = w_f;
    w = w_f;
    for i = 1 : length(p_N2O)
        w_f(:,i) = forw_rates.*[p_N2O(i);1;p_N2O(i);1;1];
        w_b(:,i) = back_rates.*[1;p_N2(i);1;p_N2(i);p_O2(i)];
        w = w_f(:,i)./w_b(:,i); 
    end
    % Resistances
    Ri = zeros(n_steps,length(p_N2O));
    Ri(1,:) = 1./w_f(1,:)*(1+1/prod(w(2:5,:))+1./prod(w(3:5,:))+1./prod(w(4:5,:))+1./w(5,:));
    Ri(2,:) = 1./w(1,:)./w_f(2,:)*(1+w(1,:)+1./prod(w(3:5,:))+1./prod(w(4:5,:))+1./w(5,:));
    Ri(3,:) = 1./w_f(3,:)*1./(prod(w(1:2,:))).*(1+w(1,:)+prod(w(1:2,:))+1./prod(w(4:5,:))+1./w(5,:));
    Ri(4,:) = 1./w_f(4,:)*1./prod(w(1:3,:)).*(1+w(1,:)+prod(w(1:2,:))+prod(w(1:3,:))+1./w(5,:));
    Ri(5,:) = 1./w_f(5,:)./prod(w(1:4,:)).*(1+w(1,:)+prod(w(1:2,:))+prod(w(1:3,:))+prod(w(1:4,:)));
    % Overall rate
    rOR =( 1-1./prod(w))./sum(Ri);
end
