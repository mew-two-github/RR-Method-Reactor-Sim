close all;clear;
p_N2O = 100*10^-6; p_N2 = 100*10^-6; p_O2 = 50*10^-6; %in ppm
T0 = 800; %K
[rOR,R] = func(p_N2O,p_O2,p_N2,T0);
%% Temperature effects
T = 600:25:1000;
rOR_T = zeros(length(T)); R_T = zeros(length(T),5);
for i = 1 : length(T)
    [rOR_T(i),R_T(i,:)] = func(p_N2O,p_O2,p_N2,T(i));
end
figure();
plot(T,rOR_T,T,rOR_T,'x'); title('T vs rOR');
figure();
plot(T,R_T(:,1),T,R_T(:,2),T,R_T(:,3),T,R_T(:,4),T,R_T(:,5)); title('R vs rOR');
legend('R1','R2','R3','R4','R5');
figure();
for i = 1:5
    subplot(2,3,i);
    plot(T,R_T(:,i));title(i);xlabel('T(K)');ylabel('Resistance');
end
%% Pressure effects
%% PN2
pN2 = linspace(10*10^-6,500*10^-6,15);
rOR_p = zeros(length(pN2)); R_p = zeros(length(pN2),5);
for i = 1 : length(pN2)
    [rOR_p(i),R_p(i,:)] = func(p_N2O,p_O2,pN2(i),T0);
end
figure();
plot(pN2,rOR_p,pN2,rOR_p,'x'); title('pN_2 vs rOR');
figure();
plot(pN2,R_p(:,1),pN2,R_p(:,2),pN2,R_p(:,3),pN2,R_p(:,4),pN2,R_p(:,5)); title('R vs rOR');
legend('R1','R2','R3','R4','R5');
figure();
for i = 1:5
    subplot(2,3,i);
    plot(pN2,R_p(:,i));title(i);xlabel('pN_2(ppm)');ylabel('Resistance');
end
%% PN2O
pN2 = linspace(10*10^-6,500*10^-6,15);
rOR_p = zeros(length(pN2)); R_p = zeros(length(pN2),5);
for i = 1 : length(pN2)
    [rOR_p(i),R_p(i,:)] = func(pN2(i),p_O2,p_N2,T0);
end
figure();
plot(pN2,rOR_p,pN2,rOR_p,'x'); title('PN_2O vs rOR');
figure();
plot(pN2,R_p(:,1),pN2,R_p(:,2),pN2,R_p(:,3),pN2,R_p(:,4),pN2,R_p(:,5)); title('R vs rOR');
legend('R1','R2','R3','R4','R5');
figure();
for i = 1:5
    subplot(2,3,i);
    plot(pN2,R_p(:,i));title(i);xlabel('pN_2O(ppm)');ylabel('Resistance');
end
% %% PO2
% pO2 = linspace(10*10^-6,500*10^-6,15);
% rOR_p = zeros(length(pO2)); R_p = zeros(length(pO2),5);
% for i = 1 : length(pO2)
%     [rOR_p(i),R_p(i,:)] = func(p_N2O,pO2(i),p_N2,T0);
% end
% figure();
% plot(pO2,rOR_p,pO2,rOR_p,'x'); title('PO2 vs rOR');
% figure();
% plot(pO2,R_p(:,1),pO2,R_p(:,2),pO2,R_p(:,3),pO2,R_p(:,4),pO2,R_p(:,5)); title('R vs rOR');
% legend('R1','R2','R3','R4','R5');
% figure();
% for i = 1:5
%     subplot(2,3,i);
%     plot(pO2,R_p(:,i));title(i);xlabel('pO2(ppm)');ylabel('Resistance');
% end
%% function to evaluate the Resistances and rOR
function [rOR, Ri] = func(p_N2O,p_O2,p_N2,T)
    % Reaction Conditions
    n_steps = 5;
    E_forward = [0; 41; 0; 135; 117]*1000;
    E_back = [16; 138; 28; 239; 0]*1000;
    pre_exp_f = [10^6; 10^13; 10^13; 4.5*10^11; 10^13];
    pre_exp_b = [10^13; 10^13; 10^13; 10^13; 10^6];
    %p_N2O = 100*10^-6; p_N2 = 100*10^-6; p_O2 = 50*10^-6; %in ppm
    R = 8.314;
    %T = 800; %K
    % Calculate rate constants using Arrhenius equation
    forw_rates = pre_exp_f.*exp(-E_forward/(R*T));
    back_rates = pre_exp_b.*exp(-E_back/(R*T));
    % Multiplying by pressure terms to get 'omega' as represented in the paper
    w_f = forw_rates.*[p_N2O;1;p_N2O;1;1];
    w_b = back_rates.*[1;p_N2;1;p_N2;p_O2];
    w = w_f./w_b; 
    % Resistances
    Ri = zeros(n_steps,1);
    Ri(1) = 1./w_f(1)*(1+1/prod(w(2:5))+1/prod(w(3:5))+1/prod(w(4:5))+1/w(5));
    Ri(2) = 1/w(1)/w_f(2)*(1+w(1)+1/prod(w(3:5))+1/prod(w(4:5))+1/w(5));
    Ri(3) = 1/w_f(3)*1/(prod(w(1:2)))*(1+w(1)+prod(w(1:2))+1/prod(w(4:5))+1./w(5));
    Ri(4) = 1/w_f(4)*1/prod(w(1:3))*(1+w(1)+prod(w(1:2))+prod(w(1:3))+1/w(5));
    Ri(5) = 1/w_f(5)/prod(w(1:4))*(1+w(1)+prod(w(1:2))+prod(w(1:3))+prod(w(1:4)));
    %% Overall rate
    rOR =( 1-1/prod(w))/sum(Ri);
end
