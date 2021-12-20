function [rOR1,rOR2,R] = rr_rOR(T,pNO,pCO,pN2,pCO2,pN2O,X)
    X1 = X(1); X2 = X(2)
    %% Compute pressures using conversion
    pNO = pNO*(1-2*X1-2*X2);
    pCO = pCO*(1-2*X1-X2);
    pN2 = pN2*(X1);
    pCO2 = pCO2*(1-2*X1-X2);
    pN2O = pN2O*X2;
    %% Forward and Backward rates
    %% elementary resistances
    %% eqvt resistances
    RA = 2*R1 + R2 + R4 + R5 + R8;
    RB = R3 + R2 + R4 + R5 + R8;
    RC = R6 + R7;
    R = [R1;R2;R3;R4;R5;R6;R7;R8];
    Rj = (RA*RB + RB*RC + RC*RA)/RB;
    Ri = (RA*RB + RB*RC + RC*RA)/RC;
    % Rk = (RA*RB + RB*RC + RC*RA)/RA;
    %% Compute zOR and rOR
    zOR1 = ( 1-1/prod(w1));
    zOR2 = ( 1-1/prod(w2));
    rOR1= (1-zOR1)/Rj;
    rOR2 = (11 - zOR2)/Ri;
    rOR = [rOR1;rOR2];
end