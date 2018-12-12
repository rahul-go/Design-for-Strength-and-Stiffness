%% Design Fatigue Calculator
% Author: Harry Whinnery
% Class: ME 328-08
%
% The critical section for fatigue failure on a shaft has been determined to be the location of
% the notch radius shown below. To aid in the design of the shaft you are asked to generate
% curves showing the effect on the factor of safety of varying the diameter of the stepped
% down portion, d, and the radius of the notch, r.

%% Clear Workspace
clear figures
clear all;
clc;



%% Inputs
D = 3;                      % diameters
L = 5;                      % step length(inches)
a = 10;                     % arm length(inches)
Pa = 750;                   % alternating load
Pm = 1500;                  % mean load

%% Constants 
Sut = 210;                   % ultimate tensile strength (ksi)
Sy = 180;                    % yield strength (ksi)
if Sut <= 200               
    SePrime = .5*Sut;        % (ksi)
else
    SePrime = 100;
end

dratio = 1.25:0.25:2.5;
rratio = 0.05:0.05:0.20;

a1 = 2.7;                    % ground surface finish
b1 = -0.265;                 % Table 6-2 (P.296)

%% Analysis

for i = 1:length(dratio)
    d = D/dratio(i);
    h = (D/2)-(d/2);
    I = pi()*d^4/64;
    J = I*2;
    for j = 1:length(rratio)
        r = rratio(j)*d;
        %% Finding Kt and Kts
        % Kt for Bending
        if(h/r)<=2
            c1 = .947 + 1.206*sqrt(h/r) - .131*h/r;
            c2 = .022 - 3.405*sqrt(h/r) + .915*h/r;
            c3 = .869 + 1.777*sqrt(h/r) - .555*h/r;
            c4 = -.810 + .422*sqrt(h/r) - .260*h/r;
        else
            c1 = 1.232 + .832*sqrt(h/r) - .008*h/r;
            c2 = -3.813 + .968*sqrt(h/r) - .260*h/r;
            c3 = 7.423 - 4.868*sqrt(h/r) + .869*h/r;
            c4 = -3.839 + 3.070*sqrt(h/r) - .600*h/r;
        end
        Kt = c1 + c2*(2*h/D) + c3*(2*h/D)^2 + c4*(2*h/D)^3;
        
        % Kts for Shear
        c1s = .905 + .783*sqrt(h/r) - .075*h/r;
        c2s = -.437 - 1.969*sqrt(h/r) + .553*h/r;
        c3s = 1.557 + 1.073*sqrt(h/r) - .578*h/r;
        c4s = -1.061 + .171*sqrt(h/r) + .086*h/r;
        Kts = c1s + c2s*(2*h/D) + c3s*(2*h/D)^2 + c4s*(2*h/D)^3;

        %% Finding Kf and Kfs
        neubBend = .246 - (3.08e-3)*Sut + (1.51e-5)*Sut^2 - (2.67e-8)*Sut^3;
        neubTors = .190 - (2.51e-3)*Sut + (1.35e-5)*Sut^2 - (2.67e-8)*Sut^3;
        Kf = 1 + ((Kt-1)/(1+neubBend/sqrt(r)));
        Kfs = 1 + ((Kts-1)/(1+neubTors/sqrt(r)));

        %% Find Knockdowns

        Ka = a1*(Sut^b1);
        if d <= 2
            Kb = (d/.3)^(-.107);
        else
            Kb = .91*d^(-.157);
        end
        Kc = 1;
        Kd = 1;
        Ke = 1;                  % Currently for 50% (Table 6-5, p.323)        

        Se = SePrime*Ka*Kb*Kc*Kd*Ke;

        %% Demand
        Ma = Pa*L/1000;
        siga = Ma*(d/2)/I;
        sigm = siga*(Pm/Pa);
        
        Ta = Pa*a/1000;
        sheara = Ta*(d/2)/J;
        shearm = sheara*(Pm/Pa);
        

        %% Factors of Safety
        sigaprime = sqrt((Kf*siga)^2 + (3*(Kfs*sheara)^2));
        sigmprime = sqrt((sigm)^2 + (3*(shearm)^2));
        
        FSlang = Sy/(sigaprime + sigmprime);
        FSgood = ((sigaprime/Se) + (sigmprime/Sut))^(-1);
        if FSlang > FSgood
            FS(i,j) = FSgood;
        else
            FS(i,j) = FSlang;
        end
    
    end
end

%% Plot
plot(rratio,FS, 'LineWidth', 2)
title('Factor of Safety vs Radius Ratio ');
xlabel('Radius Ratio');
ylabel('Factor of Safety');
legend('D/d=1.25','D/d=1.5','D/d=1.75','D/d=2','D/d=2.25','D/d=2.5');



