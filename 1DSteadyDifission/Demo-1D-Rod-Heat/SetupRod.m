function SetupRod
%Setups problem for a specific FV problem
% Customize According the problem solved
global Xmin Xmax EXX Gamma0 gStat S0 Sfunc sStat Su0 Sp0 SU SP
global Bw Be FIw FIe qw qe
%Problem:Consider the problem of source-free heat conduction in an insulated
% rod whose ends are maintained at constant temperatures of 100 oC and
% 500 oC respectively. Calculate the steady state temperature in the rod.
% Take thermal conductivity k = 1000 W/mK, cross-sectional area A = 10×10–3
% Domain is considered 0<x<0.5(m)

%==========================================================================
%Enter Domain Limits Here
Xmin=0; %Minimum X
Xmax=0.5; %Maximum X
EXX=1; %Expansion Factor ---X Direction

%=========================================================================
%Enter Cross Section Data Here
S0=10*10^-3;  %Base Cross-Section
Sfunc=@(S0,x) (S0); %Funtion for Calculation Cross-Sections

%==========================================================================
%Enter Material Properties Here
Gamma0=1000;%(k=1000 W/Mk) %Diffussion Coefficient
gStat=1; %1: Gamma is constant 2:Gamama varies with x

%==========================================================================
%Enter Source Terms Here
sStat=1; %1: No source Term 2: Constant Source Tem 3:Source Term is a function of  x
Su0=0; %Base Su
Sp0=0; %Base Sp
SU=@(Su0,x) (0);  % Enter the Exprission of Su(x) Here
SP=@(Sp0,x) (0);  % Enter the Exprission of Sp(x) Here

%**************************************************************************
%Modify Boundary Conditions below
Bw=0; % 0:Fixed FI 1:Flux Bc
Be=0; % 0:Fixed FI 1:Flux Bc
FIw=100;   % FI at The Left Boundary, Applies only if NL=0
FIe=500;   % FI at The Right Boundary, Applies only if NR=0
qw=0; % Flux at The Left Boundary, Applies only if NL=1
qe=0; % Flux at The Right Boundary, Applies only if NR=1
%**************************************************************************

end

