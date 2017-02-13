function SetupFin
%Setups problem for a specific FV problem
% Customize According the problem solved
global Xmin Xmax EXX Gamma0 gStat S0 Sfunc sStat Su0 Sp0 SU SP 
global Bw Be FIw FIe qw qe n2 Tinf
%Problem: Convection in a Cylinderal Fin with cross Section A
%Base is kept at T=100(C) and tip is insulated 
%==========================================================================
%Enter Domain Limits Here
Xmin=0; %Minimum X
Xmax=1; %Maximum X
EXX=1; %Expansion Factor ---X Direction

%=========================================================================
%Enter Cross Section Data Here
D=0.1; % Diameter of the Cylinderical Fin
S0=pi*D^2/4;  %Base Cross-Section---A=1 in y-z Plane
Sfunc=@(S0,x) (S0); %Funtion for Calculation Cross-Sections

%==========================================================================
%Enter Material Properties Here
Gamma0=1; %Diffussion Coefficient
gStat=1; %1: Gamma is constant 2:Gamama varies with x

%==========================================================================
%Enter Source Terms Here
sStat=2; %1: No source Term 2: Constant Source Tem 3:Source Term is a function of  x
n2=25; %(m^-2) (h*P/k*a)
Tinf=20; %Amibeint Temp.(C)
Su0=n2*Tinf; %Base Su
Sp0=-n2; %Base Sp
SU=@(Su0,x) (0);  % Enter the Exprission of Su(x) Here
SP=@(Sp0,x) (0);  % Enter the Exprission of Sp(x) Here

%**************************************************************************
%Modify Boundary Conditions below
Bw=0; % 0:Fixed FI 1:Flux Bc
Be=1; % 0:Fixed FI 1:Flux Bc
FIw=100; %Base is kept @ T=100(C)  % FI at The Left Boundary, Applies only if NL=0
FIe=0;   % FI at The Right Boundary, Applies only if NR=0
qw=0; % Flux at The Left Boundary, Applies only if NL=1
qe=0; % Tip is Insulated % Flux at The Right Boundary, Applies only if NR=1
%**************************************************************************

end

