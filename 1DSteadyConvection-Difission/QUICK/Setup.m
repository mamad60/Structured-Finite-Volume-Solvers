function Setup
%Setups problem for a specific FV problem
% Customize According the problem solved
global Xmin Xmax EXX Gamma0 gStat S0 Sfunc sStat Su0 Sp0 SU SP NI
global Bw Be FIw FIe qw qe
global rho u IC

%==========================================================================
%Enter Domain Limits Here
Xmin=0; %Minimum X
Xmax=1; %Maximum X
EXX=1; %Expansion Factor ---X Direction
%==========================================================================
%Material Properties and Velocities
rho=1; %Density
u=1*ones(1,NI); %Known Velociy Field
%=========================================================================
%Enter Cross Section Data Here
S0=1;  %Base Cross-Section
Sfunc=@(S0,x) (1); %Funtion for Calculation Cross-Sections
%==========================================================================
%Enter Material Properties Here
Gamma0=0.5; %Diffussion Coefficient
gStat=1; %1: Gamma is constant 2:Gamama varies with x

%==========================================================================
%Enter Source Terms Here
sStat=1; %1: No source Term 2: Constant Source Tem 3:Source Term is a function of  x
Su0=1; %Base Su
Sp0=1; %Base Sp
SU=@(Su0,x) (0);  % Enter the Exprission of Su(x) Here
SP=@(Sp0,x) (0);  % Enter the Exprission of Sp(x) Here

%**************************************************************************
%Modify Boundary Conditions below
Bw=0; % 0:Fixed FI 1:Flux Bc
Be=0; % 0:Fixed FI 1:Flux Bc
FIw=1;   % FI at The Left Boundary, Applies only if NL=0
FIe=0;   % FI at The Right Boundary, Applies only if NR=0
qw=1; % Flux at The Left Boundary, Applies only if NL=1
qe=1; % Flux at The Right Boundary, Applies only if NR=1
%**************************************************************************

end

