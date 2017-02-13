function Setup
%Setups problem for a specific FV problem
% Customize According the problem solved
global Xmin Xmax EXX Gamma0 gStat sStat Su0 Sp0 SU SP G
global Bw Be Bn Bs FIw FIe qw qe Ymin Ymax EXY  BCw BCe BCn BCs

%==========================================================================
%Enter Material Properties Here
Gamma0=1; %Diffussion Coefficient
gStat=1; %1: Gamma is constant 2:Gamama varies with x
G=@(x,y) (Gamma0);  % Enter the Exprission of Diffussion Coe. Here

%==========================================================================
%Enter Source Terms Here
sStat=1; %1: No source Term 2: Constant Source Tem 3:Source Term is a function of  x
Su0=1; %Base Su
Sp0=1; %Base Sp
SU=@(x,y) (0);  % Enter the Exprission of Su(x) Here
SP=@(x,y) (0);  % Enter the Exprission of Sp(x) Here
%**************************************************************************
%Modify Boundary Conditions below
%----------Chose Boundary Type
Bw=0; % 0:Fixed FI 1:Flux Bc --->West Boundary
Be=0; % 0:Fixed FI 1:Flux Bc--->East Boundary
Bn=1; % 0:Fixed FI 1:Flux Bc --->North Boundary
Bs=1; % 0:Fixed FI 1:Flux Bc--->South Boundary
%----------Dirichlet BCs
FIw=0;   % FI at The Left Boundary, Applies only if Bw=0
FIe=100;   % FI at The Right Boundary, Applies only if Be=0
FIn=0;   % FI at The Left Boundary, Applies only if Bn=0
FIs=100;   % FI at The Right Boundary, Applies only if Bs=0
%--------Neumann Bcs
qw=0; % Flux at The Left Boundary, Applies only if Bw=1
qe=0; % Flux at The Right Boundary, Applies only if Be=1
qn=0; % Flux at The Left Boundary, Applies only if Bn=1
qs=0; % Flux at The Right Boundary, Applies only if Bs=1
%**************************************************************************
%---BCw BCe BCn BCs are 1D Arrays containing BCs for the West,East,North and Top Boundaries, respectively
%Modify Exprission below if BCs vary with x or y
%==========================================================================
%         j=1:M;
%Left Boundary
switch Bw
    case 0
        BCw(:)=FIw;%-X(j);end
    case 1
        BCw(:)=qw;
end
%Right  Boundary
switch Be
    case 0
        BCe(:)=FIe;
    case 1
        BCe(:)=qe;
end
%i=1:NI;
%Top Boundary
switch Bn
    case 0
        BCn(:)=FIn;
    case 1
        BCn(:)=qn;
end
%Bottom Boundary
%i=1:N;
switch Bs
    case 0
        BCs(:)=FIs;
    case 1
        BCs(:)=qs;
end
%Add other Boundaries if Applies Here
%==========================================================================
end

