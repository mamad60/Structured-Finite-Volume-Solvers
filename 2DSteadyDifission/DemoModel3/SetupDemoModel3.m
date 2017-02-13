function SetupDemoModel3
%Setups problem for a specific FV problem
% Customize According the problem solved

% DemoModel3  2D heat conduction with Dirichlet BC and nonuniform conductivity
% and source term. 
%kratio = ratio of thermal conductivity in the core region to the
% thermal conductivity in outer region.  Default:  kratio = 4

global Xmin Xmax EXX Gamma0 gStat sStat Su0 Sp0  kratio
global Bw Be Bn Bs FIw FIe qw qe Ymin Ymax EXY  BCw BCe BCn BCs

%==========================================================================
%Enter Domain Limits Here
Xmin=0; %Minimum X
Xmax=1; %Maximum X
EXX=1; %Expansion Factor ---X Direction
Ymin=0; %Minimum X
Ymax=1; %Maximum X
EXY=1; %Expansion Factor ---X Direction
%==========================================================================
%Enter Material Properties Here
kratio =4; %ratio of thermal conductivity
Gamma0=1; %Diffussion Coefficient
gStat=2; %1: Gamma is constant 2:Gamama varies with x/y
%For this problem DIF file is dirctly modified
% G=@(x,y) (Gamma0);  % Enter the Exprission of Diffussion Coe. Here
%==========================================================================
%Enter Source Terms Here
sStat=3; %1: No source Term 2: Constant Source Tem 3:Source Term is a function of  x/y
Su0=1000; %Base Su
Sp0=0; %Base Sp
%For this problem Source file is dirctly modified
% SU=@(x,y) (1);  % Enter the Exprission of Su(x,y) Here
% SP=@(x,y) (0);  % Enter the Exprission of Sp(x,y) Here
%**************************************************************************
%Modify Boundary Conditions below
%----------Chose Boundary Type
Bw=0; % 0:Fixed FI 1:Flux Bc --->West Boundary
Be=0; % 0:Fixed FI 1:Flux Bc--->East Boundary
Bn=0; % 0:Fixed FI 1:Flux Bc --->North Boundary
Bs=0; % 0:Fixed FI 1:Flux Bc--->South Boundary
%----------Dirichlet BCs
FIw=20;   % FI at The Left Boundary, Applies only if Bw=0
FIe=0;   % FI at The Right Boundary, Applies only if Be=0
FIn=0;   % FI at The Left Boundary, Applies only if Bn=0
FIs=10;   % FI at The Right Boundary, Applies only if Bs=0
%--------Neumann Bcs
qw=0; % Flux at The Left Boundary, Applies only if Bw=1
qe=0; % Flux at The Right Boundary, Applies only if Be=1
qn=0; % Flux at The Left Boundary, Applies only if Bn=1
qs=0; % Flux at The Right Boundary, Applies only if Bs=1
%**************************************************************************
%---BCw BCe BCn BCs are 1D Arrays containing BCs for the West,East,North and Top Boundaries, respectively
%Modify Exprission below if BCs vary with x or y
%==========================================================================
%         j=1:NJ;
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
switch Bs
    case 0
        BCs(:)=FIs;
    case 1
        BCs(:)=qs;
end
%Add other Boundaries if Applies Here
%==========================================================================
end

