function SetupFullyDevQuarterRect
%Setups problem for a specific FV problem
% Customize According the problem solved

%DemoFullyDevQuarterRect  Fully-developed laminar flow in a rectangular duct.
%w=0 is imposed on on the north and east boundaries
%Symmetery BC is imposed on the west and south boundaries
global Xmin Xmax EXX Gamma0 gStat sStat Su0 Sp0 SU SP G
global Bw Be Bn Bs FIw FIe qw qe Ymin Ymax EXY  BCw BCe BCn BCs

%==========================================================================
%Enter Domain Limits Here
Xmin=0; %Minimum X
Xmax=0.5; %Maximum X
EXX=1; %Expansion Factor ---X Direction
Ymin=0; %Minimum X
Ymax=0.5; %Maximum X
EXY=1; %Expansion Factor ---X Direction
dpdz = -pi;  mu = sqrt(5);  %  Arbitrary pressure gradient and viscosity values
%==========================================================================
%Enter Material Properties Here
Gamma0=mu; %Diffussion Coefficient
gStat=1; %1: Gamma is constant 2:Gamama varies with x/y
G=@(x,y) (Gamma0);  % Enter the Exprission of Diffussion Coe. Here

%==========================================================================
%Enter Source Terms Here
sStat=2; %1: No source Term 2: Constant Source Tem 3:Source Term is a function of  x
Su0=-dpdz; %Base Su
Sp0=0; %Base Sp
SU=@(Su0,x,y) (0);  % Enter the Exprission of Su(x,y) Here
SP=@(Sp0,x,y) (0);  % Enter the Exprission of Sp(x,y) Here
%**************************************************************************
%Modify Boundary Conditions below
%----------Chose Boundary Type
Bw=1; % 0:Fixed FI 1:Flux Bc --->West Boundary
Be=0; % 0:Fixed FI 1:Flux Bc--->East Boundary
Bn=0; % 0:Fixed FI 1:Flux Bc --->North Boundary
Bs=1; % 0:Fixed FI 1:Flux Bc--->South Boundary
%----------Dirichlet BCs
FIw=0;   % FI at The Left Boundary, Applies only if Bw=0
FIe=0;   % FI at The Right Boundary, Applies only if Be=0
FIn=0;   % FI at The Left Boundary, Applies only if Bn=0
FIs=0;   % FI at The Right Boundary, Applies only if Bs=0
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

