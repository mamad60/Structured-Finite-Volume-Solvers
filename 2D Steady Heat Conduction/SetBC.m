function  SetBC
%In this Function Boundary Conditons are set and Modified
% Change according to the problem
%Global Variables Definition
global X Y  m n 
global NL NR NB NT BCl BCr BCb BCt Tf
%-----Boundary Conditions Input Setting
NL=1; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Left Boundary
NR=1; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Right Boundary
NB=0; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Bottom Boundary
NT=0; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Top Boundary

Tf=273.0; %Temperature @ infinity used if NL=2 Or NR=2 q=h(Tf-TB)
%-----Varibles Definition
BCl=zeros(n,1);
BCr=zeros(n,1);
BCt=zeros(n,1);
BCb=zeros(n,1);
%---Define Auxiallry Variables Here
Tl=273; % Temperature at The Left Boundary, Applies only if NL=0
Tr=283; % Temperature at The Right Boundary, Applies only if NR=0
Tb=273; % Temperature at The Bottom Boundary, Applies only if NB=0
Tt=283; % Temperature at The Top Boundary, Applies only if NT=0
Ql=0; % Flux at The Left Boundary, Applies only if NL=1
Qr=0; % Flux at The Right Boundary, Applies only if NR=1
Qb=0; % Flux at The Bottom Boundary, Applies only if NB=1
Qt=0; % Flux at The Top Boundary, Applies only if NT=1
%Flux Normal to boundary & Towards it Assumed Positive sign
hl=0; % Convection Coefficient at The Left Boundary, Applies only if NL=2
hr=0.1; % Flux at The Left Boundary, Applies only if NL=2
hb=0; % Convection Coefficient at The Bottom Boundary, Applies only if NB=2
ht=0.1; % Flux at The Top Boundary, Applies only if NT=2

%Set Boundary conditons

%Left Boundary
switch NL  %----Enter Exprission for BCs under the appropriat entry
    case 0
        BCl(:)=Tl;
    case 1
        BCl(:)=Ql;
    case 2
        BCl(:)=hl;
end

%Right  Boundary
switch NR  %----Enter Exprission for BCs under the appropriat entry
    case 0
        BCr(:)=Tr;
    case 1
        BCr(:)=Qr;
    case 2
        BCr(:)=hr;
end
%Top Boundary
switch NT  %----Enter Exprission for BCs under the appropriat entry
    case 0
        BCt(:)=Tt;
    case 1
        BCt(:)=Qt;
    case 2
        BCt(:)=ht;
end
%Bottom Boundary
switch NB  %----Enter Exprission for BCs under the appropriat entry
    case 0
        BCb(:)=Tb;
    case 1
        BCb(:)=Qb;
    case 2
        BCb(:)=hb;
end
%Add other Boundaries if Applies Here
end

