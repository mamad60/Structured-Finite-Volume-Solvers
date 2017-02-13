function  SetBC
%In this Function Boundary Conditons are set and Modified
% Change according to the problem
%Global Variables Definition
global X Y  XC YC N M NI NJ NIM NJM NIJ NM FI
global NL NR NB NT BCl BCr BCb BCt

%====Boundry Conditions Input Setting
NL=1; % 0:Fixed FI 1:Flux Bc
NR=1; % 0:Fixed FI 1:Flux Bc
NB=0; % 0:Fixed FI 1:Flux Bc
NT=0; % 0:Fixed FI 1:Flux Bc
FIl=100;   % FI at The Left Boundary, Applies only if NL=0
FIr=0;   % FI at The Right Boundary, Applies only if NR=0
FIb=0; % FI at The Bottom Boundary, Applies only if NB=0
FIt=100; % FI at The Top Boundary, Applies only if NT=0
Ql=0; % Flux at The Left Boundary, Applies only if NL=1
Qr=0; % Flux at The Right Boundary, Applies only if NR=1
Qb=0; % Flux at The Bottom Boundary, Applies only if NB=1
Qt=0; % Flux at The Top Boundary, Applies only if NT=1
%==========================================================================
BCl=zeros(NJ,1);
BCr=zeros(NJ,1);
BCt=zeros(1,NI);
BCb=zeros(1,NI);
%==========================================================================
%Set Boundary conditons
%Left Boundary
switch NL  %----Enter Exprission for BCs under the appropriat entry
    case 0
        %         j=1:NJ;
     BCl(:)=FIl;%-X(j);end
            
    case 1
        BCl(:)=Ql;
end

%Right  Boundary
switch NR  %----Enter Exprission for BCs under the appropriat entry
    case 0
        BCr(:)=FIr;
    case 1
        BCr(:)=Qr;
end
%Top Boundary
switch NT  %----Enter Exprission for BCs under the appropriat entry
    case 0
        BCt(:)=FIt;
    case 1
        BCt(:)=Qt;
end
%Bottom Boundary
switch NB  %----Enter Exprission for BCs under the appropriat entry
    case 0
        BCb(:)=FIb;
    case 1
        BCb(:)=Qb;
end
%Add other Boundaries if Applies Here
%==========================================================================
%Fill in FI Rows with Boundary Conditions
if NL==0,FI(1,:)=BCl; end    %Left BC
if NR==0,
    FI(NI,:)=BCr;
end %Right BC
% Top & Bottom
if NB==0,FI(:,1)=BCb; end  %Bottom Boundary
if NT==0,
    FI(:,NJ)=BCt;
end %Top Boundary
end

