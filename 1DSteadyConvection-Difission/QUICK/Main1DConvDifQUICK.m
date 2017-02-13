%% This Scripts Solves 1D Steady Diffussion Equation
% d/dx(rho*u*FI)=d/dx(Gamma*dphi/dx)+S,
%S=Sp*FI+Su is the source term per unit volume
%Implementation Of QUICK
%Discritization Stencil   |--->  W---P---E
%By Mohammad Aghakhani, 2014,for Teaching purposes on the CFD Class
clc
clear
close all

global X gStat XC NIM NI Xmax Xmin EXX sStat P
global DX rho u IC 
global Xw Xe Gammae Gammaw
global Bw Be FIw FIe qw qe



%--Inputs
N=5; % No. CVs in X direcition
%=====================================================
%Variable Allocation
NI=N+2;  % Two Extra points for first & last Cvs(Boundary)
NIM=NI-1; % No CVS +1, the NIM CV is on the boundary
%----------------------
X=zeros(1,NI);  %Grid Nodes
XC=zeros(1,NI); %Cell Centers
FI=zeros(1,NI); %Solution Variavle
aE=zeros(1,NIM);  %East Coe.
aW=zeros(1,NIM);  %West Coe.
aP=zeros(1,NIM);  %Point Coe.
Su=zeros(NIM,1);  % Source Term
Sp=zeros(1,NIM);  % Source Term
Xw=zeros(1,NI);
Xe=zeros(1,NI);
Gammaw=zeros(1,NI);
Gammae=zeros(1,NI);
Fe=zeros(1,NIM);
Fw=zeros(1,NIM);
Sw=zeros(1,NI);
Se=zeros(1,NI);
F=zeros(N,1); %Solution in the interior cells
B=zeros(N,5); %Auxillary Matrix for Constructing Sparse Marix
%==========================================================================
%Problem Setup
Setup
%Grid Generation
EXX=1; %QUICK On uniform Grid
[X,XC,Xw,Xe,DX] = Grid1d(Xmin,Xmax,N,EXX);
%Compute Diffussion Coe
[Gammae Gammaw] = BiHarmonic;
%Compute & Print Peclect Number
Gamma=min(Gammae);
P=rho*max(u)*(Xmax-Xmin)/Gamma; %Peclet Number
fprintf('\nThe Peclet Number is %2.2f\n',P)
%Compute Cross-Sections
[Se Sw]=CrossSection;   %Cross-Section at boundaries
%Compute Source Terms
[Su Sp]=Source(sStat);
%-----------------Store Discritisation Scheme
%% ========================================================================

%Construct Coe. for the Internal cells
for i=4:NIM-1  %If u>0, 2 cells from west boundary and last cell are excluded
    %Compute Fluxes
    [Fe(i),Fw(i),De,Dw]=Flux(i);
    %set ae & aw--QUICK Scheme
    if Fw(i)>0
        aw=1;
    else
        aw=0;
    end
    if Fe(i)>0
        ae=1;
    else
        ae=0;
    end
    %===============================================================
    %Compute Coe.
    aE(i)=(De-3/8*ae*Fe(i)-6/8*(1-ae)*Fe(i)-1/8*(1-aw)*Fw(i));  %ae
    aW(i)=(Dw+6/8*aw*Fw(i)+1/8*ae*Fe(i)+3/8*(1-aw)*Fw(i));   %aw
    aEE(i)=1/8*(1-ae)*Fe(i); %aEE
    aWW(i)=-1/8*aw*Fw(i);   %aWW
end
%

%-----Apply Boundary Conditions

%West Boundary(Left) i==2
[Fe(2),Fw(2),De,Dw]=Flux(2);
aW(2)=0;
aWW(2)=0;
aEE(2)=0;
aE(2)=De+1/3*Dw-3/8*Fe(2);
Sp(2)=-(8/3*Dw+2/8*Fe(2)+Fw(2));
Su(2)=(8/3*Dw+2/8*Fe(2)+Fw(2))*FIw;
aW(2)=0;
%West Boundary(Left) i=3
[Fe(3),Fw(3),De,Dw]=Flux(3);
aWW(3)=0;
aEE(3)=0;
aW(3)=Dw+7/8*Fw(3)+1/8*Fe(3);
aE(3)=De-3/8*Fe(3);
Sp(2)=1/4*Fw(3);
Su(2)=-(1/4*Fw(3))*FIw;
%East Boundary(Right) i=NIM
[Fe(NIM),Fw(NIM),De,Dw]=Flux(NIM);
aWW(NIM)=-1/8*Fw(NIM);
aW(NIM)=Dw+1/3*De+6/8*Fw(NIM);
aE(NIM)=0;
aEE(NIM)=0;
Sp(NIM)=-(8/3*De-Fe(NIM));
Su(NIM)=(8/3*De-Fe(NIM))*FIe;
aE(NIM)=0;

%------Construct aP
aP=aW+aE+aWW+aEE+(Fe-Fw)-Sp;

%==========================================================================
%Assemble Matrix & Solve
aP=aP(2:end);
aE=aE(2:end-1);
aW=aW(3:end);
aWW=aWW(4:end);
aEE=aEE(2:end-2);
% Sparse Matrix Construction
A=diag(aP)-diag(aE,1)-diag(aW,-1)-diag(aWW,-2)-diag(aEE,2);
b=Su(2:end);
F=A\b;
%% =========================================================================
%Post-Process
%For Dirichlet BCs and no source term equations has and exact solution
i=1:NI;
FIExact(i)=(exp(P*XC(i)/(Xmax-Xmin))-1)/(exp(P)-1)*(FIe-FIw)+FIw;
%----Fill in Solution Vector
FI(2:NIM)=F;
%----Interpolate FI @ Boundaries,IF Neumann BC is imposed
%-------West Boundary 
switch Bw
    case 0
        FI(1)=FIw;
    case 1
        FI(1)=qe*(XC(1)-XC(2))+FI(2);
end
%-------East Boundary
switch Be
    case 0
        FI(NI)=FIe;
    case 1
        FI(NI)=qw*(XC(NI)-XC(NIM))+FI(NIM);
end
% %Plot Result
if Bw==0 && Be==0
    figure
    plot(XC,FIExact,'-g',XC,FI,'o:b','LineWidth',1.5)
    xlabel('X')
    ylabel('\phi')
    xlim([Xmin Xmax])
    ylim([min(FIExact) max(FIExact)])
    title('Solution Profile')
    legend('Exact','Numerical','Location','NorthWest')
    fprintf('\nThe Norm of Error is:\t%2.2e\n',norm(FI-FIExact))
else
    plot(XC,FI,'LineWidth',1.5)
    xlabel('X')
    ylabel('\phi')
    xlim([Xmin Xmax])
    ylim([min(FI) max(FI)])
    title('Solution Profile')
end
text(0.03,0.8,strcat('Peclet=',num2str(P)),'Units','Normalized','Edge','red')    
text(0.03,0.6,strcat('N=',num2str(N)),'Units','Normalized','Edge','blue')    


% %Report Fluxes
%  qwD=Gammae(2)*(FI(1)-FI(2))/abs((XC(1)-XC(2)));
%  qeD=Gammaw(NIM)*(FI(NI)-FI(NIM))/abs((XC(NI)-XC(NIM)));
%  fprintf('\nDiffuse Flux of FI @ the West Boundary is:\t%2.4f',qwD);
%  fprintf('\nDiffuse Flux of FI @ the East Boundary is:\t%2.4f',qeD);
%  qwC=rho*u(2);
%  qeC=rho*u(NIM);
%  fprintf('\nConvective Flux of FI @ the West Boundary is:\t%2.4f',qwC);
%  fprintf('\nConvective Flux of FI @ the East Boundary is:\t%2.4f',qeC);
%  fprintf('\nFlux Inbalance @ the Boundaries is:\t%2.4f\n',qeD-qwD+qeC-qwC); 
disp('Good Lock, Mohammad Aghakhani')








