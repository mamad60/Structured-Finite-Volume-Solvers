%% This Scripts Solves 1D Steady Diffussion Equation
% d/dx(rho*u*FI)=d/dx(Gamma*dphi/dx)+S,
%S=Sp*FI+Su is the source term per unit volume
% First Order Convection Terms Discritization Schemes
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
N=20; % No. CVs in X direcition
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
B=zeros(N,3); %Auxillary Matrix for Constructing Sparse Marix
%==========================================================================
%Problem Setup
Setup
%Grid Generation
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
switch IC  
    case 0 %CDS
        Sch='Central Difference';
    case 1 %First Order Upwin(Donor Cell) Scheme
        Sch='1st Order Upwind';
    case 2 %Exponential Scheme
        Sch='Exponential';
    case 3 %Hybrid Scheme
        Sch='Hybrid';
    case 4 %Pwer Law Scheme
        Sch='Power Law';
end
%% =========================================================================
%Construct Coe. for the Internal cells
for i=2:NIM
    %Compute Fluxes
    [Fe(i),Fw(i),De,Dw]=Flux(i);
    %Compute Peclet numbers
    Pe=Fe(i)/De;
    Pw=Fw(i)/Dw;
    %===============================================================
    %Compute Coe.
    aE(i)=(De*CSch(Pe)+max(-Fe(i),0))*Se(i);  %ae
    aW(i)=(Dw*CSch(Pw)+max(Fw(i),0))*Sw(i);   %aw
end
%-----Apply Boundary Conditions
switch Bw   %West Boundary(Left)
    case 0
        Sp(2)=Sp(2)-aW(2);
        Su(2)=Su(2)+aW(2)*FIw;
        aW(2)=0;
    case 1
        Sp(2)=Sp(2)+0;
        Su(2)=Su(2)+Sw(2)*qw;
        aW(2)=0;
end
switch Be   %East Boundary(Right)
    case 0
        Sp(NIM)=Sp(NIM)-aE(NIM);
        Su(NIM)=Su(NIM)+aE(NIM)*FIe;
        aE(NIM)=0;
    case 1
        Sp(NIM)=Sp(NIM)+0;
        Su(NIM)=Su(NIM)+Se(NIM)*qe;
        aE(NIM)=0;
end

%------Construct aP
aP=aW+aE+(Fe-Fw)-Sp;

%==========================================================================
%Assemble Matrix & Solve
aP=aP(2:end);
aE=aE(2:end-1);
aW=aW(3:end);
% %Sparse Matrix Construction
B(1:N-1,1)=-aW;
B(:,2)=aP;
B(2:N,3)=-aE;
A=spdiags(B,[-1 0 1],N,N);
% A=diag(aP)-diag(aE,1)-diag(aW,-1);
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
    fprintf('\nDiscritisation Scheme= %s',Sch);
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
text(0.03,0.7,strcat('Scheme=',Sch),'Units','Normalized','Edge','blue')
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








