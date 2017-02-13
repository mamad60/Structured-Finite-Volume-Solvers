%% This Scripts Solves 1D Steady Heat Transfer Equation
% d/dx(kdT/dx)+S=0
%Problem: %Problem: A Plate with L=2(cm) and Heat Gerneration ans constant Temp. Bcs
%qdot=1000kW/m^3 Ta=100(C), Tb=200(C)
%S=Sp*FI+Su is the source term per unit volume
% Both Thermal Conductivity and Source therms can be varied with time & x
%Discritization Stencil   |--->  W---P---E
%By Mohammad Aghakhani, 2014,for Teaching purposes on the CFD Class
clc
clear
close all

%Add the Path for the Parent Folder
iFolder='../';
%Check & Set Path for Solver Directery
pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
  onPath = any(strcmpi(iFolder, pathCell));
else
  onPath = any(strcmp(iFolder, pathCell));
end
if ~onPath,addpath(genpath(iFolder)); end

global X gStat XC NIM NI Xmax Xmin EXX sStat
global Bw Be FIw FIe qw qe Gamma0 DX  n2 Tinf

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
Sw=zeros(1,NI);
Se=zeros(1,NI);
F=zeros(N,1); %Solution in the interior cells
B=zeros(N,3); %Auxillary Matrix for Constructing Sparse Marix
%==========================================================================
%Problem Setup
SetupFin
%Grid Generation
[X,XC,Xw,Xe,DX] = Grid1d(Xmin,Xmax,N,EXX);
%Compute Diffussion Coe
[Gammae Gammaw] = BiHarmonic;
%Compute Cross-Sections
[Se Sw]=CrossSection;   %Cross-Section at boundaries
%Compute Source Terms
[Su Sp]=Source(sStat);
%%==========================================================================
%Construct Coe. for the Internal cells
for i=2:NIM
    aW(i)=Gammaw(i)*Sw(i)/Xw(i);
    aE(i)=Gammae(i)*Se(i)/Xe(i);
%     aP(i)=aW(i)+aE(i)-Sp(i);
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
aP=aW+aE-Sp;

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
%%=========================================================================
%Post-Process
%----Fill in Solution Vector
FI(2:NIM)=F;
%----Interpolate FI @ Boundaries
%--------Left Boundary Interpolation, IF q is applied
switch Bw
    case 0
        FI(1)=FIw;
    case 1
        p=polyfit([XC(2) XC(3) XC(4)],[FI(2) FI(3) FI(4)],2);
        l=polyval(p,XC(1));
        FI(1)=l;
end
%--------Right Boundary
switch Be
    case 0
        FI(NI)=FIe;
    case 1
        FI(NI)=FI(NIM); % Zero Gradient
        
end

%Compute The Exact Solution & Compare Result
i=1:NI;
L=(Xmax-Xmin);
n=sqrt(n2);
ThetaExact(i)=(cosh(n*(L-XC(i))))/cosh(n*L);
Theta(i)=(FI-Tinf)/(FIw-Tinf);
figure %Comparison with the Exact Solution-Theta
plot(XC,ThetaExact,'-b',XC,Theta,'og','LineWidth',1.5)
xlabel('X((m)')
ylabel('\theta')
xlim([Xmin Xmax])
title('Solution Profile')
legend('Exact','Numerical','Location','NorthEast')
text(0.25,0.95,'Cylinderical Fin with Convection ','Units','Normalized','Edge','blue')
text(0.4,0.9,strcat('N=',num2str(N)),'Units','Normalized','Edge','black')

figure %Comparison with the Exact Solution-Temperature
TExact=ThetaExact*(FIw-Tinf)+Tinf;
plot(XC,TExact,'-b',XC,FI,'og','LineWidth',1.5)
xlabel('X((m)')
ylabel('Temperature(^oC)')
xlim([Xmin Xmax])
title('Solution Profile')
legend('Exact','Numerical','Location','NorthEast')
text(0.25,0.95,'Cylinderical Fin with Convection ','Units','Normalized','Edge','blue')
text(0.4,0.9,strcat('N=',num2str(N)),'Units','Normalized','Edge','black')

figure  %Plot Error
plot(XC(2:NIM),abs(ThetaExact(2:NIM)-Theta(2:NIM))./ThetaExact(2:NIM),'-sr','LineWidth',1.5)
xlabel('X(m)')
ylabel('Error(%)')
xlim([XC(2) XC(NIM)])
title('Error')
%Print Erroe Norm
fprintf('The Norm of Error is %2.2e\n',norm(ThetaExact-Theta))
%Plot Difiussion Coe.
if gStat==2
    figure
    plot(XC(2:NIM),Gammaw(2:NIM),'O-r','LineWidth',1.5)
    xlabel('XC')
    ylabel('\gamma_w')
    xlim([XC(2) XC(NIM)])
    ylim([min(Gammaw(2:NIM)) max(Gammaw(2:NIM))])
    title('Difiussion Coefficient Profie')
end

% %Report Fluxes
 ql=Gammae(2)*(FI(1)-FI(2))/abs((XC(1)-XC(2)));
 qr=Gammaw(NIM)*(FI(NI)-FI(NIM))/abs((XC(NI)-XC(NIM)));
 fprintf('\nFlux of FI @ the Left Boundary is:\t%2.4f',ql);
 fprintf('\nFlux of FI @ the Right Boundary is:\t%2.4f',qr);
 fprintf('\nFlux Inbalance @ the Boundaries is:\t%2.4f\n',qr+ql);
 
disp('Good Lock, Mohammad Aghakhani')








