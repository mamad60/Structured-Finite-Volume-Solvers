%% This Scripts Solves 1D Convection-Diffusion Equation
% d/dx(rho*u*FI)=d/dx(Gamma*dphi/dx)+S, From Patankar Book, Section 4-2
%Demonstates the application of the direct & iterative Matrix Solvers---Single Grid Version
%Discritization Stencil   |--->  W---P---E
%By Mohammad Aghakhani, 2010,for Teaching purposes on the CFD Class
%Feel Free using for Educaitonal & Research Purposes

clc
clear
close all

%Please set the Folder for Iterative Solvers below or Copy files to the Current directry
iFolder='../../Iterative Solvers';
%Check & Set Path for Solver Directery
pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
    onPath = any(strcmpi(iFolder, pathCell));
else
    onPath = any(strcmp(iFolder, pathCell));
end
if ~onPath,addpath(genpath(iFolder)); end

%=====Inputs
N=50; % No. CVs in X direcition
Xmin=0; %Minimum X
Xmax=1; %Maximum X
EXX=1; %Expansion Factor ---X Direction
rho=10; %Density
Gamma=1; %Diffussion Coefficient
IC=2; %Convective Flux Discritization Scheme:
%0:Centeral Difference 1:First Order Upwind  
%2:Exponential  3:Hybrid 4:Power Law
% u=zeros(N+2,1);
u=1; %Known Velociy Field
%-----------------------------------------------------
iMethod=8; % Solution Method:
%0 MATLAB Direct Solver 1: TDMA(Direct)
%2:Jacobi(Matrix) 3:Jacobi 4:Gauss-Seidel(Matrix) 5:Gauss-Seidel
%6: Symmetric Gauss-Seidel Matrix 7:Symmetric Gauss-Seidel
%8:SOR(Matrix) 9:SOR 10:Symmetric SOR Matrix
%11:Red-Black Gauss-Seidel  12:SOR Solver 4 Triangular matrix
%(13:Conjugate Gradient   14: Preconditioned Conjugate Gradient)---> Use
%this solve with causion for ,for Convection-Diffission Matrix is not
%Symmetric,Specially for high Peclet Number
%15: BiConjugate gradient
espSolver=1e-6; % Maximum Residual of Matrix Iterative Solver
MaxITSolver=100000; % Maximum Iteraton of Matrix Iterative Solver
%-----------------------------------------------------
%-----Boundary Conditions
NL=0; % 0:Fixed FI 1:Flux Bc 
NR=1; % 0:Fixed FI 1:Flux Bc 
FIl=0; % FI at The Left Boundary, Applies only if NL=0
FIr=100; % FI at The Right Boundary, Applies only if NR=0
Ql=-100; % Flux at The Left Boundary, Applies only if NL=1
Qr=1000; % Flux at The Right Boundary, Applies only if NR=1
%Flux Normal to boundary & Towards it Assumed Positive sign
%=====================================================
%Variable Allocation
NI=N+2;  % Two Extra points for first & last Cvs(Boundary)
NIM=NI-1; % No CVS +1, the NIM CV is on the boundary
%----------------------
X=zeros(1,NI);  %Grid Nodes
XC=zeros(1,NI); %Cell Centers
FI=zeros(1,NI); %Solution Variavle
Q=zeros(1,NI);  % Source Term
AE=zeros(1,NIM);  %East Coe.
AW=zeros(1,NIM);  %West Coe.
AP=zeros(1,NIM);  %Point Coe.
Q=zeros(1,NIM);  % Source Term
Xw=zeros(1,NI);
Xe=zeros(1,NI);
T=zeros(N,1); %Solution in the interior cells
B=zeros(N,3); %Auxillary Matrix for Constructing Sparse Marix
%==========================================================================
%Grid Generation
[X,XC,Xw,Xe,DX] = Grid1d(Xmin,Xmax,N,EXX);

P=rho*u*(Xmax-Xmin)/Gamma; %Peclet Number
fprintf('\nThe Peclet Number is %2.2f\n',P)
%==============Fill Boundary Values
%--------Left & Right
if NL==0,FI(1)=FIl; end    %Left BC
if NR==0, FI(NI)=FIr; end %Right BC
T(:)=0;   %Field Initialazation


%% ========================================================================
tic;
F=rho*u; %DEnsity- Velocity Product, Flux of u
%Calculate Coefficients & Assemble Solution Matrix AT=b
%-----Interior Cells
for i=2:NIM
    Fe=F; %Convection Coe.
    De=Gamma/Xe(i); %Diffission COe.
    Fw=F;
    Dw=Gamma/Xw(i);
    switch IC  % Discritize According to the choosen Scheme
        case 0 %CDS
            AE(i)=De-Fe/2;     %ae
            AW(i)=Dw+Fw/2;     %aw
            Sch='Central Difference';
        case 1 %First Order Upwin(Donor Cell) Scheme
            AE(i)=De+max(-Fe,0);   %ae
            AW(i)=Dw+max(Fw,0);    %aw
            Sch='1st Order Upwind';
        case 2 %Exponential Scheme
            AE(i)=Fe/(exp(Fe/De)-1);               %ae
            AW(i)=Fw*exp(Fw/Dw)/(exp(Fw/Dw)-1);    %aw
            Sch='Exponential';
        case 3 %Hybrid Scheme
            AE(i)=max([-Fe,De-Fe/2,0]);   %ae
            AW(i)=max([Fw,Dw+Fw/2,0]);    %aw
            Sch='Hybrid';
        case 4 %Pwer Law Scheme
            AE(i)=De*max(0,(1-0.1*abs(Fe)/De)^5)+max(0,-Fe);   %ae
            AW(i)=Dw*max(0,(1-0.1*abs(Fw)/Dw)^5)+max(0,Fw);    %aw         
            Sch='Power Law';          
    end
    %ap
    AP(i)=AE(i)+AW(i)+(Fe-Fw);
    %b
    Q(i)=0;
end
%==========================================================================
%----Apply BCs by invoking into source Terms
%Left Boundary
switch NL
    case 0  % Fixed vakue of FI
        Q(2)=Q(2)+AW(2)*FI(1);  %BCs at first CV
        AW(2)=0;
    case 1  %Fixed Flux
        Q(2)=Q(2)+Ql;  %BCs at first CV
        AP(2)=AP(2)-AW(2);
        AW(2)=0;
end

%Right Boundary
switch NR
     case 0  % Fixed vakue of FI
        Q(NIM)=Q(NIM)+AE(NIM)*FI(NI);  %BCs at first CV
        AE(NIM)=0;
    case 1  %Fixed Flux
        Q(NIM)=Q(NIM)+Qr;  %BCs at first CV
        AP(NIM)=AP(NIM)-AE(NIM);
        AE(NIM)=0;
end
%==========================================================================
%Assemble Matrix & Solve
AP=AP(2:end);
AE=AE(2:end-1);
AW=AW(3:end);
%Sparse Matrix Construction
B(1:N-1,1)=-AW;
B(:,2)=AP;
B(2:N,3)=-AE;
A=spdiags(B,[-1 0 1],N,N);
%A=diag(AP)-diag(AE,1)-diag(AW,-1);
b=Q(2:NIM)';
%Solve Equations
T=Solve(A,T,b,MaxITSolver,espSolver,iMethod); %Solution in the Interior Cells
toc

%% ==========================================================================
% -----Post-Processing
FI(2:NIM)=T(:);
%For Dirichlet BCs and no source tem equations has and exact solution
i=1:NI;
P=rho*u*(Xmax-Xmin)/Gamma; %Peclet Number
FIExact(i)=(exp(P*XC(i)/(Xmax-Xmin))-1)/(exp(P)-1)*(FIr-FIl)+FIl;
%Interpolate FI @ Boundaries
%Left Boundary Interpolation, IF q is applied
if NL==1
    p=polyfit([XC(2) XC(3) XC(4)],[FI(2) FI(3) FI(4)],2);
    l=polyval(p,XC(1));
    FI(1)=l;
end
%Right Boundary
if NR==1
    p=polyfit([XC(NIM) XC(NIM-1) XC(NIM-2)],[FI(NIM) FI(NIM-1) FI(NIM-2)],2);
    r=polyval(p,XC(NI));
    FI(NI)=r;
end
%Plot Result
if NL==0 && NR==0
    figure
    plot(XC,FIExact,'-g',XC,FI,'ob','LineWidth',1.5)
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

% %Report Fluxes
%  ql=rho*u*FI(2)-Gamma*(FI(2)-FI(1))/(XC(2)-XC(1));
%  qr=rho*u*FI(NIM)-Gamma*(FI(NI)-FI(NIM))/(XC(NI)-XC(NIM));
%  fprintf('\nFlux of FI @ the Left Boundary is:\t%2.4f',ql);
%  fprintf('\nFlux of FI @ the Right Boundary is:\t%2.4f',qr);
%  fprintf('\nFlux Inbalance @ the Boundaries is:\t%2.4f\n',qr+ql);
 
disp('Good Lock, Mohammad Aghakhani')










