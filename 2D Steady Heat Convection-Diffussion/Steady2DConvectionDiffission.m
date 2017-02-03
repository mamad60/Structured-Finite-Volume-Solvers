%% This Scripts Solves 2D Convection-Diffusion Equation
% d/dx(rho*u*FI)+d/dx(Jx)+d/dy(Jy)+S, From Patankar Book, Section 4-2
%Jx=rho*u*FI-gamma*d/dx(FI),Jy=rho*v*FI-gamma*d/dy(FI)
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
%Global varible Definition
global X Y  XC YC N M NI NJ NIM NJM NIJ NM
global NL NR NB NT BCl BCr BCb BCt FI
global rho u v Xe Xw Yn Ys Gamma DX DY
global AW AE AN AS AP Q IC 

%=====Inputs
N=10; % No. CVs in X direcition
M=10; % No. CVs in Y direcition
Xmin=0; %Minimum X
Xmax=2; %Maximum X
Ymin=0; %Minimum Y
Ymax=1; %Maximum Y
EXX=1; %Expansion Factor ---X Direction
EXY=1; %Expansion Factor ---X Direction
rho=1; %Density
Gamma=0; %Diffussion Coefficient
IC=1; %Convective Flux Discritization Scheme:
%0:Centeral Difference 1:First Order Upwind  
%2:Hybrid              3:Power Law
%-----------------------------------------------------
iMethod=0; % Solution Method:
%0 MATLAB Direct Solver 
%2:Jacobi(Matrix) 3:Jacobi 4:Gauss-Seidel(Matrix) 5:Gauss-Seidel
%6: Symmetric Gauss-Seidel Matrix 7:Symmetric Gauss-Seidel
%8:SOR(Matrix) 9:SOR 10:Symmetric SOR Matrix
%11:Red-Black Gauss-Seidel  
%(13:Conjugate Gradient   14: Preconditioned Conjugate Gradient)---> Use
%this solve with causion for ,for Convection-Diffission Matrix is not
%Symmetric,Specially for high Peclet Number
%15: BiConjugate gradient
espSolver=1e-6; % Maximum Residual of Matrix Iterative Solver
MaxITSolver=100000; % Maximum Iteraton of Matrix Iterative Solver
%-----------------------------------------------------
%-----Boundary Conditions
% Set & Modify Boundary Conditions in SetBC file
%=====================================================
%Variable Allocation
NI=N+2;  % Two Extra points for first & last Cvs(Boundary)
NIM=NI-1; % No CVS +1, the NIM CV is on the boundary
NJ=M+2;  % Two Extra points for first & last Cvs(Boundary)
NJM=NJ-1; % No CVS +1, the NIM CV is on the boundary
NIJ=NI*NJ;  %Number of cells including the boundary cells
NM=N*M; %Number of internal Cells
%----------------------
X=zeros(1,NI);  %Grid Nodes
XC=zeros(1,NI); %Cell Centers
Y=zeros(1,NJ);  %Grid Nodes
YC=zeros(1,NJ); %Cell Centers
FI=zeros(NI,NJ); %Solution Variavle
Q=zeros(NM,1);  % Source Term
AE=zeros(NM,1);  %East Coe.
AW=zeros(NM,1);  %West Coe.
AS=zeros(NM,1);  %East Coe.
AN=zeros(NM,1);  %West Coe.
AP=zeros(NM,1);  %Point Coe.
Q=zeros(NM,1);  % Source Term
Xw=zeros(1,NI);
Xe=zeros(1,NI);
Ys=zeros(1,NJ);
Yn=zeros(1,NJ);
u=zeros(NIM,NJM); %Known Velociy Field @ the x direction
v=zeros(NIM,NIM); %Known Velociy Field @ the y direction
T=zeros(NM,1); %Solution in the interior cells
B=zeros(NM,5); %Auxillary Matrix for Constructing Sparse Marix
%==========================================================================
%Grid Generation
[X,XC,Xw,Xe,DX,Y,YC,Ys,Yn,DY] = Grid2d(Xmin,Xmax,Ymin,Ymax,N,M,EXX,EXY);
% V=sqrt(u^2+v^2);
% P=rho*V*(Xmax-Xmin)*(Ymax-Ymin)/Gamma; %Peclet Number
% fprintf('\nThe Peclet Number is %2.2f\n',P)
%Calaculate the know velocity filed u & v
Velocity
%==============Set BC and Fill into Fi
%-----------Modify SetBC file for adjasting boundary conditions
SetBC;
T(:)=0;   %Field Initialazation


%% ========================================================================
tic;
%Calculate Coefficients & Assemble Solution Matrix AT=b
%-----Interior Cells
for I=2:NIM
    for J=2:NJM       
        Discritize(I,J)
    end
end

%==========================================================================
%----Apply BCs by invoking into source Terms
BCs
%==========================================================================
%Assemble Matrix & Solve
AE=AE(1:end-1);
AW=AW(2:end);
AN=AN(1:end-N);
AS=AS(N+1:end);
%Sparse Matrix Construction
B(1:NM-1,2)=-AW;
B(:,3)=AP;
B(2:NM,4)=-AE;
B(1:NM-N,1)=-AN;
B(N+1:NM,5)=-AS;

A=spdiags(B,[-N -1 0 1 N],NM,NM);
% A=diag(AP)-diag(AE,1)-diag(AW,-1)-diag(AN,N)-diag(AS,-N);
%================================================================
%Solve Algebraic system with any Direct or Iterative Solver
[T] = Solve(A,T,Q,MaxITSolver,espSolver,iMethod); %Solution in the Interior Cells
toc
T=reshape(T,N,M); %Convert to 2D Array

%% ==========================================================================
T=reshape(T,N,M); %Convert to 2D Array
%Fill into solution Varibles and Interpolate IF Neaded
FI(2:NIM,2:NJM)=T;
% Interpolate FI @ Boundaries with specified q
%left & Right Boundary interpolation(IF Flux BCs specified)
for j=2:NJM;
    %Left Boundary
    if NL==1, FI(1,j)=BCl(j)*(XC(1)-XC(2))+FI(2,j);end
    %Right Boundary
    if NR==1, FI(NI,j)=BCr(j)*(XC(NI)-XC(NIM))+FI(NIM,j);end
end
%Top & Bottom Boundary interpolation(IF Flux BCs specified)
for i=2:NIM;
    %Top Boundary
    if NT==1, FI(i,1)=BCt(i)*(YC(1)-YC(2))+FI(i,2);end
    %Bottom Boundary
    if NB==1, FI(i,NJ)=BCb(i)*(YC(NJ)-XC(NJM))+FI(i,NJM);end
end

%--------Top & Bottom Boundary
% %Plot Result
[xc yc]=ndgrid(XC,YC);  %2D mesh grid for contour
figure
contourf(xc,yc,FI);
xlabel('X')
ylabel('Y')
title('Contour of \phi')
colorbar
figure
surf(xc,yc,FI);
xlabel('X')
ylabel('Y')
title('Surface of \phi')

% %---Compute Fluxes @ Boundary Points
% %Left & Right Boundaries
% for j=1:M;
%     q_l(j)=(FI(1,j)-FI(2,j))/abs(X(1)-X(2));   %Flux @ the Left Boundary
%     q_r(j)=(FI(N,j)-FI(N-1,j))/abs(X(N-1)-X(N));   %Flux @ the Left Boundary
% end
% %Top & Bottom Boundaries
% for i=1:N;
%     q_t(i)=(FI(i,M)-FI(i,M-1))/abs(Y(M)-Y(M-1));   %Flux @ the Bottom Boundary
%     q_b(i)=(FI(i,1)-FI(i,2))/abs(Y(1)-Y(2));   %Flux @ the Top Boundary
% end
% %-----Compute and Plot Fluxes @ East & South CV Faces
% for j=1:M;
%     for i=1:N;
%         qx(i,j)=(FI(i,j)-FI(i+1,j))/abs(X(i+1)-X(i));
%     end
%     qx(N,j)=-q_r(j);
% end
% for i=1:N;
%     for j=1:M;
%         qy(i,j)=(FI(i,j)-FI(i,j+1))/abs(Y(j+1)-Y(j));
%     end
%     qy(i,M)=-q_b(i);   
% end
% %Report Flux in inbalence @ Boundaries
% fprintf('\nFlux Imbalance @ Right & Left Boundaries is %2.2e\n',sum(q_r)+sum(q_l));
% fprintf('Flux Imbalance @ Top & Bottom Boundaries is %2.2e\n',sum(q_t)+sum(q_b));
% 
% %Plot Flux vectors
% figure('Name','Flux Vectors');
% contour(xc',yc',FI)
% hold on
% quiver(xc(2:NJM,2:NJM)',yc(2:NIM,2:NJM)',qx,qy)
% hold off
% % 
disp('Good Lock, Mohammad Aghakhani')






