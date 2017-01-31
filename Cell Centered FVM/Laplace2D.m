%This Program Solves 2D Steady Laplace Eq. by node Centerd FVM
clear
clc
close all

%-----------------Inputs
N=10; % No. CVs in X direcition
M=10; % No. CVs in Y direcition
XMIN=0; %Minimum X
XMAX=1; %Maximum X
EXX=1; %Expansion Factor ------X Direction
YMIN=0; %Minimum X
YMAX=1; %Maximum X
EXY=1; %Expansion Factor ------X Direction

%====Boundry Conditions Input Setting
NL=0; % 0:Fixed Temp. 1:Flux Bc
NR=0; % 0:Fixed Temp. 1:Flux Bc
NB=1; % 0:Fixed Temp. 1:Flux Bc
NT=1; % 0:Fixed Temp. 1:Flux Bc
FIl=100;   % FI at The Left Boundary, Applies only if NL=0
FIr=0;   % FI at The Right Boundary, Applies only if NR=0
FIb=0; % FI at The Bottom Boundary, Applies only if NB=0
FIt=1; % Temperature at The Top Boundary, Applies only if NT=0
Ql=0; % Flux at The Left Boundary, Applies only if NL=1
Qr=0; % Flux at The Right Boundary, Applies only if NR=1
Qb=0; % Flux at The Bottom Boundary, Applies only if NB=1
Qt=0; % Flux at The Top Boundary, Applies only if NT=1
%=======================================================================

NI=N+2;  % Two Extra points for first & last Cvs(Boundary)
NIM=NI-1; % No CVS +1, the NIM CV is on the boundary
NJ=M+2;  % Two Extra points for first & last Cvs(Boundary)
NJM=NJ-1; % No CVS +1, the NIM CV is on the boundary
NIJ=NI*NJ;
NM=N*M; %Number of internal Cells
%----------------------
X=zeros(1,NI);  %Grid Nodes
XC=zeros(1,NI); %Cell Centers
Y=zeros(1,NJ);  %Grid Nodes
YC=zeros(1,NJ); %Cell Centers
FI=zeros(NI,NJ); %Solution Variable
Q=zeros(NM,1);  % Source Term
AE=zeros(NM,1);  %East Coe.
AW=zeros(NM,1);  %West Coe.
AN=zeros(NM,1);  %East Coe.
AS=zeros(NM,1);  %West Coe.
AP=zeros(NM,1);  %Point Coe.
B=zeros(NM,5); %Auxillary Matrix for Constructing Sparse Marix
%==========================================================
%----------Grid Generation
%DEFINE GRID IN X-DIRECTION
if EXX==1
    DX=(XMAX-XMIN)/N;
else
    DX=(XMAX-XMIN)*(1-EXX)/(1-EXX^N);
end

X(1)=XMIN;
for I=2:NIM
    X(I)=X(I-1)+DX;
    DX=DX*EXX;
end
X(NI)=X(NIM);

XC(1)=X(1); %first and last CV center lay on the boundaries
XC(NI)=X(NIM);
for I=2:NIM
    XC(I)=0.5*(X(I)+X(I-1));
end

%DEFINE GRID IN Y-DIRECTION
if EXY==1
    DY=(YMAX-YMIN)/M;
else
    DY=(YMAX-YMIN)*(1-EXY)/(1-EXY^M);
end

Y(1)=YMIN;
for I=2:NJM
    Y(I)=Y(I-1)+DY;
    DY=DY*EXY;
end
Y(NJ)=Y(NJM);

YC(1)=Y(1); %first and last CV center lay on the boundaries
YC(NJ)=Y(NJM);
for I=2:NJM
    YC(I)=0.5*(Y(I)+Y(I-1));
end
%=====================================================================
%Fill Boundary Values
j=1:NJ; %--------Left & Right
if NL==0,FI(1,j)=FIl; end    %Left BC
if NR==0,
%       FI(NI,j)=X(NI)*Y(j); For Validation
    FI(NI,j)=FIr; 
end %Right BC
i=1:NI;  % Top & Bottom
if NB==0,FI(i,1)=FIb; end  %Bottom Boundary
if NT==0,
%       FI(i,NJ)=X(i)*Y(NJ); For Validation
    FI(i,NJ)=FIt;
end %Top Boundary
%=====================================================================
%CALCULATE ELEMENTS OF MATRIX [A]
for I=2:NIM
    for J=2:NJM
        %Auxillary for internal cells
        i=I-1;
        j=J-1;
        IJ=(j-1)*N+i; %Convert 2D to 1D index
        AE(IJ)=(Y(J)-Y(J-1))/(XC(I+1)-XC(I)); %Dy/Xw
        AW(IJ)=(Y(J)-Y(J-1))/(XC(I)-XC(I-1)); %Dy/Xe
        AN(IJ)=(X(I)-X(I-1))/(YC(J+1)-YC(J)); %Dy/Xw
        AS(IJ)=(X(I)-X(I-1))/(YC(J)-YC(J-1)); %Dy/Xe
        AP(IJ)=(AE(IJ)+AW(IJ)+AN(IJ)+AS(IJ));
        VOL=(X(I)-X(I-1))*(Y(J)-Y(J-1));
        Q(IJ)=0;%sin(2*pi*XC(I)/XMAX)*sin(2*pi*YC(J)/YMAX)*VOL
    end
end
%======================================================================
%IF source Term is not zero plot it
if any(Q) %If Source Term is non-zero plot it
    figure
    [xc yc]=meshgrid(XC(2:NIM),YC(2:NJM));  %2D mesh grid 
    surf(xc',yc',reshape(Q,N,M));
    xlabel('X')
    ylabel('Y')
    title('Surface of Source Term')
end
%======================================================================
%Apply BCs
%--------Left & Right
for j=1:M;
    i=1; %Left Boundary
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch NL %Dirichlet BC
        case 0
            Q(IJ)=Q(IJ)+AW(IJ)*FI(1,j);  %BCs at first CV
            AW(IJ)=0;
        case 1    %Neuman BC on the Left
            Q(IJ)=Q(IJ)+Ql;  %BCs at first CV
            AP(IJ)=AP(IJ)-AW(IJ);
            AW(IJ)=0;
    end
    i=N;
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch NR
        case 0
            Q(IJ)=Q(IJ)+AE(IJ)*FI(NI,j);  %BCs at first CV
            AE(IJ)=0;
        case 1    %Neuman BC on the Right
            Q(IJ)=Q(IJ)+Qr;  %BCs at first CV
            AP(IJ)=AP(IJ)-AE(IJ);
            AE(IJ)=0;
    end
end
%--------Bottom & Top
for i=1:N;
    %Bottom Boundary
    j=1;
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch NB  %Dirichlet BC
        case 0
            Q(IJ)=Q(IJ)+AS(IJ)*FI(i,1);  %BCs at first CV
            AS(IJ)=0;
        case 1 %Neuman BC on the Bottom
            Q(IJ)=Q(IJ)+Qb;  %BCs at first CV
            AP(IJ)=AP(IJ)-AS(IJ);
            AS(IJ)=0;
    end
    %Top Boundary
    j=M;
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch NT
        case 0 %Dirichlet BC
            Q(IJ)=Q(IJ)+AN(IJ)*FI(i,NJ);  %BCs at first CV
            AN(IJ)=0;
        case 1 %Neuman BC on the Top
            Q(IJ)=Q(IJ)+Qt;  %BCs at first CV
            AP(IJ)=AP(IJ)-AN(IJ);
            AN(IJ)=0;
    end
end
%================================================================
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
T=A\Q; %Solution in the Interior Cells
%==============================================================
%Post-Proceccing
T=reshape(T,N,M); %Convert to 2D Array
%Fill into solution Varibles and Interpolate IF Neaded
FI(2:NIM,2:NJM)=T;
% Interpolate FI @ Boundaries with specified q
%left & Right Boundary interpolation(IF Flux BCs specified)
j=2:NJM;
%Left Boundary
if NL==1, FI(1,j)=Ql*(XC(1)-XC(2))+FI(2,j);end
%Right Boundary
if NR==1, FI(NI,j)=Qr*(XC(NI)-XC(NIM))+FI(NIM,j);end
%Top & Bottom Boundary interpolation(IF Flux BCs specified)
i=2:NIM;
%Top Boundary
if NT==1, FI(i,1)=Qt*(YC(1)-YC(2))+FI(i,2);end
%Bottom Boundary
if NB==1, FI(i,NJ)=Qb*(YC(NJ)-XC(NJM))+FI(i,NJM);end

%--------Top & Bottom Boundary
% %Plot Result
[xc yc]=meshgrid(XC,YC);  %2D mesh grid for contour
figure
contourf(xc',yc',FI);
xlabel('X')
ylabel('Y')
title('Contour of \phi')
colorbar
figure
surf(xc',yc',FI);
xlabel('X')
ylabel('Y')
title('Surface of \phi')
%==========================================Validation
%Exact Solution: exact solution
% Fi = x * y when Dirichlet boundary condition is applied at all boundaries 
% In that case, FI=0 at south and west boundary 
%and at east and north boundaries  must be prescribed:
% for i=1:NI;
%     for j=1:NJ;
%         FIExact(i,j)=XC(i)*YC(j);
%     end
% end
% figure
% contourf(xc',yc',FIExact);
% xlabel('X')
% ylabel('Y')
% title('Contour of Exact \phi')
% colorbar
% error=FIExact-FI;
% figure
% surf(error)
% xlabel('X')
% ylabel('Y')
% title('Sourface of Error')
% fprintf('\nThe Norm of Error is %2.2e\n',norm(error))
%================================================================
%---Compute Fluxes @ Boundary Points
%Left & Right Boundaries
for j=1:M;
    q_l(j)=(FI(1,j)-FI(2,j))/abs(X(1)-X(2));   %Flux @ the Left Boundary
    q_r(j)=(FI(N,j)-FI(N-1,j))/abs(X(N-1)-X(N));   %Flux @ the Left Boundary
end
%Top & Bottom Boundaries
for i=1:N;
    q_t(i)=(FI(i,M)-FI(i,M-1))/abs(Y(M)-Y(M-1));   %Flux @ the Bottom Boundary
    q_b(i)=(FI(i,1)-FI(i,2))/abs(Y(1)-Y(2));   %Flux @ the Top Boundary
end
%-----Compute and Plot Fluxes @ East & South CV Faces
for j=1:M;
    for i=1:N;
        qx(i,j)=(FI(i,j)-FI(i+1,j))/abs(X(i+1)-X(i));
    end
    qx(N,j)=-q_r(j);
end
for i=1:N;
    for j=1:M;
        qy(i,j)=(FI(i,j)-FI(i,j+1))/abs(Y(j+1)-Y(j));
    end
    qy(i,M)=-q_b(i);   
end
%Report Flux in inbalence @ Boundaries
fprintf('\nFlux Imbalance @ Right & Left Boundaries is %2.2e\n',sum(q_r)+sum(q_l));
fprintf('Flux Imbalance @ Top & Bottom Boundaries is %2.2e\n',sum(q_t)+sum(q_b));

%Plot Flux vectors
figure('Name','Flux Vectors');
contour(xc',yc',FI)
hold on
quiver(xc(2:NJM,2:NJM)',yc(2:NIM,2:NJM)',qx,qy)
hold off
% 
disp('Good Lock, Mohammad Aghakhani')

