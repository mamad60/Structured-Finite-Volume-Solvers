%This Program Solves 1D Steady Laplace Eq. by node Centerd FVM
clear
clc
close all

%-----------------Inputs
N=5; % No. CVs in X direcition
XMIN=0; %Minimum X
XMAX=1; %Maximum X
EXX=1.1; %Expansion Factor ------X Direction
%-----------------------------------------------
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
B=zeros(N,3); %Auxillary Matrix for Constructing Sparse Marix

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

% FI(1)=0; %Left BC
FI(NI)=100; %Right BC
%.....CALCULATE ELEMENTS OF MATRIX [A] (For Laplace equation, Q(IJ)=0.;

for I=2:NIM
        AE(I)=1/(XC(I+1)-XC(I)); %1/Xw
        AW(I)=1/(XC(I)-XC(I-1)); %1/Xe
        AP(I)=(AE(I)+AW(I));   
        Q(I)=0;%sin(2*pi*XC(I));
end

if any(Q) %If Source Term is non-zero plot it
    figure
    plot(XC,Q)
    xlabel('X')
    ylabel('Source Term')
    title('Source Term')
end

%Neuman BC on the Left
q=-100;
Q(2)=Q(2)+q;  %BCs at first CV
AP(2)=AP(2)-AW(2);
% Q(2)=Q(2)+AW(2)*FI(1);  %BCs at first CV
% AW(2)=0;
Q(NIM)=Q(NIM)+AE(NIM)*FI(NI);  %BCs at Last CV by Invoking to Source Term
AE(NIM)=0;   %Break Link to the Right

%Assemble Matrix & Solve
AP=AP(2:end);
AE=AE(2:end-1);
AW=AW(3:end);
%Sparse Matrix Construction
B(1:N-1,1)=-AW;
B(:,2)=AP;
B(2:N,3)=-AE;
A=spdiags(B,[-1 0 1],N,N);
%  A=diag(AP)-diag(AE,1)-diag(AW,-1);
b=Q(2:NIM)';
T=A\b; %Solution in the Interior Cells
FI(2:NIM)=T(:);
%Interpolate FI @ Boundaries
%Left Boundary Interpolation, IF q is applied
% p=polyfit([XC(2) XC(3) XC(4)],[FI(2) FI(3) FI(4)],2);
% l=polyval(p,XC(1));
% FI(1)=l;
% %Right Boundary
% p=polyfit([XC(NIM) XC(NIM-1) XC(NIM-2)],[FI(NIM) FI(NIM-1) FI(NIM-2)]2)
% r=polyval(p,XC(NI))
% FI(NI)=r;
%Plot Result
figure
plot(XC,FI)
xlabel('X')
ylabel('\phi')
title('Solution Profile')
%Report Fluxes
ql=(T(1)-T(2))/(XC(3)-XC(2));
qr=(T(end)-T(end-1))/(XC(end-1)-XC(end-2));
fprintf('\nFlux of FI @ the Left Boundary is:\t%2.4f',ql);
fprintf('\nFlux of FI @ the Right Boundary is:\t%2.4f',qr);
fprintf('\nFlux Inbalance @ the Boundaries is:\t%2.4f\n',qr+ql);

disp('Good Lock, Mohammad Aghakhani')


