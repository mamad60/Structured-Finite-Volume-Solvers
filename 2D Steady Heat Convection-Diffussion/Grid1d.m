function [X,XC,Xw,Xe,DX] = Grid1d(Xmin,Xmax,N,EXX)
% Generates A 1D Grid N is Number of Cells, for Cell-Centered FVM
%Returns X=Grid Coordinate, XC=Center of Control Volume Cooor. EX=Grid Expasion Coe.
%XW,XE are West and East X of each point, Dx is length of Eac CV
%First & Last XC lays on the boundaries. 
%2010, Mohammad Aghakhani

NI=N+2;  % Two Extra points for first & last Cvs(Boundary)
NIM=NI-1; % No CVS +1, the NIM CV is on the boundary
%-----------------------------------------
X=zeros(1,NI);
XC=zeros(1,NI);
Xw=zeros(1,NI);
Xe=zeros(1,NI);
DX=zeros(1,NIM);
%----------Grid Generation
%DEFINE GRID IN X-DIRECTION
if EXX==1
    DX(2)=(Xmax-Xmin)/N;
else
    DX(2)=(Xmax-Xmin)*(1-EXX)/(1-EXX^N);
end
DX(1)=DX(2);
X(1)=Xmin;
for I=2:NIM
    X(I)=X(I-1)+DX(I);
    DX(I+1)=DX(I)*EXX;
end
X(NI)=X(NIM);
DX(1)=0;
DX(NI)=0;
XC(1)=X(1); %first and last CV center lay on the boundaries
XC(NI)=X(NIM);
for I=2:NIM
    XC(I)=0.5*(X(I)+X(I-1));
end

% Calculate Xe and Xw
i=2:NI;
Xw(1)=0;
Xw(i)=abs(XC(i)-XC(i-1));
i=1:NIM;
Xe(i)=abs(XC(i+1)-XC(i));
Xe(N+2)=0;


end

