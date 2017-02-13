%% This Scripts Solves 2D Steady Diffussion Equation
% d/dx(kdT/dx)+S=0, From Patankar Book, Section 4-2
%S=Sp*FI+Su is the source term per unit volume
% Both Thermal Conductivity and Source therms can be varied with time & x
%                                    N
%Discritization Stencil   |--->  W---P---E
%                                    S
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


global X Y gStat XC YC Xmax Xmin EXX sStat
global DX DY Ymin Ymax EXY
global BCw BCe BCn BCs
global Gammae Gammaw Gammas Gamman
global Su Sp aE aW aN aS FI F
global N M NI NIM NJ NJM NM

%--Inputs
N=20; % No. CVs in X direcition
M=20; % No. CVs in Y direcition
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
FI=zeros(NI,NJ); %Solution Variable
Su=zeros(NM,1);  % Source Term---Su
Sp=zeros(NM,1);  % Source Term---Sp
aE=zeros(NM,1);  %East Coe.
aW=zeros(NM,1);  %West Coe.
aS=zeros(NM,1);  %East Coe.
aN=zeros(NM,1);  %West Coe.
aP=zeros(NM,1);  %Point Coe.
Xw=zeros(1,NI);
Xe=zeros(1,NI);
Ys=zeros(1,NJ);
Yn=zeros(1,NJ);
Gammae=zeros(NI,NJ);
Gammaw=zeros(NI,NJ);
Gammas=zeros(NI,NJ);
Gamman=zeros(NI,NJ);
FI=zeros(NI,NJ);
BCw=zeros(NJ,1);
BCe=zeros(NJ,1);
BCn=zeros(1,NI);
BCs=zeros(1,NI);
F=zeros(NM,1); %Temporary Solution in the interior cells
B=zeros(NM,5); %Auxillary Matrix for Constructing Sparse Marix
%==========================================================================
%Problem Setup
Setup1Dcx
%Grid Generation
[X,XC,Xw,Xe,DX,Y,YC,Ys,Yn,DY] = Grid2d(Xmin,Xmax,Ymin,Ymax,N,M,EXX,EXY);
%Compute Diffussion Coe
BiHarmonic;

%Compute Source Terms
[Su Sp]=Source(sStat);
%%==========================================================================
%Construct Coe. for the Internal cells
for i=2:NIM
    for j=2:NJM
        I=i-1;
        J=j-1;
        IJ=(J-1)*N+I; %Convert 2D to 1D index
        aW(IJ)=Gammaw(i,j)/Xw(i)*DY(j);
        aE(IJ)=Gammae(i,j)/Xe(i)*DY(j);
        aN(IJ)=Gamman(i,j)/Yn(j)*DX(i);
        aS(IJ)=Gammas(i,j)/Ys(j)*DX(i);
    end
end
%--------Apply Boundary Conditions
AdjustBC;
%------Construct aP
aP=aW+aE+aS+aN-Sp;

%==========================================================================
%Assemble Matrix & Solve
aE=aE(1:end-1);
aW=aW(2:end);
aN=aN(1:end-N);
aS=aS(N+1:end);
%Sparse Matrix Construction
B(1:NM-1,2)=-aW;
B(:,3)=aP;
B(2:NM,4)=-aE;
B(1:NM-N,1)=-aN;
B(N+1:NM,5)=-aS;
A=spdiags(B,[-N -1 0 1 N],NM,NM);
%+++++++++Solve Equations
F=A\Su;
%%=========================================================================
%Post-Process
%-------------Adjust FI for Post-Processing
AdjustPost
%--------------------------------------------------------------------------
%Plot Result
[xc yc]=ndgrid(XC,YC);  %2D mesh grid for contour1
scrsz = get(0,'ScreenSize');
%----Contour Plot
figure('Name','Contour Plot','NumberTitle','off','OuterPosition',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
contourf(xc,yc,FI);
xlabel('X')
ylabel('Y')
title('Contour of \phi')
colorbar
%----Surface Plot
figure('Name','Surface Plot','NumberTitle','off','OuterPosition',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
surf(xc,yc,FI);
xlabel('X')
ylabel('Y')
title('Surface of \phi')
grid on;
xlim([Xmin Xmax])
ylim([Ymin Ymax])
shading interp;
%Plot profile of FI @  mid-section of the domain
figure('Name','Profile Plot(Y)','NumberTitle','off','OuterPosition',[scrsz(3)/2  1 scrsz(3)/2 scrsz(4)/2]);
%Y Direction
plot(XC,FI(:,floor(M/8)),XC,FI(:,floor(M/4)),XC,FI(:,floor(M/2)),XC,FI(:,floor(3*M/4)),XC,FI(:,floor(7*M/8)));  
legend('M/8','M/4','M/2','3*M/4','7*M/8','Location','SouthEast')
xlabel('X')
ylabel('\phi')
xlim([Xmin Xmax])
title('Profile of \phi on Horizontal Lines')
figure('Name','Profile Plot(X)','NumberTitle','off','OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)/2]);
%X Direction
plot(YC,FI(floor(N/8),:),YC,FI(floor(N/4),:),YC,FI(floor(N/2),:),YC,FI(floor(3*N/4),:),YC,FI(floor(7*N/8),:));  
legend('N/8','N/4','N/2','3*N/4','7*N/8','Location','NorthEast')
xlabel('Y')
ylabel('\phi')
title('Profile of \phi on Vertical Lines')
xlim([Ymin Ymax])
%Plot Difiussion Coe.
if gStat==2
    figure
    contourf(xc(2:NIM,2:NJM),yc(2:NIM,2:NJM),Gammaw(2:NIM,2:NJM))
    xlabel('X')
    ylabe('Y')
    title('Difiussion Coefficient Profie(\gamma_w)')
end
%--------------------------------------------------------------------------
%-------Compute Boundary Fluxes
for j=1:NJ;
    Qw(j)=Gammaw(2,j)*(FI(1,j)-FI(2,j))/abs((XC(1)-XC(2)))*DY(j);
    Qe(j)=Gammae(NIM,j)*(FI(NI,j)-FI(NIM,j))/abs((XC(NI)-XC(NIM)))*DY(j);
end
for i=1:NI;
    Qs(i)=Gammas(i,2)*(FI(i,1)-FI(i,2))/abs((YC(1)-YC(2)))*DX(i);
    Qn(i)=Gamman(i,NJM)*(FI(i,NJ)-FI(i,NJM))/abs((XC(NJ)-XC(NJM)))*DX(i);
end
%------Compute & Plot  Fluxes
for j=1:NJ;
    for i=1:NIM
        Qx(i,j)=Gammaw(i,j)*(FI(i,j)-FI(i+1,j))/abs((XC(i)-XC(i+1)))*DY(j);
    end
    Qx(NI,j)=Qw(j);
    Qx(1,j)=-Qe(j);
end
for i=1:NI;
    for j=1:NJM
        Qy(i,j)=Gammas(i,j)*(FI(i,j)-FI(i,j+1))/abs((YC(j)-YC(j+1)))*DX(i);
    end
    Qy(i,NJ)=Qs(i);
    Qy(i,1)=-Qn(i);
end

figure('Name','Profile Plot','NumberTitle','off','OuterPosition',[1  1 scrsz(3)/2 scrsz(4)/2]);
quiver(xc(2:NIM,2:NIM),yc(2:NIM,2:NIM),Qx(2:NIM,2:NIM),Qy(2:NIM,2:NIM))
xlabel('X')
ylabel('Y')
title('Flux Vector Plot of \phi')
axis([Xmin Xmax Ymin Ymax])
grid on;


%------% Report Fluxes
fprintf('\nTotal Flux of FI @ the West Boundary is:\t%2.4f',sum(Qw));
fprintf('\nTotal Flux of FI @ the East Boundary is:\t%2.4f',sum(Qe));
fprintf('\nTotal Flux of FI @ the South Boundary is:\t%2.4f',sum(Qs));
fprintf('\nTotal Flux of FI @ the North Boundary is:\t%2.4f',sum(Qn));
fprintf('\nFlux Inbalance @ the Boundaries is:\t%2.4f\n',sum(Qw)+sum(Qe)+sum(Qs)+sum(Qn));

disp('Good Lock, Mohammad Aghakhani')








