function [ ]= Grid2d(Lx,Ly,m,n)
% Generates A 2D m x n Grid m,n is Number of Cells,
%Lx the lenght of the Domain y direction and Ly the lenght of the Domain in y direction
%Returns X,Y=Vertices Coordinate, Xvol,YVol=Control Volume Coor.
%XW and XE are West and East X of each point
%YN and YS are North and South Y of each point
%Vetices are Placed on Boundaries & Finite Volume Cells Are in between primary Veticis
%2010, Mohammad Aghakhani
%----Global variables
global X Y Xvol Yvol Xe Xw Yn Ys Dx Dy

%Generate Vertices
DX=Lx/(m-1);
DY=Ly/(n-1);
[X,Y]=ndgrid(0:DX:Lx,0:DY:Ly);

%Locations of the Control Volumes
j=1:n;
i=2:m;
    Xvol(i-1,j)=0.5*(X(i-1,j)+X(i,j));

i=1:m;
j=2:n;
    Yvol(i,j-1)=0.5*(Y(i,j-1)+Y(i,j));

%Calculate Xe and Xw
j=1:n;
i=2:m;
    Xw(1,j)=0;
    Xw(i,j)=abs(X(i,j)-X(i-1,j));

i=1:m-1;
    Xe(m,j)=0;
    Xe(i,j)=abs(X(i+1,j)-X(i,j));


%Calculate Xn and Xs
i=1:m;
j=2:n;
    Yn(i,1)=0;
    Yn(i,j)=abs(Y(i,j-1)-Y(i,j));


j=1:n-1;
    Ys(i,n)=0;
    Ys(i,j)=abs(Y(i,j+1)-Y(i,j));


%Calculate Dx of CVs
j=1:n;
Dx(1,j)=Xvol(1,j);
i=2:m-1;
    Dx(i,j)=abs(Xvol(i,j)-Xvol(i-1,j));

Dx(m,j)=X(m,j)-Xvol(m-1,j);
%Calculate Dy of CVs
i=1:m;
Dy(i,1)=Yvol(i,1);
j=2:n-1;
    Dy(i,j)=abs(Yvol(i,j)-Yvol(i,j-1));

Dy(i,n)=abs(Y(i,n)-Yvol(i,n-1));


