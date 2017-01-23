function [X,Xvol,Xw,Xe] = Grid1d(L,m)
% Generates A 1D Grid m is Number of Cells, L the lenght of the Domain 
%Returns X=Vertices Coordinate, Xvol=Control Volume Cooor. XW,XE are West
%and East X of each point
% Vetices are Placed on Boundaries & Finite Volume Cells Are in between
% primary Veticis
%2010, Mohammad Aghakhani

%Generate Vertices
X(1)=0;
Dx=L/(m-1);
for i=2:m;
    X(i)=X(i-1)+Dx;
end
%Locations of the Control Volumes
i=2:m;
Xvol(i-1)=0.5*(X(i-1)+X(i));
% Calculate Xe and Xw
i=2:m;
Xw(1)=0;
Xw(i)=abs(X(i-1)-X(i));
i=1:m-1;
Xe(i)=abs(X(i+1)-X(i));
Xe(m)=0;


end

