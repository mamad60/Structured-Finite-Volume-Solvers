function [X,XC,Xw,Xe,DX,Y,YC,Ys,Yn,DY] = Grid2d(Xmin,Xmax,Ymin,Ymax,N,M,EXX,EXY)
% Generates A 2D Grid N is Number of Cells, for Cell-Centered FVM
%Returns X=Grid Coordinate, XC=Center of Control Volume Cooor. EX=Grid Expasion Coe.
%XW,XE are West and East X of each point, Dx is length of Eac CV
%First & Last XC lays on the boundaries. 
%For Y direction definitions are the same, M is No. CVs in the y direction
%2010, Mohammad Aghakhani

%Define Grid by calling Grid1d two time for each direction
%----------Grid Generation

%DEFINE GRID IN X-DIRECTION
[X,XC,Xw,Xe,DX] = Grid1d(Xmin,Xmax,N,EXX);
%DEFINE GRID IN Y-DIRECTION
[Y,YC,Ys,Yn,DY] = Grid1d(Ymin,Ymax,M,EXY);

end

