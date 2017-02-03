function Velocity
%Sets Velocity known @ Each at each point of the grid
% u and v are NIM(N+1) and NJM(M+1)
%U  at cell face 'e'; V  at cell face 'n'
global X Y  XC YC N M NI NJ NIM NJM NIJ NM u v
%--------------------------------------
u=zeros(NIM,NJM); %Known Velociy Field @ the x direction
v=zeros(NIM,NIM); %Known Velociy Field @ the y direction

u0=1;
v0=10;
%Enter the Expression for u and v velocity here
for i=1:NIM
    for j=1:NJM
        u(i,j)=u0*cos(45*pi/180);
        v(i,j)=u0*sin(45*pi/180);
    end
end


end

