function [ U ] = RK4(Uold,Dt,i,j)
%Explit time Iteration by 4th Order Rung-Kutta Iteration
% Uold is the value @ old itertioan & Dt is time step

k1=RHS(Uold,i,j)*Dt;
k2=RHS(Uold+k1/2,i,j)*Dt;
k3=RHS(Uold+k2/2,i,j)*Dt;
k4=RHS(Uold+k3,i,j)*Dt;
U=Uold(i,j)+k1/6+k2/3+k3/3+k4/6;
end

