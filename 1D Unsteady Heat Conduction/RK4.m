function [ U ] = RK4(Uold,Dt,i)
%Explit time Iteration by 4th Order Rung-Kutta Iteration
% Uold is the value @ old itertioan & Dt is time step

k1=RHS(Uold,i)*Dt;
k2=RHS(Uold+(RHS(Uold,i)*Dt)/2,i)*Dt;
k3=RHS(Uold+k2/2,i)*Dt;
k4=RHS(Uold+k3,i)*Dt;
U=Uold(i)+k1/6+k2/3+k3/3+k4/6;
end

