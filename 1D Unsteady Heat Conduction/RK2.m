function [ U ] = RK2(Uold,Dt,i)
%Explit time Iteration by 2nd Order Rung-Kutta Iteration
% Uold is the value @ old itertioan & Dt is time step

%k1=RHS(Uold,i)*Dt;
%k2=RHS(Uold+k1/2)*Dt
%k2=RHS(Uold+(RHS(Uold,i)*Dt)/2,i)*Dt;
U=Uold(i)+RHS(Uold+(RHS(Uold,i)*Dt)/2,i)*Dt;
end

