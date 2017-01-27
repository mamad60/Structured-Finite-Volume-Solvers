function [rhs] = RHS(T,i)
%Returs Right Hand Side of EQ. for time iteration
%------Global Varibles
global ke kw rho Cp Xe Xw Dx Sp Sc
%--------------
rhs=((ke(i)*(T(i+1)-T(i))/Xe(i))-(kw*(T(i)-T(i-1))/Xw(i))+(Sp+Sc)*Dx(i))/(rho*Cp*Dx(i));

end

