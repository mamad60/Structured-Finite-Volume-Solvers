function [rhs] = RHS(T,i,j)
%Returs Right Hand Side of EQ. for time iteration
%------Global Varibles
global ke kw rho Cp Xe Xw Dx Dy kn ks Yn Ys  Sc Sp
%---------------------------------
rhs=( ((ke(i,j)*(T(i+1,j)-T(i,j))/Xe(i,j))-(kw(i,j)*(T(i,j)-T(i-1,j))/Xw(i,j)))*Dy(i,j)+...
    ((ks(i,j)*(T(i,j+1)-T(i,j))/Ys(i,j))-(kn(i,j)*(T(i,j)-T(i,j-1))/Yn(i,j)))*Dx(i,j)+...+
    (Sc(i,j)+Sp(i,j)*T(i,j))*Dx(i,j)*Dy(i,j) )/(rho*Cp*Dx(i,j)*Dy(i,j));

end

