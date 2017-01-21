function [ k ] = ThermalConductivity(status,x,T)
%Returns Thermal Condunctivity of the Material @ x position & T temperature
% Status=0:Constant 1:Function of X 2:Function of Temperature 2: Function
% of both Temerature and Time
switch status
    case 0
        %Enter Fixed Thermal Condunctivity here
        k=1e-3;
    case 1
        %If Thermal Condunctivity varies with x enter expresssion here
        % k=f(x)
        k=1e-3+x*1e-3;
    case 2
        %If Thermal Condunctivity varies both with T enter expresssion here
        %k=f(T)
        k=1e-3+1e+4*abs(T-273);
    case 3
        %If Thermal Condunctivity varies both with x & T enter expresssion here
        %k=f(x,T)
        k=1e-3+(x*1e-3)*+1e-4+1e-4*abs(T-273);
end
      
end

