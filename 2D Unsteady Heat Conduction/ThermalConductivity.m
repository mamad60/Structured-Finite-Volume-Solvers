function [ k ] = ThermalConductivity(status,x,y,T)
%Returns Thermal Condunctivity of the Material @ x position & T temperature
% Status=0:Constant 1:Function of X,y 2:Function of Temperature 2: Function
% of Temerature and Space-----%2d Version
switch status
    case 0
        %Enter Fixed Thermal Condunctivity here
        k=1e-3;
    case 1
        %If Thermal Condunctivity varies with x(and/Or)y enter expresssion here
        % k=f(x,y)
        k=1e-3+(x+y)*1e-3;
    case 2
        %If Thermal Condunctivity varies both with T enter expresssion here
        %k=f(T)
        k=1e-3+1e-2*abs(T-290);
    case 3
        %If Thermal Condunctivity varies  with x,y & T enter expresssion here
        %k=f(x,y,T)
        k=1e-3+(x+y)*1e-3+1e-2*abs(T-290);
end
      
end

