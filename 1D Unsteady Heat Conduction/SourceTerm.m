function [Sc,Sp] = SourceTerm(status,x,T)
%Returns Source Term @ x position & T temperature
% Status=0: No Source Term 1:Constant 2:Function of X 3:Function of
% Temperature 4: Function of Both
%S=Sc+Sp*T  |------> Linerized Source term
%For Boundness od the solution Sp must be Negative
switch status
    case 0
        %No Source Term
        Sc=0;
        Sp=0;
    case 1
        %Enter Fixed Source Term here
        Sc=1;%1e-3;
        Sp=0;
        
    case 2
        %If SourceTerm varies with x enter expresssion here
        % Sc,Sp=f(x)
%         if 0.2<=x<=0.6
%             Sc=1;
%         else
%             Sc=0;
%         end
        Sc=abs(x-0.5)*1;%1e-3+x*1e-3;
        Sp=0;
    case 3
        %If SourceTerm varies both with T enter expresssion here
        %Sc,Sp=f(T)
        Sp=-2*1e-3;
        Sc=(3+9*T)*1e-5;
    case 4
        %If Source Term varies both with x & T enter expresssion here
        %Sc, Sp=f(x,T)
        Sc=1e-4+x*1e-3;
        Sp=-3*1e-5;
        
end

end

