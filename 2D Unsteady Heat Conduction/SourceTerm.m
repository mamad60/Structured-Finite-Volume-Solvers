function [Sc,Sp] = SourceTerm(status,x,y,T)
%Returns Source Term @ x position & T temperature
% Status=0: No Source Term 1:Constant 2:Function of X(and/or) Y 3:Function of
% Temperature 4: Function of Both
%S=Sc+Sp*T  |------> Linerized Source term
%For Boundness of the solution Sp must be Negative
global m n

Sc=zeros(size(x));
Sp=zeros(size(x));
%Check iFf Source Term is Zero or Constant
switch status
    case 0
        %No Source Term
        Sc(:)=0;
        Sp(:)=0;
        return
    case 1
        %Enter Fixed Source Term here
        Sc(:)=1;%1e-3;
        Sp(:)=0;
        return
end
for i=1:m
    for j=1:n
        switch status
            case 2
                %If SourceTerm varies with x,y enter expresssion here
                %Sc,Sp=f(x,y)
                if (x(i,j)<=0.6 && x(i,j)>=1.2)
                  if (y(i,j)<=0.3 && y(i,j)>=0.6)
                    Sc(i,j)=1;
                  end
                else
                    Sc(i,j)=0;
                end
                %         Sc=abs(x-0.5)*(y-.05)*1;%1e-3+x*1e-3;
                Sp(i,j)=0;
            case 3
                %If SourceTerm varies both with T enter expresssion here
                %Sc,Sp=f(T)
                Sp(i,j)=-2*1e-3;
                Sc(i,j)=(3+9*T(i,j))*1e-5;
            case 4
                %If Source Term varies both with x,y & T enter expresssion here
                %Sc, Sp=f(x,T)
                Sc(i,j)=1e-4+(x(i,j)+y(i,j))*1e-3+(3+9*T(i,j))*1e-4;
                Sp(i,j)=-3*1e-5;
        end
    end
end

end

