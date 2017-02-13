function [gamma]=DIF(i)
%returns the diffussion coefficient gamam @ ith cv centeroid
% i is the number of cv

global gStat Gamma0 XC Xmax Xmin
switch gStat
    case 1 %contant
        gamma=Gamma0;
    case 2 %function of x
         %enter the expersion for gamma here gamm=f(xc)
        if XC(i)>=(Xmin+0.45*(Xmax-Xmin)) && XC(i)<=(Xmin+0.55*(Xmax-Xmin))
            gamma=10*Gamma0;
        else
            gamma=Gamma0;
        end
end

end

