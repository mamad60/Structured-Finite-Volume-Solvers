function [gamma]=DIF3(i,j)
%Returns the diffussion coefficient gamam @ ith cv centeroid
% i is the number of cv

global gStat Gamma0 XC YC Xmax Xmin Ymax Ymin kratio

switch gStat
    case 1 %contant
        gamma=Gamma0;
    case 2 %function of x/y
        %enter the expersion for gamma here gamm=f(Xc,Yc)
        %         if XC(i)>=(Xmin+0.45*(Xmax-Xmin)) && XC(i)<=(Xmin+0.55*(Xmax-Xmin))
        %             gamma=10*Gamma0;
        %         else
        %             gamma=Gamma0;
        %         end
        Lx=abs(Xmax-Xmin);
        Ly=abs(Ymax-Ymin);
        if (XC(i)>=0.25*Lx && XC(i)<=0.75*Lx) && (YC(j)>=0.25*Ly && YC(j)<=0.75*Ly)
            gamma=Gamma0*kratio;
        else
            gamma=Gamma0;
        end
        
end

  