function [Su Sp]=Source(sStat)
%Return Source Term S=Sp*FIp*Su
global XC  YC Su0 Sp0 DX DY NM N  NIM NJM NI NJ
global Xmin Xmax  Ymin Ymax

U=zeros(NI,NJ);  % Source Term---Su
P=zeros(NI,NJ);  % Source Term---Sp
Su=zeros(NM,1);  % Source Term---Su
Sp=zeros(NM,1);  % Source Term---Sp
Lx=abs(Xmax-Xmin);
Ly=abs(Ymax-Ymin);
switch sStat
    case 1 %zero Source Term
        %Do nothing
        return
    case 2 %Constant Sourc Term;
        for  i=1:NI;
            for j=2:NJ
                U(i,j)=Su0*DX(i)*DY(j);
                P(i,j)=Sp0*DX(i)*DY(j);
            end
        end
    case 3  % Su or Sp is function of X/Y
        for  i=1:NI;
            for j=1:NJ
                if (XC(i)>=0.25*Lx && XC(i)<=0.75*Lx) && (YC(j)>=0.25*Ly && YC(j)<=0.75*Ly)
                    U(i,j)=Su0*DX(i)*DY(j);
                end
                P(i,j)=0;
            end
        end
end
%Fill into 1D Array
for i=2:NIM
    for j=2:NJM
        I=i-1;
        J=j-1;
        IJ=(J-1)*N+I; %Convert 2D to 1D index
        Sp(IJ)=P(i,j);
        Su(IJ)=U(j,j);
    end
end
if sStat==3
    figure
    [x1 , y1]=ndgrid(XC,YC);
    contourf(x1,y1,U,2)
    title('Source Term')
end
