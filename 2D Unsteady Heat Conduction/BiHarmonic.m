function [ ke kw ks kn ] = BiHarmonic(kStat)
%Calculates k @ the boundary by Biharmonic Averaging
%Returns k at the CV Faces
global X Y Xvol Yvol Xe Xw Yn Ys   m n  Told 

ke=zeros(size(X));
kw=zeros(size(X));
kn=zeros(size(X));
ks=zeros(size(X));
%Check if k or Constant, So return it Quickly
if kStat==0
    k=ThermalConductivity(kStat,X(2,2),Y(2,2),Told(2,2));
    ke(1:m-1,:)=k;
    kw(2:m,:)=k;
    ks(:,1:n-1)=k;
    kn(:,2:n)=k;
    return
end

for i=1:size(X,1)
    for j=1:size(X,2)
        %------Compute k @ Grid Points
        kP=ThermalConductivity(kStat,X(i,j),Y(i,j),Told(i,j));
        if i~=1
            kW=ThermalConductivity(kStat,X(i-1,j),Y(i-1,j),Told(i-1,j));
        end
        if i~=m
            kE=ThermalConductivity(kStat,X(i+1,j),Y(i+1,j),Told(i+1,j));
        end
        if j~=1
            kN=ThermalConductivity(kStat,X(i,j-1),Y(i,j-1),Told(i,j-1));
        end
        if j~=n
            kS=ThermalConductivity(kStat,Y(i,j+1),Y(i,j+1),Told(i,j+1));
        end
        
        %Compute k @ CV Faces by BiHarmonic Averaging
        %------------------------ke
        if i~=m
            Xep=abs(Xvol(i,j)-X(i,j));
            fe=Xep/Xe(i,j);
            ke(i,j)=1/((1-fe)/kP+fe/kE);
        else
            ke(i,j)=0;
        end
        %------------------------kw
        if i~=1
            Xwp=abs(X(i,j)-Xvol(i-1,j));
            fw=Xwp/Xw(i,j);
            kw(i,j)=1/((1-fw)/kP+fw/kW);
        else
            kw(i,j)=0;
        end
        %------------------------kn
        if j~=1
            Ynp=abs(Y(i,j)-Yvol(i,j-1));
            fn=Ynp/Yn(i,j);
            kn(i,j)=1/((1-fn)/kP+fn/kN);
        else
            kn(i,j)=0;
        end
        %------------------------ks
        if j~=n
            Ysp=abs(Y(i,j)-Yvol(i,j));
            fs=Ysp/Ys(i,j);
            ks(i,j)=1/((1-fs)/kP+fs/kS);
        else
            ks(i,j)=0;
        end
    end

end
end

