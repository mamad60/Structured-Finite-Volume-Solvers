function  BiHarmonic3
%Calculates k @ the boundary by Biharmonic Averaging
%Returns k at the CV Faces

global  NIM NJM gStat Gamma0 Gammae Gammaw Gammas Gamman

%Check if Gamma Constant, So return it Quickly
if gStat==1
    g=Gamma0;
    Gammae(:)=g;
    Gammaw(:)=g;
    Gammas(:)=g;
    Gamman(:)=g;
    return
end
%X-Direction
for j=2:NJM
    for i=2:NIM     %Biharmonic Interpolation-Interior CVs
        gP=DIF3(i,j);
        gW=DIF3(i-1,j);
        gE=DIF3(i+1,j);
        Gammae(i,j)=(2*gP*gE)/(gP+gE);
        Gammaw(i,j)=(2*gP*gW)/(gP+gW);
    end
end
%Y-Direction
for i=2:NIM
    for j=2:NJM     %Biharmonic Interpolation-Interior CVs
        gP=DIF3(i,j);
        gS=DIF3(i,j-1);
        gN=DIF3(i,j+1);
        Gammas(i,j)=(2*gP*gS)/(gP+gS);
        Gamman(i,j)=(2*gP*gN)/(gP+gN);
    end
end
end

