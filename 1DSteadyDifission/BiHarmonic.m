function [ Gammae Gammaw ] = BiHarmonic()
%Calculates k @ the boundary by Biharmonic Averaging
%Returns k at the CV Faces
global NIM NI gStat
Gammae=zeros(1,NI);
Gammaw=zeros(1,NI);
%Check if Gamma Constant, So return it Quickly
if gStat==1
    g=DIF(1);
    Gammae(:)=g;
    Gammaw(:)=g;
    return
end

for i=3:NIM-1     %Biharmonic Interpolation-Interior CVs
    gP=DIF(i);
    gW=DIF(i-1);
    gE=DIF(i+1);
    Gammae(i)=(2*gP*gE)/(gP+gE);
    Gammaw(i)=(2*gP*gW)/(gP+gW);
end
%Two Boubdary CVS
%West CV
gP=DIF(2);
gE=DIF(3);
Gammaw(2)=gP;
Gammae(2)=(2*gP*gE)/(gP+gE);
%East CV
gP=DIF(NIM);
gW=DIF(NIM-1);
Gammae(NIM)=gP;
Gammaw(NIM)=(2*gP*gW)/(gP+gW);
end

