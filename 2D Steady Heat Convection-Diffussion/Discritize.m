function  Discritize(I,J)
%Discritize the Convection-Diffiusion  EQ. @ interior cell
%Scheme(P)---> is a function of Peclet used for Convecftive Term Scheme Implemetation
global AW AE AN AS AP Q N DX DY
%=============================
%Auxillary for internal cells
i=I-1;
j=J-1;
IJ=(j-1)*N+i; %Convert 2D to 1D index
%==============================
%Compute Fluxes
[Fe,Fw,Fs,Fn,De,Dw,Dn,Ds]=Flux(I,J);  
%Compute Peclet numbers
Pe=Fe/De;
Pw=Fw/Dw;
Ps=Fs/Ds;
Pn=Fn/Dn;
%===============================================================
%Compute Coe.
AE(IJ)=(De*CSch(Pe)+max(-Fe,0))*DY(J);  %ae
AW(IJ)=(Dw*CSch(Pw)+max(Fw,0))*DY(J);   %aw
AN(IJ)=(Dn*CSch(Pn)+max(Fn,0))*DX(I);  %ae
AS(IJ)=(Ds*CSch(Ps)+max(-Fs,0))*DX(I);   %aw

%ap----It's Assumed tha Velocity field is divergence free(Continiuity EQ. applies)
AP(IJ)=AE(IJ)+AW(IJ)+AN(IJ)+AS(IJ)+(Fe-Fw+Fn-Fs);
%b
Q(IJ)=0; %Source Term


end

