function AdjustPost
%This function Adjusts and reshapes the solution for Post-Processing
global FI F N M NI NJ NIM NJM
global Bw Be Bs Bn
global BCw BCe BCn BCs XC YC

%==========================================================================
F=reshape(F,N,M); %Convert to 2D Array
%Fill into solution Varibles and Interpolate IF Neaded
FI(2:NIM,2:NJM)=F;
%==========================================================================
%Fill in FI Rows with Boundary Conditions
if Bw==0,FI(1,:)=BCw; end  %West BC
if Be==0,FI(NI,:)=BCe;end  %East BC
if Bs==0,FI(:,1)=BCs; end  %South Boundary
if Bn==0,FI(:,NJ)=BCn;end  %North Boundary
%==========================================================================
% Interpolate FI @ Boundaries with specified q
%West & East Boundary interpolation(IF Flux BCs specified)
for j=1:NJ;
    %West Boundary
    if Bw==1, FI(1,j)=BCw(j)*(XC(1)-XC(2))+FI(2,j);end
    %East Boundary
    %modified for the problem
    if j>0.5*NJ, FI(NI,j)=FI(NIM,j);end %Insulated
end
%Top & Bottom Boundary interpolation(IF Flux BCs specified)
for i=1:NI;
    %South Boundary
    if Bs==1, FI(i,1)=BCs(i)*(YC(1)-YC(2))+FI(i,2);end
    %North Boundary
    %modified for the problem
    if i>0.5*NI, FI(i,NJ)=FI(i,NJM);end %insulated
end
%==========================================================================
% %Correct 4 Corners
% %---South-West Corner
% FI(1,1)=0.5*(FI(1,2)+FI(2,1));
% %---South-West Corner
% FI(NI,1)=0.5*(FI(NIM,1)+FI(NI,2));
% %---North-West Corner
% FI(1,NJ)=0.5*(FI(2,NJ)+FI(1,NJM));
% %---North-East Corner
% FI(NI,NJ)=0.5*(FI(NIM,NJ)+FI(NI,NJM));
% % 









end