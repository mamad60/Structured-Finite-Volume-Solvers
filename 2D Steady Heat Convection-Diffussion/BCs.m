function  BCs
%Invokes Boundary Conditions by Modifying Soure Term & Coefficient.


global  N M 
global NL NR NB NT BCl BCr BCb BCt 
global AW AE AN AS AP Q 

%--------Left & Right
for j=1:M;
    i=1; %Left Boundary
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch NL %Dirichlet BC
        case 0
            Q(IJ)=Q(IJ)+AW(IJ)*BCl(j); %FI(1,j);  BCs at first CV
            AW(IJ)=0;
        case 1    %Neuman BC on the Left
            Q(IJ)=Q(IJ)+BCl(j);%Ql  BCs at first CV
            AP(IJ)=AP(IJ)-AW(IJ);
            AW(IJ)=0;
    end
    i=N; %right Boundary
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch NR
        case 0
            Q(IJ)=Q(IJ)+AE(IJ)*BCr(j);  %FI(NI,j)BCs at first CV
            AE(IJ)=0;
        case 1    %Neuman BC on the Right
            Q(IJ)=Q(IJ)+BCr(j);  %Qr BCs at first CV
            AP(IJ)=AP(IJ)-AE(IJ);
            AE(IJ)=0;
    end
end
%--------Bottom & Top
for i=1:N;
    j=1;     %Bottom Boundary
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch NB  %Dirichlet BC
        case 0
            Q(IJ)=Q(IJ)+AS(IJ)*BCb(i);%FI(i,1) BCs at first CV
            AS(IJ)=0;
        case 1 %Neuman BC on the Bottom
            Q(IJ)=Q(IJ)+BCb(i);%   Qb  BCs at first CV
            AP(IJ)=AP(IJ)-AS(IJ);
            AS(IJ)=0;
    end
    j=M; %Top Boundary
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch NT
        case 0 %Dirichlet BC
            Q(IJ)=Q(IJ)+AN(IJ)*BCt(i);  %FI(i,NJ) BCs at first CV
            AN(IJ)=0;
        case 1 %Neuman BC on the Top
            Q(IJ)=Q(IJ)+BCt(i);   %Qt  %BCs at first CV
            AP(IJ)=AP(IJ)-AN(IJ);
            AN(IJ)=0;
    end
end

end

