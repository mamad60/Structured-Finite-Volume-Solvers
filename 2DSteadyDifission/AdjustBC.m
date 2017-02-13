function  AdjustBC
%Invokes Boundary Conditions by Modifying Soure Term & Coefficient.

global NIM NJM M N
global Su Sp aE aW aN aS
global Bw Be Bn Bs BCw BCe BCn BCs DX DY

%West & East Boundaries
for j=1:M;
    i=1; %Left Boundary
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch Bw   %West Boundary(Left)
        case 0
            Sp(IJ)=Sp(IJ)-aW(IJ);
            Su(IJ)=Su(IJ)+aW(IJ)*BCw(j);
            aW(IJ)=0;
        case 1
            Sp(IJ)=Sp(IJ)+0;
            Su(IJ)=Su(IJ)+BCw(j)*DY(j+1);
            aW(IJ)=0;
    end
    i=N; %Right Boundary
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch Be   %East Boundary(Right)
        case 0
            Sp(IJ)=Sp(IJ)-aE(IJ);
            Su(IJ)=Su(IJ)+aE(IJ)*BCe(j);
            aE(IJ)=0;
        case 1
            Sp(IJ)=Sp(IJ)+0;
            Su(IJ)=Su(IJ)+BCe(j)*DY(j+1);
            aE(IJ)=0;
    end
end
% %North & South Boundaries
for i=1:N;
    j=M;     %North Boundary
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch Bn   %North Boundary(Top)
        case 0
            Sp(IJ)=Sp(IJ)-aN(IJ);
            Su(IJ)=Su(IJ)+aN(IJ)*BCn(i);
            aN(IJ)=0;
        case 1
            Sp(IJ)=Sp(IJ)+0;
            Su(IJ)=Su(IJ)+BCn(i)*DX(i+1);
            aN(IJ)=0;
    end
    j=1; %South Boundary
    IJ=(j-1)*N+i; %Convert 2D to 1D index
    switch Bs   %South Boundary(Bottom)
        case 0
            Sp(IJ)=Sp(IJ)-aS(IJ);
            Su(IJ)=Su(IJ)+aS(IJ)*BCs(i);
            aS(IJ)=0;
        case 1
            Sp(IJ)=Sp(IJ)+0;
            Su(IJ)=Su(IJ)+BCs(i)*DX(i+1);
            aS(IJ)=0;
    end
end

end

