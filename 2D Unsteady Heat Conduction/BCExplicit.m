function BCExplicit
%Boundary Conditons are applied in this functions-Explicit Unsteady Version
%Customize According to the problem

%Global Variables Definition
%Global Variables Definition
global X Y Dx Dy m n Sc Sp  b ke kw kn ks
global NL NR NB NT BCl BCr BCb BCt Tf T

%Apply BCs
for j=1:n %Left & Right Boundaries
    %Left Boundary
    switch NL
        case 0  % Fixed Temperature
            T(1,j)=BCl(j);
        case 1  %Fixed Flux
            %Calcalate k Source Terms Boundies
            %----------Set Coefficients
            Dxl=abs(X(2,j)-X(1,j));   %Distance First Two Points Adjacent to Left Boundary
            aI=ke(1,j)/Dxl; % Point next to Boundary
            aB=aI-Sp(1,j)*Dx(1,j);  % Boundary Point
            b=Sc(1,j)*Dx(1,j)+BCl(j);
            %Solve for T  aBTB=AITI+b
            T(1,j)=(aI*T(2,j))/aB;
        case 2  %Convection Coeffient is Given
            %Calcalate k Source Terms Boundies
            %----------Set Coefficients
            Dxl=abs(X(2,j)-X(1,j));   %Distance First Two Points Adjacent to Left Boundary
            aI=ke(1,j)/Dxl; % Point next to Boundary
            aB=aI-Sp(1,j)*Dx(1,j)+BCl(j);  % Boundary Point
            b=Sc(1,j)*Dx(1,j)+BCl(j)*Tf;
            %Solve for T  aBTB=AITI+b
            T(1,j)=(aI*T(2,j))/aB;
    end
    
    %Right Boundary
    switch NR
        case 0  % Fixed Temperature
            T(m,j)=BCr(j);
        case 1  %Fixed Flux
            %Calcalate k Source Terms Boundies
            %----------Set Coefficients
            Dxr=abs(X(m,j)-X(m-1,j)); %Distance First Two Points Adjacent to Right Boundary
            aI=kw(m,j)/Dxr; % Point next to Boundary
            aB=aI-Sp(m,j)*Dx(m,j);  % Boundary Point
            b=Sc(m,j)*Dx(m,j)+BCr(j);
            %Solve for T  aBTB=AITI+b
            T(m,j)=(aI*T(m-1,j))/aB;
        case 2  %Convection Coeffient is Given
            %----------Set Coefficients
            Dxr=abs(X(m,j)-X(m-1,j)); %Distance First Two Points Adjacent to Right Boundary
            aI=kw(m,j)/Dxr; % Point next to Boundary
            aB=aI-Sp(m,j)*Dx(m,j)+BCr(j);  % Boundary Point
            b=Sc(m,j)*Dx(m,j)+BCr(j)*Tf;
            %Solve for T  aBTB=AITI+b
            T(m,j)=(aI*T(m-1,j))/aB;
    end
end

for i=2:m-1; %Top & Bottom Boudaries
    %Bottom Boundary
    switch NB
        case 0  % Fixed Temperature
            T(i,1)=BCb(i);
        case 1  %Fixed Flux
            %----------Set Coefficients
            Dyt=abs(Y(i,2)-Y(i,1));   %Distance First Two Points Adjacent to Top Boundary
            aI=ks(i,1)/Dyt*Dx(i,1); % Point next to Boundary
            aB=aI-Sp(i,1)*Dx(i,1)*Dy(i,1);  % Boundary Point
            b=Sc(i,1)*Dx(i,1)*Dy(i,1)+BCb(i);
            %Solve for T  aBTB=AITI+b
            T(i,1)=(aI*T(i,2))/aB;
        case 2  %Convection Coeffient is Given
            %----------Set Coefficients
            Dyt=abs(Y(i,2)-Y(i,1));   %Distance First Two Points Adjacent to Left Boundary
            aI=ks(i,1)/Dyt*Dx(i,1); % Point next to Boundary
            aB=aI-Sp(i,1)*Dx(i,1)*Dy(i,1)+BCb(i);  % Boundary Point
            b=Sc(i,1)*Dx(i,1)*Dy(i,1)+BCb(i)*Tf;
            %Solve for T  aBTB=AITI+b
            T(i,1)=(aI*T(i,2))/aB;
    end
    %Top Boundary
    switch NT
        case 0  % Fixed Temperature
            T(i,n)=BCt(i);
        case 1  %Fixed Flux
            %----------Set Coefficients
            Dyb=abs(Y(i,n)-Y(i,n-1)); %Distance First Two Points Adjacent to Top Boundary
            aI=kn(i,n)/Dyb*Dx(i,n); % Point next to Boundary
            aB=aI-Sp(i,n)*Dx(i,n)*Dy(i,n);  % Boundary Point
            b=Sc(i,n)*Dx(i,n)*Dy(i,n)+BCt(i);
            %Solve for T  aBTB=AITI+b
            T(i,n)=(aI*T(i,n-1))/aB;
        case 2  %Convection Coeffient is Given
            Dyb=abs(Y(i,n)-Y(i,n-1)); %Distance First Two Points Adjacent to Top Boundary
            aI=kn(i,n)/Dyb*Dx(i,n); % Point next to Boundary
            aB=aI-Sp(i,n)*Dx(i,n)*Dy(i,n)+BCt(i);  % Boundary Point
            b=Sc(i,n)*Dx(i,n)*Dy(i,n)+BCt(i)*Tf;
            %Solve for T  aBTB=AITI+b
            T(i,n)=(aI*T(i,n-1))/aB;
    end
end

end

