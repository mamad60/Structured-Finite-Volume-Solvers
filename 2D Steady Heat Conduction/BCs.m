function BCs
%Sets Soluton th Matrix for the boundary point
%Customize for other Geometries

%Global Variables Definition
global X Y Dx Dy m n Sc Sp A b ke kw kn ks
global NL NR NB NT BCl BCr BCb BCt Tf

%Apply BCs
for j=2:n-1 %Left & Right Boundaries
    %Left Boundary
    np=(j-1)*m+1; %i=1
    switch NL
        case 0  % Fixed Temperature
            b(np)=BCl(j);
            A(np,np)=1;
        case 1  %Fixed Flux
            %aw=0;  %only for clarification
            %Calcalate k Source Terms Boundies
            %----------Set Coefficients
            Dxl=abs(X(2,j)-X(1,j));   %Distance First Two Points Adjacent to Left Boundary
            aI=ke(1,j)/Dxl; % Point next to Boundary
            aB=aI-Sp(1,j)*Dx(1,j);  % Boundary Point
            b(np)=Sc(1,j)*Dx(1,j)+BCl(j);
            %Assemble Soloution Matrix
            A(np,np)=aB;
            A(np,np+1)=-aI;
        case 2  %Convection Coeffient is Given
            %aw=0;  %only for clarification
            %Calcalate k Source Terms Boundies
            %----------Set Coefficients
            Dxl=abs(X(2,j)-X(1,j));   %Distance First Two Points Adjacent to Left Boundary
            aI=ke(1,j)/Dxl; % Point next to Boundary
            aB=aI-Sp(1,j)*Dx(1,j)+BCl(j);  % Boundary Point
            b(np)=Sc(1,j)*Dx(1,j)+BCl(j)*Tf;
            %Assemble Soloution Matrix
            A(np,np)=aB;
            A(np,np+1)=-aI;
    end

    %Right Boundary
    np=(j-1)*m+m; %i=m
    switch NR
        case 0  % Fixed Temperature
            b(np)=BCr(j);
            A(np,np)=1;
        case 1  %Fixed Flux
            %ae=0;  %only for clarification
            %Calcalate k Source Terms Boundies
            %----------Set Coefficients
            Dxr=abs(X(m,j)-X(m-1,j)); %Distance First Two Points Adjacent to Right Boundary
            aI=kw(m,j)/Dxr; % Point next to Boundary
            aB=aI-Sp(m,j)*Dx(m,j);  % Boundary Point
            b(np)=Sc(m,j)*Dx(m,j)+BCr(j);
            %Assemble Soloution Matrix
            A(np,np)=aB;
            A(np,np-1)=-aI;
        case 2  %Convection Coeffient is Given
            %ae=0;  %only for clarification
            %----------Set Coefficients
            Dxr=abs(X(m,j)-X(m-1,j)); %Distance First Two Points Adjacent to Right Boundary
            aI=kw(m,j)/Dxr; % Point next to Boundary
            aB=aI-Sp(m,j)*Dx(m,j)+BCr(j);  % Boundary Point
            b(np)=Sc(m,j)*Dx(m,j)+BCr(j)*Tf;
            %Assemble Soloution Matrix
            A(np,np)=aB;
            A(np,np-1)=-aI;
    end
end
     
    for i=1:m; %Top & Bottom Boudaries
        %Bottom Boundary
        np=(1-1)*m+i; %j=1
        switch NB
            case 0  % Fixed Temperature
                b(np)=BCb(i);
                A(np,np)=1;
            case 1  %Fixed Flux
                %an=0;  %only for clarification
                %----------Set Coefficients
                Dyt=abs(Y(i,2)-Y(i,1));   %Distance First Two Points Adjacent to Bottom Boundary
                aI=ks(i,1)/Dyt; % Point next to Boundary
                aB=aI-Sp(i,1)*Dy(i,1);  % Boundary Point
                b(np)=Sc(i,1)*Dy(i,1)+BCb(i);
                %Assemble Soloution Matrix
                A(np,np)=aB;
                A(np,np+m)=-aI;
            case 2  %Convection Coeffient is Given
                %an=0;  %only for clarification
                %----------Set Coefficients
                Dyt=abs(Y(i,2)-Y(i,1));   %Distance First Two Points Adjacent to Bottom Boundary
                aI=ks(i,1)/Dyt; % Point next to Boundary
                aB=aI-Sp(i,1)*Dy(i,1)+BCb(i);  % Boundary Point
                b(np)=Sc(i,1)*Dy(i,1)+BCb(i)*Tf;
                %Assemble Soloution Matrix
                A(np,np)=aB;
                A(np,np+m)=-aI;
        end
        %Top Boundary
        np=(n-1)*m+i; %j=n
        switch NT
            case 0  % Fixed Temperature
                b(np)=BCt(i);
                A(np,np)=1;
            case 1  %Fixed Flux
                %as=0;  %only for clarification
                %----------Set Coefficients
                Dyb=abs(Y(i,n)-Y(i,n-1)); %Distance First Two Points Adjacent to Top Boundary
                aI=kn(i,n)/Dyb; % Point next to Boundary
                aB=aI-Sp(i,n)*Dy(i,n);  % Boundary Point
                b(np)=Sc(i,n)*Dy(i,n)+BCt(i);
                %Assemble Soloution Matrix
                A(np,np)=aB;
                A(np,np-m)=-aI;
            case 2  %Convection Coeffient is Given
                %as=0;  %only for clarification
                Dyb=abs(Y(i,n)-Y(i,n-1)); %Distance First Two Points Adjacent to Top Boundary
                aI=kn(i,n)/Dyb; % Point next to Boundary
                aB=aI-Sp(i,n)*Dy(i,n)+BCt(i);  % Boundary Point
                b(np)=Sc(i,n)*Dy(i,n)+BCt(i)*Tf;
                %Assemble Soloution Matrix
                A(np,np)=aB;
                A(np,np-m)=-aI;
        end
    end

end

