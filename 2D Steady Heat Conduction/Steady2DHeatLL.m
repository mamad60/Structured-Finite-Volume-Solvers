% This Scripts Solves 2D Steady Heat Transfer Equation
% d/dx(kdT/dx)+d/dy(kdT/dy)+S=0, From Patankar Book, Section 4-4
%Line by Line Solotion is Employed
% Both Thermal Conductivity and Source therms can be varied with time & X,Y
%Application of the iterative Matrix Solvers---Single Grid Version
%Line by line Sweeping               N
%Discritization Stencil   |--->  W---P---E
%                                    S
%By Mohammad Aghakhani, 2010,for Teaching purposes on the CFD Class
%Feel Free using for Educaitonal & Research Purposes

clc
clear
close all

%Please set the Folder for Iterative Solvers below or Copy files to the Current directry
iFolder='../../Iterative Solvers';
%Check & Set Path for Solver Directery
pathCell = regexp(path, pathsep, 'split');
if ispc  % Windows is not case-sensitive
    onPath = any(strcmpi(iFolder, pathCell));
else
    onPath = any(strcmp(iFolder, pathCell));
end
if ~onPath,addpath(genpath(iFolder)); end
%Global Varibles Definintion
global X Y Xvol Yvol Xe Xw Yn Ys Dx Dy m n T Told kStat
%-------------------------------------

%Inputs
Lx=2;  % Lenght of Domain in X direction(i)
Ly=1;  % Lenght of Domain in X direction(i)
m=20; % Number of Control Volumes in X direction(i)
n=20; % Number of Control Volumes in Y direction(j)
kStat=3;  % Themal Conductivity Calculation Method
% KStat=0:Constant 1:Function of X(and/Or)Y 2:Function of Temperature 3: Function of Both Space & Temp.
sStat=0;  %Source Term Calculation Method
%0: No Source Term 1:Constant 2:Function of X(and/Or)Y 3:Function of Temperature 4: Function of Both
iMethod=1; % Solution Method:
%0 MATLAB Direct Solver 1: TDMA(Direct)
%2:Jacobi(Matrix) 3:Jacobi 4:Gauss-Seidel(Matrix) 5:Gauss-Seidel
%6: Symmetric Gauss-Seidel Matrix 7:Symmetric Gauss-Seidel
%8:SOR(Matrix) 9:SOR 10:Symmetric SOR Matrix
%11:Red-Black Gauss-Seidel  12:SOR Solver 4 Triangular matrix
%13:Conjugate Gradient   14: Preconditioned Conjugate Gradient
%15: BiConjugate gradient
espSolver=1e-3; % Maximum Residual of Matrix Iterative Solver
MaxITSolver=100000; % Maximum Iteraton of Matrix Iterative Solver

%-----Boundary Conditions
NL=1; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Left Boundary
NR=1; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Right Boundary
NB=0; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Bottom Boundary
NT=0; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Top Boundary
Tl=373; % Temperature at The Left Boundary, Applies only if NL=0
Tr=273; % Temperature at The Right Boundary, Applies only if NR=0
Tb=273; % Temperature at The Bottom Boundary, Applies only if NB=0
Tt=373; % Temperature at The Top Boundary, Applies only if NT=0
Ql=0; % Flux at The Left Boundary, Applies only if NL=1
Qr=0; % Flux at The Right Boundary, Applies only if NR=1
Qb=0; % Flux at The Bottom Boundary, Applies only if NB=1
Qt=0; % Flux at The Top Boundary, Applies only if NT=1
%Flux Normal to boundary & Towards it Assumed Positive sign
hl=0; % Convection Coefficient at The Left Boundary, Applies only if NL=2
hr=0.1; % Flux at The Left Boundary, Applies only if NL=2
hb=0; % Convection Coefficient at The Bottom Boundary, Applies only if NB=2
ht=0.1; % Flux at The Top Boundary, Applies only if NT=2
Tf=273.0; %Temperature @ infinity used if NL=2 Or NR=2 q=h(Tf-TB)
%------------------Sweep Settings
MaxSi=1000; % Maximum  Number of Sweeps in X direction
epsSi=1e-6; %Convergence critreia for Sweeps in X direction

%------------Only if Thermal Conductivity Or Source Term are function of
%Temperature ,variables Defined in this Section are used
MaxIT=100;  %Maximum Allowed Iteration
epsT=1e-6;   %Convergence Criteria in Two Concecutive Iterations
%-----------------------------------------------------------
%Variable Allocation
X=zeros(m,n); %X Location of Verices
Y=zeros(m,n); %Y Location of Verices
Xvol=zeros(m,n); % Location of Control Volume Faces
Yvol=zeros(m,n); % Location of Control Volume Faces
T=zeros(m,n); %Variable for Storing Solution
Told=zeros(m,n); %Variable for Storing Old Iteration Solution
A=zeros(n); % Coefficient Matrix---Line by line Solution-X Direction
b=ones(n,1); % Right hand side Matrix
Dx=zeros(m,n); % Control Volume lenght for Source Terms Evaluation
Dy=zeros(m,n); % Control Volume lenght for Source Terms Evaluation
Xe=zeros(m,n); %Distance Between P and E point
Xw=zeros(m,n); %Distance Between P and W point
Yn=zeros(m,n); %Distance Between P and N point
Ys=zeros(m,n); %Distance Between P and S point
ke=zeros(m,n); % k @ boundaries
kw=zeros(m,n); 
kn=zeros(m,n);
ks=zeros(m,n); 
qx=zeros(m,n);
Sc=zeros(m,n);  %Source Term Sc  @  Grid Nodes
Sp=zeros(m,n);  %Source Term Sp @  Grid Nodes
qy=zeros(m,n);
q_l=zeros(n,1);  %Flux @ Boundary Points
q_r=zeros(n,1);
q_t=zeros(1,m);
q_b=zeros(1,m);
err=1000; % Varible for Storing Error for Sweep loop
err_out=1000; % Varible for Storing Error for Sweep loop
great=1e10; % Assumsion for infinity
%------------------------------------------------------------

%Grid Generation & Geometrical Calculations
Grid2d(Lx,Ly,m,n);
%First check if an iterative solution must be applied
%for Temperature Depencency of Thermal Conductinvity or Source Term
tIT=MaxIT;
MaxIT=1;
iflag=0;
if kStat==2 || kStat==3 || sStat==3 || sStat==4
    MaxIT=tIT;
    iflag=1;
end
%Allocate Matrix for storing Errors for Post-Processing and Initiale T
T(:)=0.5*(Tl+Tr);
%----------------------------------------------------------------------
error=zeros(1,MaxSi);
scrsz = get(0,'ScreenSize');

nS=1; %Number of Sweeps
IT=1; %Outer loop Counter
%% Outer Loop : Perfomed more than one if k(t) or s(x)
tic;
while (IT<=MaxIT) && (err_out>epsT)
    if iflag
        Told_out=T; 
        fprintf('\nBegining Outer Loop Iteration=%i\n\n',IT);
    end
    while (nS<=MaxSi) && (err>epsSi) %Sweep Itration
        Told=T;
        %Compute Source Terms
        [Sc,Sp] = SourceTerm(sStat,X,Y,Told);
        %Calcalate Thermal Conductivity by BiHarmonic Averaging
        [ ke kw ks kn ] = BiHarmonic(kStat);
        %-Apply Bcs @ Left & Right vertical line by the Explicit Method
        for   j=1:n;
            %Left Boundary
            switch NL
                case 0  % Fixed Temperature
                    T(1,j)=Tl;
                case 1  %Fixed Flux
                    Dxl=abs(X(2,j)-X(1,j));   %Distance First Two Points Adjacent to Left Boundary
                    T(1,j)=(Dxl*(Ql+Sc(1,j)*Dx(1,j))+ke(1,j)*Told(2,j))/(ke(1,j)-Sp(1,j)*Dx(1,j)*Dxl);
                case 2  %Convection Coeffient is Given
                    Dxl=abs(X(2,j)-X(1,j));   %Distance First Two Points Adjacent to Left Boundary
                    T(1,j)=(Dxl*(hl*Tf+Sc(1,j)*Dx(1,j))+ke(1,j)*Told(2,j))/(ke(1,j)+Dxl*(hl-Sp(1,j)*DX(1,j)));
            end
            %Right Boundary
            switch NR
                case 0  % Fixed Temperature
                    T(m,j)=Tr;
                case 1  %Fixed Flux
                    Dxr=abs(X(m,j)-X(m-1,j));   %Distance First Two Points Adjacent to Right Boundary
                    T(m,j)=(Dxr*(Qr+Sc(m,j)*Dx(m,j))+kw(m,j)*Told(m-1,j))/(kw(m,j)-Sp(m,j)*Dx(m,j)*Dxr);
                case 2  %Convection Coeffient is Given
                    Dxr=abs(X(m,j)-X(m-1,j));   %Distance First Two Points Adjacent to Right Boundary
                    T(m,j)=(Dxr*(hr*Tf+Sc(m,j)*Dx(m,j))+kw(m,j)*Told(m-1,j))/(kw(m,j)+Dxr*(hr-Sp(m,j)*Dx(m,j)));
            end
        end
        
        
        %Left To Right Swip(EQ. is solved @ vertical lines and moved to Right)
        for i=2:m-1
            %Set Bcs on Two end points of the line j=1 & j=n
            %Bottom Boundary
            switch NB
                case 0  % Fixed Temperature
                    b(1)=Tt;
                    A(1,1)=1;
                case 1  %Fixed Flux
                    an=0;  %only for clarification
                    %----------Set Coefficients
                    Dyt=abs(Y(i,2)-Y(i,1));   %Distance First Two Points Adjacent to Top Boundary
                    aI=ks(i,1)/Dyt*Dx(i,1); % Point next to Boundary
                    aB=aI-Sp(i,1)*Dx(i,1)*Dy(i,1);  % Boundary Point
                    b(1)=Sc(i,1)*Dx(i,1)*Dy(i,1)+Qt;
                    %Assemble Soloution Matrix
                    A(1,1)=aB;
                    A(1,2)=-aI;
                case 2  %Convection Coeffient is Given
                    an=0;  %only for clarification
                    %----------Set Coefficients
                    Dyt=abs(Y(i,2)-Y(i,1));   %Distance First Two Points Adjacent to Left Boundary
                    aI=ks(i,1)/Dyt*Dx(i,1); % Point next to Boundary
                    aB=aI-Sp(i,1)*Dx(i,1)*Dy(i,1)+ht;  % Boundary Point
                    b(1)=Sc(i,1)*Dx(i,1)*Dy(i,1)+ht*Tf;
                    %Assemble Soloution Matrix
                    A(1,1)=aB;
                    A(1,2)=-aI;
            end
            %Top Boundary
            switch NT
                case 0  % Fixed Temperature
                    b(n)=Tb;
                    A(n,n)=1;
                case 1  %Fixed Flux
                    as=0;  %only for clarification
                    %----------Set Coefficients
                    Dyb=abs(Y(i,n)-Y(i,n-1)); %Distance First Two Points Adjacent to Right Boundary
                    aI=kn(i,n)/Dyb*Dx(i,n); % Point next to Boundary
                    aB=aI-Sp(i,n)*Dx(i,n)*Dy(i,n);  % Boundary Point
                    b(n)=Sc(i,n)*Dx(i,n)*Dy(i,n)+Qb;
                    %Assemble Soloution Matrix
                    A(n,n)=aB;
                    A(n,n-1)=-aI;
                case 2  %Convection Coeffient is Given
                    as=0;  %only for clarification
                    Dyb=abs(Y(i,n)-Y(i,n-1)); %Distance First Two Points Adjacent to Right Boundary
                    aI=kn(i,n)/Dyb*Dx(i,n); % Point next to Boundary
                    aB=aI-Sp(i,n)*Dx(i,n)*Dy(i,n)+hb;  % Boundary Point
                    b(n)=Sc(i,n)*Dx(i,n)*Dy(i,n)+hr*Tf;
                    %Assemble Soloution Matrix
                    A(n,n)=aB;
                    A(n,n-1)=-aI;
            end
            for j=2:n-1
                %Discritize on the vertical line
                %Compute Coefficients
                %ae
                ae=ke(i,j)/Xe(i,j)*Dy(i,j);
                %aw
                aw=kw(i,j)/Xw(i,j)*Dy(i,j);
                %an
                an=kn(i,j)/Yn(i,j)*Dx(i,j);
                %as
                as=ks(i,j)/Yn(i,j)*Dx(i,j);
                %ap
                ap=an+as+ae+aw-Sp(i,j)*Dx(i,j)*Dy(i,j);
                %In Line by line Iteration T from old Iteration is used @ E & W
                c=ae*T(i+1,j)+aw*T(i-1,j);
                %b
                b(j)=Sc(i,j)*Dx(i,j)*Dy(i,j)+c;
                %Assemble Soloution Matrix
                A(j,j)=ap;
                A(j,j-1)=-an;
                A(j,j+1)=-as;
            end
            %Solve Equations
            %A=sparse(A); Improves solution time for Matrix Solvers
            T1=Solve(A,T(i,:)',b,MaxITSolver,espSolver,iMethod);
            T(i,:)=T1;
            %Show how Sweep Propegate information form boundaries
            %         contourf(T)
            %         colorbar
        end
        err=abs(norm((T-Told)));
        fprintf('\nNumber of Sweeps=%i\t Error=%2.6e\n',nS,err);
        if err>=great,break; end
        error(nS)=err;
        %         %Show how Sweep Propegate information form boundaries
        %         contourf(T)
        nS=nS+1;
    end
    if iflag
        err_out=abs(norm((T-Told_out)));
        fprintf('\n\nCompleted Outer Loop Iteration=%i\t Error=%2.6e\n\n',IT,err_out);
        err=1000;
        nS=1;
        if err_out>=great,break; end
    end
    IT=IT+1;
end

toc
% Check Convergence & Plot Convergence History
nS=nS-1;
IT=IT-1;
if(err>great || err_out>great)
    fprintf(1,'Iterations Diverged\n');
    fprintf(1,'Please Consider Change in Solution Parmeters and Run the Code Again\n');
    disp('Press any key')
    pause
    return;
end
figure('Name','Sweep Convergence History','NumberTitle','off','OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)/2]);
fprintf(1,'Calculated by %i Sweeps\n',nS);
semilogy(1:length(error),error,'- r');
xlabel('Sweeps');
ylabel('Error');
title('Sweep Convergence History');


%% Post-Processing-----------------------------------------------------------

%---Report Average Temeratures @ Right , Left, Top & Bottom Boundaries
fprintf('\nThe Average Temperature @ the Left Boundary is: %2.2e\n',mean(T(1,:)));
fprintf('The Average Temperature @ the Right Boundary is: %2.2e\n',mean(T(m,:)));
fprintf('The Average Temperature @ the Top Boundary is: %2.2e\n',mean(T(:,n)));
fprintf('The Average Temperature @ the Bottom Boundary is: %2.2e\n',mean(T(1,:)));

%Plot results
%----Plot  T Contour
figure('Name','Temperature Countor','NumberTitle','off','OuterPosition',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
contourf(X,Y,T,'LineWidth',2);
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Temperature Contour');
colorbar
%----Surface of Temperature
figure('Name','Surface of Temperature','NumberTitle','off','OuterPosition',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
surf(X,Y,T);
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('Surface of Temperature');

%---Compute Fluxes @ Boundary Points
%Left & Right Boundaries
for j=1:n;
    q_l(j)=ke(1,j)*(T(1,j)-T(2,j))/abs(X(1,j)-X(2,j))*Dy(1,j);   %Flux @ the Left Boundary
    q_r(j)=kw(m,j)*(T(m,j)-T(m-1,j))/abs(X(m-1)-X(m,j))*Dy(m,j);   %Flux @ the Left Boundary
end
%Top & Bottom Boundaries
for i=1:m;
    q_t(i)=kn(i,n)*(T(i,n)-T(i,n-1))/abs(Y(i,n)-Y(i,n-1))*Dx(i,n);   %Flux @ the Bottom Boundary
    q_b(i)=ks(i,1)*(T(i,1)-T(i,2))/abs(Y(i,1)-Y(i,2))*Dx(i,2);   %Flux @ the Top Boundary
end
%-----Compute and Plot Fluxes @ East & South CV Faces
for j=1:n;
    for i=1:m-1;
        qx(i,j)=ke(i,j)*(T(i,j)-T(i+1,j))/abs(X(i+1,j)-X(i,j));
    end
    qx(m,j)=-q_r(j)/Dy(m,j);
end
for i=1:m;
    for j=1:n-1;
        qy(i,j)=ks(i,j)*(T(i,j)-T(i,j+1))/abs(Y(i,j+1)-Y(i,j));
    end
    qy(i,n)=-q_b(i)/Dx(i,2);    
end


%Plot Flux vectors
figure('Name','Flux Vectors','NumberTitle','off','OuterPosition',[1 1 scrsz(3)/2 scrsz(4)/2]);
contour(X,Y,T)
hold on
quiver(X,Y,qx,qy)
hold off

%Report Flux in inbalence @ Boundaries
fprintf('\nFlux Inbalance @ Right & Left Boundaries is %2.2e\n',sum(q_r)+sum(q_l));
fprintf('Flux Inbalance @ Top & Bottom Boundaries is %2.2e\n',sum(q_t)+sum(q_b));


%---if k is a function of x Plot its Variation @ CV Faces
if kStat==1 || kStat==2 || kStat==3
    figure('Name','Thermal Conductivity Contour','NumberTitle','off');
    contourf(X,Y,ThermalConductivity(kStat,X,Y,T),'LineWidth',2);
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Thermal Counductiviy ')
    colorbar
end


fprintf('\nGood Lock, Mohammad Aghakhani,2012\n');










