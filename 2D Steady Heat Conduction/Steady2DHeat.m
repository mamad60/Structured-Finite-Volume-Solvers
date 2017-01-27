% This Scripts Solves 2D Steady Heat Transfer Equation
% d/dx(kdT/dx)+d/dy(kdT/dy)+S=0, From Patankar Book, Section 4-4
%Direct Solution via Relabeling Method is Employed
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
global X Y Xvol Yvol Xe Xw Yn Ys Dx Dy m n T Told kStat Sc Sp A b ke kw kn ks
%-------------------------------------

%Inputs
Lx=2;  % Lenght of Domain in X direction(i)
Ly=1;  % Lenght of Domain in X direction(i)
m=10; % Number of Control Volumes in X direction(i)
n=10; % Number of Control Volumes in Y direction(j)
kStat=3;  % Themal Conductivity Calculation Method
% KStat=0:Constant 1:Function of X(and/Or)Y 2:Function of Temperature 3: Function of Both Space & Temp.
sStat=4;  %Source Term Calculation Method
%0: No Source Term 1:Constant 2:Function of X(and/Or)Y 3:Function of Temperature 4: Function of Both
iMethod=13; % Solution Method:
%0 MATLAB Direct Solver 
%2:Jacobi(Matrix) 3:Jacobi 4:Gauss-Seidel(Matrix) 5:Gauss-Seidel
%6: Symmetric Gauss-Seidel Matrix 7:Symmetric Gauss-Seidel
%8:SOR(Matrix) 9:SOR 10:Symmetric SOR Matrix
%11:Red-Black Gauss-Seidel  12:SOR Solver 4 Triangular matrix
%13:Conjugate Gradient   14: Preconditioned Conjugate Gradient
%15: BiConjugate gradient
espSolver=1e-6; % Maximum Residual of Matrix Iterative Solver
MaxITSolver=100000; % Maximum Iteraton of Matrix Iterative Solver

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
T1=zeros(m*n,1); % T in Relabled form
%A=spalloc(m*n,n*m,5*(m-1)*(n-1)+4*(m+n)); % Coefficient Matrix---Sparse Matrix Allocation
A=sparse(n*m); % Coefficient Matrix
b=zeros(n*m,1); % Right hand side Matrix
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
SetBC   % Set boundary Conditions for the problem
%Allocate Matrix for storing Errors for Post-Processing and Initiale T
T(:)=273;
%----------------------------------------------------------------------
error=zeros(1,MaxSi);
scrsz = get(0,'ScreenSize');

nS=1; %Number of Sweeps
IT=1; %Outer loop Counter
%% Outer Loop : Perfomed more than one if k(t) or S(x,y)
tic;
while (IT<=MaxIT) && (err>epsT)
    Told=T;
    %Compute Source Terms
    [Sc,Sp] = SourceTerm(sStat,X,Y,Told);
    %Calcalate Thermal Conductivity by BiHarmonic Averaging
    [ ke kw ks kn ] = BiHarmonic(kStat);

    %Apply BCs
    BCs()   %Customize this section for changing Problem
    %Interior Cells
    for i=2:m-1
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
            as=ks(i,j)/Ys(i,j)*Dx(i,j);
            %ap
            ap=an+as+ae+aw-Sp(i,j)*Dx(i,j)*Dy(i,j);
            % Set Labaling Coe.
            np=(j-1)*m+i; 
            %b
            b(np)=Sc(i,j)*Dx(i,j)*Dy(i,j);
            %Assemble Soloution Matrix
            A(np,np)=ap;
            A(np,np+1)=-ae;
            A(np,np-1)=-aw;
            A(np,np-m)=-an;
            A(np,np+m)=-as;
        end
    end
    %Solve Equations
    %A=sparse(A); Improves solution time for Matrix Solvers
    T1=reshape(T,m*n,1);
    T1=Solve(A,T1,b,MaxITSolver,espSolver,iMethod);
    T=reshape(T1,m,n);
    err=abs(norm((T-Told)));
    if err>=great,break; end
    error(IT)=err;
    if iflag
        fprintf('\n\nCompleted Outer Loop Iteration=%i\t Error=%2.6e\n\n',IT,err);
        if IT==1,err=1000;end
    end
    IT=IT+1;
end

toc
% Check Convergence & Plot Convergence History
IT=IT-1;
if(err>great)
    fprintf(1,'Iterations Diverged\n');
    fprintf(1,'Please Consider Change in Solution Parmeters and Run the Code Again\n');
    disp('Press any key')
    pause
    return;
end
if iflag
    figure('Name','Outer Loop Convergence History','NumberTitle','off','OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)/2]);
    fprintf(1,'Calculated by %i Iterations\n',IT);
    semilogy(1:IT,error(1:IT),'- r');
    xlabel('Iteration');
    ylabel('Error');
    title('Outer Loop Convergence History');
end


%% Post-Processing-----------------------------------------------------------

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
    q_t(i)=kn(i,n)*(T(i,n-1)-T(i,n))/abs(Y(i,n-1)-Y(i,n))*Dx(i,n);   %Flux @ the Bottom Boundary
    q_b(i)=ks(i,1)*(T(i,1)-T(i,2))/abs(Y(i,1)-Y(i,2))*Dx(i,2);   %Flux @ the Top Boundary
end
%-----Compute and Plot Fluxes @ East & South CV Faces
for j=1:n;
    for i=1:m-1;
        qx(i,j)=ke(i,j)*(T(i,j)-T(i+1,j))/abs(X(i+1,j)-X(i,j))*Dy(i,j);
    end
    qx(m,j)=-q_r(j);
end
for i=1:m;
    for j=1:n-1;
        qy(i,j)=ks(i,j)*(T(i,j)-T(i,j+1))/abs(Y(i,j+1)-Y(i,j))*Dx(i,j);
    end
    qy(i,n)=-q_b(i);    
end

%Plot Flux vectors
figure('Name','Flux Vectors','NumberTitle','off','OuterPosition',[1 1 scrsz(3)/2 scrsz(4)/2]);
contour(X,Y,T)
hold on
quiver(X,Y,qx,qy)
hold off

%Report Flux in inbalence @ Boundaries
fprintf('\nFlux inblalence @ Right & Left Boundaries is %2.2e\n',sum(q_r)+sum(q_l));
fprintf('Flux inblalence @ Top & Bottom Boundaries is %2.2e\n',sum(q_t)+sum(q_b));


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










