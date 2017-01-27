% This Scripts Solves 2D UnSteady Heat Transfer Equation-Implicit Time Stepping
% rho*cp*DT/Dt=d/dx(kdT/dx)+d/dy(kdT/dy)+S=0, From Patankar Book, Section 4-4
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
m=20; % Number of Control Volumes in X direction(i)
n=20; % Number of Control Volumes in Y direction(j)
rho=1; % Density
Cp=1; %Thermal Capacity
Dt=10; % Time Step
NumberIt=10;  % Number of time steps
epsT=1e-6;  %Convergence Criteria(If Steady Solution is Intended)
Ti=273; %Initial Value of T at t=0
%------------
kStat=0;  % Themal Conductivity Calculation Method
% KStat=0:Constant 1:Function of X(and/Or)Y 2:Function of Temperature 3: Function of Both Space & Temp.
sStat=0;  %Source Term Calculation Method
%0: No Source Term 1:Constant 2:Function of X(and/Or)Y 3:Function of Temperature 4: Function of Both
iMethod=4; % Solution Method:
%0 MATLAB Direct Solver
%2:Jacobi(Matrix) 3:Jacobi 4:Gauss-Seidel(Matrix) 5:Gauss-Seidel
%6: Symmetric Gauss-Seidel Matrix 7:Symmetric Gauss-Seidel
%8:SOR(Matrix) 9:SOR 10:Symmetric SOR Matrix
%11:Red-Black Gauss-Seidel  12:SOR Solver 4 Triangular matrix
%13:Conjugate Gradient   14: Preconditioned Conjugate Gradient
%15: BiConjugate gradient
espSolver=1e-6; % Maximum Residual of Matrix Iterative Solver
MaxITSolver=100000; % Maximum Iteraton of Matrix Iterative Solver
%-----------------------------------------------------------

%---------------------Post-Processing Parmeters
animateT=1; %If set to 1 Frames for Temperatue Countour are animated every intAnimT Time Step
animateQ=1; %If set to 1 Frames for Temperatue Countor are animated every intAnimQ Time Step
animateSurfT=1; %If set to 1 Frames for Temperatue Surface are animated every intAnimT Time Step
intAnimT=5; %Temperature Contour Results are animated every intanimT Time Steps
intAnimQ=5; %Heat Flux Results are animated every intanimT Time Steps--May Be Slow Due to Calcualtions
intAnimST=5; %Temperature Surface are animated every intanimT Time Steps
AVIOutputT=0; %1: Animation of  T Contour is Saved in AVI Format
AVIOutputQ=0; %1: Animation of  Q Contour is Saved in AVI Format
AVIOutputST=0; %1: Animation of  T Surface is Saved in AVI Format
intSaveT=0;  %1: Temperature Contour  are Saved each intanimT Time Steps(if animateT=1), File Format Tprofile+TimeStep--- Slow
intSaveQ=0;  %1: Heat Flux Contour & Surface are Saved each intanimT Time Steps(if animateQ=1), File Format Qprofile+TimeStep----Very Slow
intSaveST=0; %1: Temperature Surface are Saved each intanimST Time Steps(if animateST=1), File Format Tprofile+TimeStep---Slow
saveRestart=0; %1: Tempearature,Iteration number, L & m Data are saved to a file called Restart.mat to be used later
restart=0; %1: Tempearature, L & m Data are read are read from file Restart.mat and overwrited to privious values
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
error=zeros(NumberIt);
IT=1;   %Current Iteration Counter
resIt=1; %Initial Iteration number of restart
%------------------------------------------------------------

%Grid Generation & Geometrical Calculations
Grid2d(Lx,Ly,m,n);
%Show Time Step & Mesh Fourier Number
DX=min(min(Dx))*2;
DY=min(min(Dy))*2;
alpha=max(max(ThermalConductivity(kStat,X,Y,Ti))/(rho*Cp)); %Fourier Number
F0m=(alpha*Dt)/(DX*DY); %Mesh Fourier Number
fprintf('\nCurrent Time Step is %2.3f\nMesh Fourier Number is=%2.3e\nEnjoy!!!\n',Dt,F0m);
%Set Initial Value of T Or read Restart File
if restart
    %Read Restart File
    load('Restart','T','Lx','Ly','m','n','IT')
    disp('Restart File loaded')
    resIt=NumberIt+1;
    NumberIt=IT+NumberIt-1;
else
    %Set Initial Values
    T(:)=Ti;
end

%Allocate Matrix for storing Errors for Post-Processing and Initiale T
T(:)=273;
%Set Boundary Conditions
SetBC
%Set Figures for Animation
scrsz = get(0,'ScreenSize');
if animateT
    mCountT=1; % Index for Movie Frames(T)
    nFrames=floor(NumberIt/intAnimT)+1;
    % Preallocate movie structure.
    movT(1:nFrames) = struct('cdata', [],'colormap', []);
    hT=figure('Name','Animation of the Temperature Contours','NumberTitle','off','OuterPosition',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
end
if animateSurfT
    mCountST=1; % Index for Movie Frames(T)
    nFrames=floor(NumberIt/intAnimST)+1;
    % Preallocate movie structure.
    movST(1:nFrames) = struct('cdata', [],'colormap', []);
    hST=figure('Name','Animation of the Temperature Surface','NumberTitle','off','OuterPosition',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
end
if animateQ
    mCountQ=1; % Index for Movie Frames(Q)
    % Preallocate movie structure.
    nFrames=floor(NumberIt/intAnimQ)+1;
    movQ(1:nFrames) = struct('cdata', [],'colormap', []);
    hQ=figure('Name','Animation of the Heat Flux Vectors','NumberTitle','off','OuterPosition',[1 1 scrsz(3)/2 scrsz(4)/2]);
end


%% Time Stepping
tic;
while (IT<=NumberIt) && (err>epsT)
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
            %ap0-----------------Unsteady Coe. Applies The effect of Previous Time Step
            ap0=(rho*Cp*Dx(i)*Dy(i))/Dt;
            %ap
            ap=an+as+ae+aw+ap0-Sp(i,j)*Dx(i,j)*Dy(i,j);
            % Set Labaling Coe.
            np=(j-1)*m+i;
            %b
            b(np)=Sc(i,j)*Dx(i,j)*Dy(i,j)+ap0*Told(i,j);
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
    %Compute Steady Errors
    err=abs(norm((T-Told)));
    fprintf('IT=%i\t Steady Error=%2.6e  Time=%2.3f\n',IT,err,IT*Dt);
    %Animation
    %     if animateT || animateQ, Animation(IT,hT,hQ,Dt), end
    %---------------------------------------------------------------------
    %Animation
    %----Plot  T Countours
    if animateT
        if mod(IT,intAnimT)==0 || IT==1
            figure(hT);
            %----Plot  T Contour
            contourf(X,Y,T,'LineWidth',2);
            xlabel('X Coordinate');
            ylabel('Y Coordinate');
            title(strcat('Temperature Contour  ','  (','Time Step=',num2str(IT),' ,Time= ',num2str(Dt*(IT)),')'));
            drawnow
            pause(1);
            movT(mCountT)=getframe(gcf);
            mCountT=mCountT+1;
            if intSaveT, print(hT,'-dtiff',strcat('Tconntour',num2str(IT))) ,end
        end
    end
    if animateSurfT
        if mod(IT,intAnimST)==0 || IT==1
            figure(hST);
            surf(X,Y,T);
            xlabel('X Coordinate');
            ylabel('Y Coordinate');
            title(strcat('Surface of Temperature @ Final time ','  (','Time Step=',num2str(IT),' ,Time= ',num2str(Dt*(IT)),')'));
            drawnow
            pause(1);
            movST(mCountST)=getframe;
            mCountST=mCountST+1;
            if intSaveST, print(hST,'-dtiff',strcat('TSurface',num2str(IT))) ,end
        end
    end
    %----Plot  Q Profile
    if animateQ % Use with precuation, Slow down the code very much
        if mod(IT,intAnimQ)==0 || IT==1
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
            figure(hQ);
            %Plot Flux vectors
            contour(X,Y,T)
            hold on
            quiver(X,Y,qx,qy)
            colormap hsv
            hold off
            title(strcat('Profile of Heat Flux  ','  (','Time Step=',num2str(IT),' ,Time= ',num2str(Dt*(IT)),')'));
            drawnow
            pause(1);
            movQ(mCountQ)=getframe(gcf);
            mCountQ=mCountQ+1;
            if intSaveQ, print(hQ,'-dtiff',strcat('Qprofile',num2str(IT))) ,end
        end
    end
    %End of Animation--------------------------------------------------------------------
    if err>=great,break; end
    IT=IT+1;
    error(IT-resIt)=err;
end

toc
% Check Convergence & Plot Convergence History
if(err>great)
    fprintf(1,'Iterations Diverged\n');
    fprintf(1,'Please Consider Change in Solution Parmeters and Run the Code Again\n');
    disp('Press any key')
    pause
    return;
end
figure('Name','Steady Error Convergence History','NumberTitle','off','OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)/2]);
semilogy(resIt:NumberIt,error,'- g','LineWidth',1.5);
xlabel('Iteration');
ylabel('Steady Error');
xlim([resIt IT-1])
title('Convergence History');


%% Post-Processing-----------------------------------------------------------

%---Report Average Temeratures @ Right , Left, Top & Bottom Boundaries
fprintf('\nThe Average Temperature @ the Left Boundary at Final time is: %2.2e\n',mean(T(1,:)));
fprintf('The Average Temperature @ the Right Boundary is: at Final time %2.2e\n',mean(T(m,:)));
fprintf('The Average Temperature @ the Top Boundary is at Final time: %2.2e\n',mean(T(:,n)));
fprintf('The Average Temperature @ the Bottom Boundary at Final time is: %2.2e\n',mean(T(1,:)));


%Save Restart File
if saveRestart
    save('Restart','T','Lx','Ly','m','n','IT')
    disp('Restart File Saved')
end

%Plot results
%----Plot  T Contour
figure('Name','Temperature Contour at Final time','NumberTitle','off','OuterPosition',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
contourf(X,Y,T,'LineWidth',2);
xlabel('X Coordinate');
ylabel('Y Coordinate');
title('');
title(strcat('Temperature Contour @ Final Time  ','  (','Time Step=',num2str(IT-1),' ,Time= ',num2str(Dt*(IT-1)),')'));

colorbar
%----Surface of Temperature
figure('Name','Surface of Temperature at Final time','NumberTitle','off','OuterPosition',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
surf(X,Y,T);
xlabel('X Coordinate');
ylabel('Y Coordinate');
title(strcat('Surface of Temperature @ Final time ','  (','Time Step=',num2str(IT-1),' ,Time= ',num2str(Dt*(IT-1)),')'));

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
figure('Name','Flux Vectors at Final time','NumberTitle','off','OuterPosition',[1 1 scrsz(3)/2 scrsz(4)/2]);
contour(X,Y,T)
hold on
quiver(X,Y,qx,qy)
hold off
xlabel('X Coordinate');
ylabel('Y Coordinate');
title(strcat('Flux Vectors @ Final time ','  (','Time Step=',num2str(IT-1),' ,Time= ',num2str(Dt*(IT-1)),')'));



%Report Flux in inbalence @ Boundaries
fprintf('\nFlux inblalence @ Right & Left Boundaries at Final time is %2.2e\n',sum(q_r)+sum(q_l));
fprintf('Flux inblalence @ Top & Bottom Boundaries at Final time is %2.2e\n\n',sum(q_t)+sum(q_b));


%---if k is a function of x Plot its Variation @ CV Faces
if kStat==1 || kStat==2 || kStat==3
    figure('Name','Thermal Conductivity Contour at Final time','NumberTitle','off');
    contourf(X,Y,ThermalConductivity(kStat,X,Y,T),'LineWidth',2);
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    title('Thermal Conductivity @ Final time')
    colorbar
end

%--------------Play Movie
if AVIOutputT && animateT
    movie2avi(movT, 'TemperatureContour','fps',2,'quality', 100,'compression','None')
    display('AVI Movie for Temperature Contour exported to the Current Directory')
end
if AVIOutputQ && animateQ
    movie2avi(movQ, 'HeatFluxVector','fps',2,'quality', 100,'compression','None')
    display('AVI Movie for Heat Flux Vector exported to the Current Directory')
end
if AVIOutputST && animateSurfT
    movie2avi(movST, 'Temerature Surface Plot','fps',2,'quality', 100,'compression','None')
    display('AVI Movie for Temerature Surface Plot exported to the Current Directory')
end


fprintf('\nGood Lock, Mohammad Aghakhani,2012\n');










