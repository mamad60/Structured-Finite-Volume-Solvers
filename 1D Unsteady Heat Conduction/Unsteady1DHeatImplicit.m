% This Scripts Solves 1D Unsteady Heat Conduction Equation
% rho*cp*DT/Dt=d/dx(kdT/dx)+S, From Patankar Book, Section 4-2-By Implicit
% Method
%Density(rho) and Thermal Capacity(cp) are constant
% Both Thermal Conductivity(k)and Source therm(S) can be varied with time & x
%Discritization Stencil   |--->  W---P---E
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


%Inputs
L=1;  % Lenght of Domain
m=100; % Number of Control Volumes
Ti=273; %Initial Value of T at t=0
rho=1; % Density
Cp=1; %Thermal Capacity
Dt=2; % Time Step
NumberIt=10;  % Number of time steps
epsT=1e-6;  %Convergence Criteria(If Steady Solution is Intended)
%-----------

kStat=0;  % Themal Conductivity Calculation Method
% KStat=0:Constant 1:Function of X 2:Function of Temperature 3: Function of Both
sStat=0;  %Source Term Calculation Method
%sStat=0: No Source Term 1:Constant 2:Function of X 3:Function of Temperature 4: Function of Both
iMethod=0; % Solution Method:
%0 MATLAB Direct Solver 1: TDMA(Direct) 2:Conjugate Gradient

%-----Boundary Conditions
NL=0; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Left Boundary
NR=0; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Right Boundary
Tl=0; % Temperature at The Left Boundary, Applies only if NL=0
Tr=0; % Temperature at The Right Boundary, Applies only if NR=0
Ql=1; % Flux at The Left Boundary, Applies only if NL=1
Qr=1; % Flux at The Right Boundary, Applies only if NR=1
%Flux Normal to boundary & Towards it Assumed Positive sign
hl=1; % Convection Coefficient at The Left Boundary, Applies only if NL=2
hr=1; % Flux at The Left Boundary, Applies only if NL=2
Tf=293.0; %Temperature @ infinity used if NL=2 Or NR=2 q=h(Tf-TB)

%---------------------Post-Processing Parmeters
animateT=1; %If set to 1 Frames for Temperatue Profile are animated every intAnimT Time Step
animateQ=1; %If set to 1 Frames for Temperatue Profile are animated every intAnimQ Time Step
intAnimT=5; %Temperature Results are animated every intanimT Time Steps
intAnimQ=5; %Heat Flux Results are animated every intanimT Time Steps--May Be Slow Due to Calcualtions
AVIOutputT=0; %1: Animation of  T Profile is Saved in AVI Format
AVIOutputQ=0; %1: Animation of  Q Profile is Saved in AVI Format
intSaveT=0;  %1: Temperature Profile are Saved each intanimT Time Steps(if animateT=1), File Format Tprofile+TimeStep---Very Slow
intSaveQ=0;  %1: Heat Flux Profile are Saved each intanimT Time Steps(if animateQ=1), File Format Qprofile+TimeStep----Very Slow
saveRestart=0; %1: Tempearature,Iteration number, L & m Data are saved to a file called Restart.mat to be used later
restart=0; %1: Tempearature, L & m Data are read are read from file Restart.mat and overwrited to privious values
%-----------------------------------------------------------

%Variable Allocation
X=zeros(m,1); %Location of Verices
Xvol=zeros(m-1,1); % Location of Control Volume Faces
T=zeros(m,1); %Variable for Storing Solution
Told=zeros(m,1); %Variable for Storing T @ Privious Time Step
A=zeros(m); % Coefficient Matrix
b=ones(m,1); % Right hand side Matrix
Dx=zeros(m,1); % Control Volume lenght for Source Terms Evaluation
ke=zeros(size(Xvol)); % k @ East Boundary of the CV
% kw=zeros(size(Xvol)); % k @ West Boundary of the CV
q=zeros(size(ke));
error=zeros(1,NumberIt);
err=1000; % Varible for Storing Error
great=1e10; % Assumsion for infinity
IT=1;   %Current Iteration Counter
resIt=1; %Initial Iteration number of restart
if animateT || animateQ
    mCountT=1; % Index for Movie Frames(T)
    mCountQ=1; % Index for Movie Frames(Q)
    % Preallocate movie structure.
    nFrames=floor(NumberIt/intAnimT)+1;
    movT(1:nFrames) = struct('cdata', [],'colormap', []);
    nFrames=floor(NumberIt/intAnimQ)+1;
    movQ(1:nFrames) = struct('cdata', [],'colormap', []);
end
%------------------------------------------------------------

%Grid Generation
[X,Xvol,Xw,Xe] = Grid1d(L,m);
% Calculate Dx of CVs
Dx(1)=Xvol(1);
for i=2:m-1;
    Dx(i)=abs(Xvol(i)-Xvol(i-1));
end
Dx(m)=X(m)-Xvol(m-1);
Ti=sin(pi*X/L);
%Show Time Step & Mesh Fourier Number
DX=min(Dx)*2;
alpha=max(ThermalConductivity(kStat,X,Ti))/(rho*Cp); %Fourier Number
F0m=(alpha*Dt)/(DX*DX); %Mesh Fourier Number
fprintf('\nCurrent Time Step is %2.3f\nMesh Fourier Number is=%2.3f\nEnjoy!!!\n',Dt,F0m);

%Set Initial Value of T Or read Restart File
if restart
    %Read Restart File
    load('Restart','T','L','m','IT')
    disp('Restart File loaded')
    resIt=NumberIt+1;
    NumberIt=IT+NumberIt-1;
else
    %Set Initial Values of T
    T(:)=Ti;
end

Texact=@(x,t) (sin(pi*x/L)*exp(-alpha*pi*pi*t/L^2)); %Validation

%Set Figures for Animation
scrsz = get(0,'ScreenSize');
if animateT
    hT=figure('Name','Animation of the Temperature Profile','NumberTitle','off','OuterPosition',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
end
if animateQ
    hQ=figure('Name','Animation of the Heat Flux Profile','NumberTitle','off','OuterPosition',[scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
end


%% Time Stepping
while (IT<=NumberIt) && (err>epsT)
    Told=T;
    %--------------------------------------------------------------
    %-------------------Calcalate k @ Boundies
    %The Left Boundary
    %--------Calculate ke
    kP=ThermalConductivity(kStat,X(1),Told(1));
    kE=ThermalConductivity(kStat,X(2),Told(2));
    Xep=abs(Xvol(1)-X(1));
    fe=Xep/Xe(1);
    ke(1)=1/((1-fe)/kP+fe/kE); %Compute ke by Biharmonic Averaging
    kl=ke(1);  %Store k for future calculations
    %The Right Boundary
    ae=0;  %only for clarification
    %--------Calculate ke
    kP=ThermalConductivity(kStat,X(m),Told(m));
    kW=ThermalConductivity(kStat,X(m-1),Told(m-1));
    Xwp=abs(X(m)-Xvol(m-1));
    fw=Xwp/Xw(m);
    kw=1/((1-fw)/kP+fw/kW);%Compute kw by Biharmonic Averaging
    kr=kw;  %Store k for future calculations
    %----------------------------------------------------------------
    %---------------Boundary Cells & Application of BCs
    %Left Boundary
    switch NL
        case 0  % Fixed Temperature
            b(1)=Tl;
            A(1,1)=1;
        case 1  %Fixed Flux
            aw=0;  %only for clarification
            %Calcalate k Source Terms Boundies
            [Sc,Sp] = SourceTerm(sStat,X(1),Told(1));
            %----------Set Coefficients
            Dxl=abs(X(2)-X(1));   %Distance First Two Points Adjacent to Left Boundary
            aI=-ke(1)/Dxl; % Point next to Boundary
            aB=aI-Sp*Dx(1);  % Boundary Point
            b(1)=Sc*Dx(1)+Ql;
            %Assemble Soloution Matrix
            A(1,1)=aB;
            A(1,2)=-aI;
        case 2  %Convection Coeffient is Given
            aw=0;  %only for clarification
            %Calcalate k Source Terms Boundies
            [Sc,Sp] = SourceTerm(sStat,X(1),Told(1));
            %----------Set Coefficients
            Dxl=abs(X(2)-X(1));   %Distance First Two Points Adjacent to Left Boundary
            aI=ke(1)/Dxl; % Point next to Boundary
            aB=aI-Sp*Dx(1)+hl;  % Boundary Point
            b(1)=Sc*Dx(1)+hl*Tf;
            %Assemble Soloution Matrix
            A(1,1)=aB;
            A(1,2)=-aI;
    end
    %Right Boundary
    switch NR
        case 0  % Fixed Temperature
            b(m)=Tr;
            A(m,m)=1;
        case 1  %Fixed Flux
            ae=0;  %only for clarification
            %Calcalate k Source Terms Boundies
            [Sc,Sp] = SourceTerm(sStat,X(m),Told(m));
            %----------Set Coefficients
            Dxr=abs(X(m)-X(m-1)); %Distance First Two Points Adjacent to Right Boundary
            aI=kw/Dxr; % Point next to Boundary
            aB=aI-Sp*Dx(m);  % Boundary Point
            b(m)=Sc*Dx(m)+Qr;
            %Assemble Soloution Matrix
            A(m,m)=aB;
            A(m,m-1)=-aI;
        case 2  %Convection Coeffient is Given
            Dxr=abs(X(m)-X(m-1)); %Distance First Two Points Adjacent to Right Boundary
            ae=0;  %only for clarification
            [Sc,Sp] = SourceTerm(sStat,X(m),Told(m));
            %----------Set Coefficients
            aI=kw/Dxr; % Point next to Boundary
            aB=aI-Sp*Dx(m)+hr;  % Boundary Point
            b(m)=Sc*Dx(m)+hr*Tf;
            %Assemble Soloution Matrix
            A(m,m)=aB;
            A(m,m-1)=-aI;
    end
   
    %------------Interior Cells
    for i=2:m-1
        %Compute k @ CV Faces
        kP=ThermalConductivity(kStat,X(i),Told(i));
        kW=ThermalConductivity(kStat,X(i-1),Told(i-1));
        kE=ThermalConductivity(kStat,X(i+1),Told(i+1));
        %------Compute by Biharmonic Averaging
        %------------------------ke
        Xep=abs(Xvol(i)-X(i));
        fe=Xep/Xe(i);
        ke(i)=1/((1-fe)/kP+fe/kE);
        %------------------------kw
        Xwp=abs(X(i)-Xvol(i-1));
        fw=Xwp/Xw(i);
        kw=1/((1-fw)/kP+fw/kW);
        %Compute Source Terms
        [Sc,Sp] = SourceTerm(sStat,X(i),Told(i));
        %Compute Coefficients
        %ae
        ae=ke(i)/Xe(i);
        %aw
        aw=kw/Xw(i);
        %ap0-----------------Unsteady Coe. Applies The effect of Previous Time Step
        ap0=(rho*Cp*Dx(i))/Dt;
        %ap
        ap=ae+aw+ap0-Sp*Dx(i);
        %b
        b(i)=Sc*Dx(i)+ap0*Told(i);
        %Assemble Soloution Matrix
        A(i,i)=ap;
        A(i,i-1)=-aw;
        A(i,i+1)=-ae;  
    end
    %Solve Equations
    A=sparse(A); %Improves solution time for Matrix Solvers
    switch iMethod
        case 0
            T=A\b;
        case 1
            T=TDMA(A,b);
        case 2
            T=CG(A,T,b,10000,1e-6);
        otherwise
            display('Please First Select a solver')
            pause
            return;
    end
    %Animation
    %----Plot  T Profile
    if animateT
        if mod(IT,intAnimT)==0 || IT==1
            figure(hT);
            plot(X,T,'LineWidth',2);
            xlabel('X Coordinate');
            ylabel('T');
            xlim([0 L]);
            title(strcat('Profile of Temperature  ','  (','Time Step=',num2str(IT),' ,Time= ',num2str(Dt*(IT)),')'))
            drawnow
            pause(0.5);
            movT(mCountT)=getframe(gcf);
            mCountT=mCountT+1;
            if intSaveT, print(hT,'-dtiff',strcat('Tprofile',num2str(IT))) ,end
        end
    end
    %----Plot  Q Profile
    if animateQ 
        if mod(IT,intAnimQ)==0 || IT==1
            %-----Compute  Fluxes @ East CV Faces
            for i=1:m-1
                q(i)=ke(i)*(T(i+1)-T(i))/abs(X(i+1)-X(i));
            end
            figure(hQ);
            plot(Xvol,q,'g','LineWidth',2);
            xlabel('X Coordinate of CV Faces');
            h2=ylabel('q_{e}=-q_{w}');
            set(h2, 'interpreter', 'tex');
            if max(q)~=min(q)
                ylim([min(q) max(q)]);
            end
            title(strcat('Profile of Heat Flux  ','  (','Time Step=',num2str(IT-1),' ,Time= ',num2str(Dt*(IT-1)),')'));
            drawnow
            pause(0.5);
            movQ(mCountQ)=getframe(gcf);
            mCountQ=mCountQ+1;
            if intSaveQ, print(hQ,'-dtiff',strcat('Qprofile',num2str(IT))) ,end
        end
    end
    %Compute Steady Errors
    err=abs(norm((T-Told)));
    fprintf('IT=%i\t Steady Error=%2.6e  Time=%2.3f\n',IT,err,IT*Dt);
    if err>=great
        break;
    end
    IT=IT+1;
    error(IT-resIt)=err;
end
% Check Convergence & Plot Convergence History
if(error(IT-resIt)>1000)
    fprintf(1,'Iterations Diverged\n');
    fprintf(1,'Please Consider Change in Time Step Or Solution Parmeters and Run the Code Again\n');
    disp('Press any key')
    pause
    return;
end
figure
plot(resIt:NumberIt,log(error),'- g','LineWidth',1.5);
xlabel('Iteration');
ylabel('Log(Steady Error)');
xlim([resIt IT])
title('Convergence History');

%% Post-Processing-----------------------------------------------------------

%---Report Fluxes @ Right & Left Boundries
ql=kl*(T(2)-T(1))/abs(X(1)-X(2));   %Flux @ the Left Boundary
qr=kr*(T(m)-T(m-1))/abs(X(m)-X(m-1));   %Flux @ the Left Boundary
fprintf('\nThe Flux @ the Left Boundary is: %2.4e\n',ql);
fprintf('The Flux @ the Right Boundary is: %2.4e\n',qr);
%---Report Temeratures @ Right & Left Boundries
fprintf('The Temperature @ the Left Boundary is: %2.4e\n',T(1));
fprintf('The Temperature @ the Right Boundary is: %2.4e\n',T(m));

%Save Restart File
if saveRestart
    save('Restart','T','L','m','IT')
    disp('Restart File Saved')
end
%Plot results
%----Plot  T Profile
figure('Name','Temperature Profile at the Final Time Step','NumberTitle','off','OuterPosition',[1 1 scrsz(3)/2 scrsz(4)/2]);plot(X,T,'LineWidth',2);
xlabel('X Coordinate');
ylabel('T');
xlim([0 L]);
title(strcat('Profile of Temperature @ Final Time','  (','Time Step=',num2str(IT-1),' ,Time= ',num2str(Dt*(IT-1)),')'));
%-----Compute and Plot Fluxes @ East CV Faces

for i=1:m-1
    q(i)=ke(i)*(T(i+1)-T(i))/abs(X(i+1)-X(i));
end
%------------------Plot Fluxes
figure('Name','Heat Flux Profile at the Final Time Step','NumberTitle','off','OuterPosition',[scrsz(3)/2 1 scrsz(3)/2 scrsz(4)/2]);plot(X,T,'LineWidth',2);
plot(Xvol,q,'g','LineWidth',2);
xlabel('X Coordinate of CV Faces');
h2=ylabel('q_{e}=-q_{w}');
set(h2, 'interpreter', 'tex');
if max(q)~=min(q)
    ylim([min(q) max(q)]);
end
title(strcat('Profile of Heat Flux @ Final Time ','  (','Time Step=',num2str(IT-1),' ,Time= ',num2str(Dt*(IT-1)),')'));

%---if k is a function of x Plot its Variation @ CV Faces
if  kStat==1 || kStat==2 || kStat==3
    figure
    plot(Xvol,ke,'m','LineWidth',2);
    h1=xlabel('X Coordinate of CV Faces');
    ylabel('k_{e}');
    set(h1, 'interpreter', 'tex');
    title('Thermal Conductiviy at CV Faces(East) @ Final Time')
end

%--------------Play Movie
if AVIOutputT && animateT
    movie2avi(movT, 'TemperatureProfile','fps',5,'quality', 100,'compression','None')
    display('AVI Movie for Temperature Profile exported to the Current Directory')
end
if AVIOutputQ && animateQ
    movie2avi(movQ, 'HeatFluxProfile','fps',5,'quality', 100,'compression','None')
    display('AVI Movie for Heat Flux Profile exported to the Current Directory')
end

%Compare Exact & Numerical Solution
figure('Name','Comparison Between Exact & Numerical Solution','NumberTitle','off');
plot(X,T,X,Texact(X,(IT-1)*Dt),'or')
title(strcat('Comparison Between Exact & Numerical Solution  ','  (','Time Step=',num2str(IT-1),' ,Time= ',num2str(Dt*(IT-1)),')'));
legend('Numerical','Exact',2);


fprintf('\nGood Lock, Mohammad Aghakhani,2010\n');











