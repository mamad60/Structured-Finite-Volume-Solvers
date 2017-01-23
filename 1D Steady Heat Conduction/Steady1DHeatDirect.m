% This Scripts Solves 1D Steady Heat Transfer Equation
% d/dx(kdT/dx)+S=0, From Patankar Book, Section 4-2
% Both Thermal Conductivity and Source therms can be varied with time & x
%Discritization Stencil   |--->  W---P---E
%By Mohammad Aghakhani, 2010,for Teaching purposes on the CFD Class

clc
clear
close all

%Inputs
L=1;  % Lenght of Domain
m=100; % Number of Control Volumes
kStat=0;  % Themal Conductivity Calculation Method
% KStat=0:Constant 1:Function of X 2:Function of Temperature 3: Function of Both
sStat=2;  %Source Term Calculation Method
%0: No Source Term 1:Constant 2:Function of X 3:Function of Temperature 4: Function of Both
sMethod=1; % Solution Method: 0 MATLAB Direct Solver 1: TDMA(Thomas)
%-----Boundary Conditions
NL=0; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Left Boundary
NR=0; % 0:Fixed Temp. 1:Flux Bc 2:Convection on the Right Boundary
Tl=0; % Temperature at The Left Boundary, Applies only if NL=0
Tr=0; % Temperature at The Right Boundary, Applies only if NR=0
Ql=0.1; % Flux at The Left Boundary, Applies only if NL=1
Qr=1; % Flux at The Right Boundary, Applies only if NR=1
%Flux Normal to boundary & Towards it Assumed Positive sign
hl=1; % Convection Coefficient at The Left Boundary, Applies only if NL=2
hr=0.1; % Flux at The Left Boundary, Applies only if NL=2
Tf=293.0; %Temperature @ infinity used if NL=2 Or NR=2 q=h(Tf-TB)

%------------Only if Thermal Conductivity Or Source Term are function of
%Temperature ,variables Defined in this Section are used
IT=1;   %Current Iteration Counter
MaxIT=10000;  %Maximum Allowed Iteration
epsT=1e-6;   %Convergence Criteria in Two Concecutive Iterations
%-----------------------------------------------------------
%Variable Allocation
X=zeros(m,1); %Location of Verices
Xvol=zeros(m-1,1); % Location of Control Volume Faces
T=zeros(m,1); %Variable for Storing Solution
Told=zeros(m,1); %Variable for Storing Old Iteration Solution
A=zeros(m); % Coefficient Matrix
b=ones(m,1); % Right hand side Matrix
Dx=zeros(m,1); % Control Volume lenght for Source Terms Evaluation
ke=zeros(size(Xvol)); % Used For Flux Calculation
err=1000; % Varible for Storing Error for Outer loop
great=1e10; % Assumsion for infinity
%------------------------------------------------------------

%Grid Generation
[X,Xvol,Xw,Xe] = Grid1d(L,m);
% Calculate Dx of CVs
Dx(1)=Xvol(1);
for i=2:m-1;
    Dx(i)=abs(Xvol(i)-Xvol(i-1));
end
Dx(m)=X(m)-Xvol(m-1);
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
T(:)=min(Tl,Tr);
if iflag
    error=zeros(1,MaxIT);
end

%% Outer Loop : Perfomed more than one if k(t) or s(x)
tic;
while (IT<=MaxIT) && (err>epsT)
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
    
    %Calculate Coefficients & Assemble Solution Matrix AT=b
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
            aI=ke(1)/Dxl; % Point next to Boundary
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
            ae=0;  %only for clarification
            [Sc,Sp] = SourceTerm(sStat,X(m),Told(m));
            %----------Set Coefficients
            Dxr=abs(X(m)-X(m-1)); %Distance First Two Points Adjacent to Right Boundary
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
        %ap
        ap=ae+aw-Sp*Dx(i);
        %b
        b(i)=Sc*Dx(i);
        %Assemble Soloution Matrix
        A(i,i)=ap;
        A(i,i-1)=-aw;
        A(i,i+1)=-ae;
    end
    %Solve Equations
    switch sMethod
        case 0
            %         t1=tic;
            T=A\b;
            %             display('Direct Matlab Solver')
            %         toc(t1)
        case 1
            %         t2=tic;
            T=TDMA(A,b);
            %             display('TDMA Method')
            %         toc(t2)
    end
%     T=uRelax*T-(1-uRelax)*Told;
    if iflag
        err=abs(norm((T-Told)));
        fprintf('IT=%i\t Error=%2.6e\n',IT,err);
        if err>=great
            break;
        end
        error(IT)=err;
    end
    IT=IT+1;
end
toc;
% Check Convergence & Plot Convergence History
if iflag
    if(error(IT-1)>1000)
        fprintf(1,'Iterations Diverged\n');
        fprintf(1,'Please Consider Change in Solution Parmeters and Run the Code Again\n');
        disp('Press any key')
        pause
        return;
    end
    if(error(IT)<eps)
        fprintf(1,'Converged in %i Iterations\n',IT);
        plot(1:IT,log(error(1:IT)),'- r');
        xlabel('Iteration');
        ylabel('Log(Error)');
        title('Convergence History');
    else
        disp('Maximum Iteration Number Reached');
        plot(1:IT,lg(error(1:IT)),'- r');
        xlabel('Iteration');
        ylabel('Log(Error)');
        title('Convergence History');
        pause;
        return;
    end
end


%% Post-Processing-----------------------------------------------------------

%---Report Fluxes @ Right & Left Boundries
ql=kl*(T(1)-T(2))/(X(2)-X(1));   %Flux @ the Left Boundary
qr=kr*(T(m)-T(m-1))/(X(m)-X(m-1));   %Flux @ the Left Boundary
fprintf('The Flux @ the Left Boundary is: %2.4e\n',ql);
fprintf('The Flux @ the Right Boundary is: %2.4e\n',qr);
%---Report Temeratures @ Right & Left Boundries
fprintf('The Temperature @ the Left Boundary is: %2.4e\n',T(1));
fprintf('The Temperature @ the Right Boundary is: %2.4e\n',T(m));


%Plot results
%----Plot  T Profile
figure
plot(X,T,'LineWidth',2);
xlabel('X Coordinate');
ylabel('T');
xlim([0 L]);
if T(m)>T(1)
    ylim([T(1) T(m)]);
end
title('Profile of T');
%-----Compute and Plot Fluxes @ East CV Faces

for i=1:m-1
    q(i)=ke(i)*(T(i+1)-T(i))/abs(X(i+1)-X(i));
end
%------------------Plot Fluxes
figure
plot(Xvol,q,'g','LineWidth',2);
xlabel('X Coordinate of CV Faces');
h2=ylabel('q_{e}=-q_{w}');
set(h2, 'interpreter', 'tex');
if max(q)~=min(q)
    ylim([min(q) max(q)]);
end
title('Profile of Heat Flux');

%---if k is a function of x Plot its Variation @ CV Faces
if kStat==1 || kStat==3
    figure
    plot(Xvol,ke,'m','LineWidth',2);
    h1=xlabel('X Coordinate of CV Faces');
    ylabel('k_{e}');
    set(h1, 'interpreter', 'tex');
    title('Thermal Counductiviy at CV Faces(East)')
end


fprintf('Good Lock, Mohammad Aghakhani,2010\n');











