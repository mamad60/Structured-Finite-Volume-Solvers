function [v_h] = MV(A,v,f,n_level,nmax)
%V-cycle (recursive form),For solution of Au=f--> ,Ae=r is a form of it
%n_level is number Current levels-1: nmax:Finest Grid 1:Coarset Grid
%Initial Guess is stored under v

n=length(f);
I2h_h=Inject(n); % h to 2h space restriction by injection(Fine to coarse)
Ih_2h=ProlongI(floor(n/2)); % 2h to h prolongation (Coarse to fine)
alpha1=2; % Number of Initial Relaxations on the fine grid
alpha2=5;  %Number of Final Relaxations on the fine grid
%Pre smoothing: Apply the smoother Alpha1 times to A_h*u_h = f_h with the initial guess v_h
v=CG1(A,v,f,alpha1,1e-6);
%If h is the coarsest grid
%? solve the problem
if n_level~=1  %----If we are in the coarsest grid-Just solve for erros
    %Restrict to the next coarser grid: f_2h = Ih_2h(f_h = A_h*v_h)
    %-----------Compute f_2h and  A_2h
    %-----------Compute Residual on the fine grid
    f_h=f-A*v;
    f_2h=I2h_h*f_h;
    %---------------Compute A_2h----->Coarse Grid Operator
    if(~mod(n,2))
        A_2h=I2h_h*A*Ih_2h;
    else
        A_2h=I2h_h(1:end,1:end-1)*A(1:end-1,1:end-1)*Ih_2h;
    end
    % Set initial iterate on the next coarser grid: v_2h = 0
    % If h is the finest grid, set lambda = 1.
    %lambda=2; %set 1: V-Cycle 2:W-Cycle
    % if n==nmax, lambda=1;end
    % Call the V-cycle scheme lambda times for the next coarser grid:
    % v_2h=MV_ 2h(v_2h, f_2h)
    %for i=1:lambda
    v_2h=MV(A_2h,zeros(size(f_2h)),f_2h,n_level-1,nmax);
    %end
    %3 Correct with the prolongated update: v_h = v_h + Ih_2h*v_2h
    if(~mod(n,2))
        v=v+Ih_2h*v_2h;
    else
        v(1:end-1)=v(1:end-1)+Ih_2h*v_2h;
    end
end
%4. Post smoothing: Apply the smoother ?2 times to A_h*u_h = f h with the initial guess v_h
if n_level==1
    v=A\f;
else
    v=CG1(A,v,f,alpha2,1e-6);
end
v_h=v;
end

