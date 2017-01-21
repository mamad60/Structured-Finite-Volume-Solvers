function [x] = CoarseGridCorrection(A,X,b,MaxIT,tol)
%Coarse Grid Correction Scheme For Solotion of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X
%Conjugate Gradiant Method is used for relaxing on the fine grid-Two Grid
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
r=b-A*X;
res=norm(r); %Compute Residual
IT=1;
display('Please wait..... Calculating the solution')
alpha1=1;  %Number of Initial Relaxations on the fine grid
alpha2=1;  %Number of Final Relaxations on the fine grid
%Compute Restriction and Prolongation Marices
n=length(b);
I2h_h=Inject(n); % h to 2h space restriction by injection(Fine to coarse)
Ih_2h=ProlongI(floor(n/2)); % 2h to h prolongation (Coarse to fine)
while(res>=tol) && (IT<=MaxIT)
    %1-Relax alpha1 times on the fine grid with initial guess of X
    X=CG1(A,X,b,alpha1,tol);
    %2-Compute Residual on the fine grid
    r_h=b-A*X;
    %3-Compute Residual on the Coarse grid by restriction I(h,2h)-Inject
    r_2h=I2h_h*r_h;
    %4-Solve for Error in coarse grid
    %---------------Compute A_2h----->Coarse Grid Operator
    if(~mod(n,2))
        A_2h=I2h_h*A*Ih_2h;
    else
        A_2h=I2h_h(1:end,1:end-1)*A(1:end-1,1:end-1)*Ih_2h;
    end
    %-------------- Solve Rsidual Equaiton by Dirct method
    e_2h=A_2h\r_2h;
    %5-Correct fine-grid solution
    if(~mod(n,2))
        X=X+Ih_2h*e_2h;
    else
        X(1:end-1)=X(1:end-1)+Ih_2h*e_2h;
    end
    %6-Relax alpha2 times on the fine grid with initial guess of X_h-Final
    X=CG1(A,X,b,alpha2,tol);
    res=norm(b-A*X); %Compute Residual
    fprintf('MGC Iteration=%i\tResidual=%2.6e\n',IT,res);
    IT=IT+1;
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. MGC Solver Diverged\n');
    fprintf('MGC Iteration=%i\tResidual=%2.6e\n',IT-1,res);
    x=X;
else
    fprintf('\nMGC Solver Converged');
    fprintf('\nMGC Iteration=%i\tResidual=%2.6e\n',IT-1,res);
end
x=X;
end

