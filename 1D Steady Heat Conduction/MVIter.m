function [x] = MVIter(A,X,b,MaxIT,tol,nmax)
% Scheme For Solotion of a linear system Ax=b
% MaxIt=Maxium number of Iteration tol=Maximun allowable residual
%X must contain initial guess for X, nmax= Total number of Grid Levels
%Conjugate Gradiant Method is used for relaxing on the fine grid-Two Grid
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
%Construct Restriction & Prolongation Matrices

r=b-A*X;
res=norm(r); %Compute Residual
IT=1;
display('Please wait..... Calculating the solution')
%Compute Restriction and Prolongation Marices
while(res>=tol) && (IT<=MaxIT)
    %1-Call V-Cycel of the Finest(Top) level
    X = MV(A,X,b,nmax,nmax);
    %2-Compute Residual on the fine grid
    res=norm(b-A*X); %Compute Residual
    fprintf('V Cycle MG Iteration=%i\tResidual=%2.6e\n',IT,res);
    IT=IT+1;
end
if IT>MaxIT && res>tol
    fprintf('\nMaximum Iteratons Reached. V Cycle MG Solver Diverged\n');
    fprintf('V Cycle MG Iteration=%i\tResidual=%2.6e\n',IT-1,res);
    x=X;
else
    fprintf('\nV Cycle MG Solver Converged');
    fprintf('\nV Cycle MG Iteration=%i\tResidual=%2.6e\n',IT-1,res);
end
x=X;
end

