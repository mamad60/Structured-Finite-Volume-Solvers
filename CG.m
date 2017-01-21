function [x] = CG(A,X,b,MaxIT,tol )
%Conjugate Gradient Method For Solotion of a linear system Ax=b-Matrix Form
%By Mohammad Aghakhani, Feel Free using for Educaitonal & Research Purposes
r=b-A*X;
p=r;
r2=dot(r,r);
IT=1;
display('Please wait..... Calculating the solution')
while(sqrt(r2)>=tol) && (IT<=MaxIT) 
    Ap=A*p;
    alpha=r2/dot(Ap,p);
    %alpha=r'*r/Ap'*p;%------Removed--Slower in MATLAB
    X=X+alpha*p;
    rold=r;
    r=r-alpha*Ap;
    r2=dot(r,r);
    betha=r2/dot(rold,rold);
    p=r+betha*p;
    %fprintf('Conjugate Gradien Iteration=%i\tResidual=%2.6e\n',IT,res);
    IT=IT+1;
end
if IT>MaxIT && r2>tol
    fprintf('\nMaximum Iteratons Reached. Gauss-Seidel Solver Diverged\n');
    fprintf('Conjugate Gradient Iteration=%i\tResidual=%2.6f\n',IT-1,sqrt(r2));
    x=X;
else
    fprintf('\nConjugate Gradient Solver Converged');
    fprintf('\nConjugate Gradient Iteration=%i\tResidual=%2.6e\n',IT-1,sqrt(r2));
end
x=X;
end
