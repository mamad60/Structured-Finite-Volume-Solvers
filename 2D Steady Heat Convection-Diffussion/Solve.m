function [x_out] = Solve(A,x,b,MaxITSolver,espSolver,iMethod)
%Solves the Equation AX=b by the Method given in iMethod, 
%----------Ax=b--------

    switch iMethod
        case 0
            x_out=A\b;
        case 1
            x_out=TDMA(A,b);
        case 2
            x_out=JacobiMatrix(A,x,b,MaxITSolver,espSolver);
        case 3
            x_out=Jacobi(A,x,b,MaxITSolver,espSolver);
        case 4
            x_out=GaussSeidelMatrix(A,x,b,MaxITSolver,espSolver);
        case 5
            x_out=GaussSeidel(A,x,b,MaxITSolver,espSolver);
        case 6
            x_out=SymmetricGaussSeidelMatrix(A,x,b,MaxITSolver,espSolver);
        case 7
            x_out=SymmetricGaussSeidel(A,x,b,MaxITSolver,espSolver);
        case 8
            % Get the Over-Relaxation Coe. From User.
            omega=input('Please Enter Over-Relaxation(Omega) for the SOR Solver: ');
            x_out=SORMatrix(A,x,b,omega,MaxITSolver,espSolver);
        case 9
            % Get the Over-Relaxation Coe. From User.
            omega=input('Please Enter Over-Relaxation(Omega) for the SOR Solver: ');
            x_out=SOR(A,x,b,omega,MaxITSolver,espSolver);
        case 10
            % Get the Over-Relaxation Coe. From User.
            omega=input('Please Enter Over-Relaxation(Omega) for the SOR Solver: ');
            x_out=SymmetricSORMatrix(A,x,b,omega,MaxITSolver,espSolver);
        case 11
            x_out=GaussSeidelRedBlack(A,x,b,MaxITSolver,espSolver);
        case 12
            % Get the Over-Relaxation Coe. From User.
            omega=input('Please Enter Over-Relaxation(Omega) for the SOR Solver: ');
            x_out=SOR4Tri(A,x,b,omega,MaxITSolver,espSolver);
        case 13
            x_out=CG(A,x,b,MaxITSolver,espSolver);
        case 14
            B=diag(diag(A)); %Preconditioning Matrix--Simplest
            x_out=PCG(A,x,b,B,MaxITSolver,espSolver);
        case 15
            x_out=BCG(A,x,b,MaxITSolver,espSolver);
        otherwise
            display('Please First Select a solver')
            Pause
            return;
    end
end

