function [ q ] = Euler( input, Dt, RHS )
%Explicit Euler Iteration
%input is value of q in current time
q=input+RHS*Dt;
end

