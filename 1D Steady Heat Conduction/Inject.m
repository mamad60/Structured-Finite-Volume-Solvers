function [ output ] = Inject( n )
%Restricts(From fine to Coarse grid) by injection-1D
%n is lenght of the input row vector
%Output is the matrix for injection
% output size is half of input--->Ih_2h
nn=floor(n/2);
Ih_2h=zeros(nn,n);
%Set the Restriction Matrix
j=1;
for i=1:nn
    Ih_2h(i,j)=1;
    j=j+2;
end
output=Ih_2h;
end

