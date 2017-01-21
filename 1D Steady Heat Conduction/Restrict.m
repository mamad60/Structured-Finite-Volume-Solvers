function [ output ] = Restrict( n )
%Restricts(From fine to Coarse grid) by Full Weight Interpolation-1D
%n is lenght of the input row vector
%Output is the matrix for injection
% output size is half of input--->Ih_2h
nn=floor(n/2);
Ih_2h=zeros(nn,n);
%Set the Restriction Matrix
Ih_2h(1,1)=1;  %Unchanged Boundaries
Ih_2h(nn,n)=1; %Unchanged Boundaries
j=2;
for i=2:nn-1
    Ih_2h(i,j)=0.5;
    Ih_2h(i,j+1)=0.25;
    Ih_2h(i,j-1)=0.25;   
    j=j+2;
end
output=Ih_2h;
end

