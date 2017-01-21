function [ output ] = ProlongI( n )
%Prolongation(From Coarse to Fine grid) by Interpolation-1D
%input is a row vector-size n---> returns Matrix for I(2h,h)
% output size is half of input--->Ih_2h
I2h_h=zeros(2*n,n);
%Set the Restriction Matrix
i=2;
for j=2:n
    I2h_h(i,j-1)=0.5;
    I2h_h(i,j)=0.5;
    I2h_h(i-1,j-1)=1;
    i=i+2;
end
I2h_h(2*n-1,n)=1;
output=I2h_h;

end

