close all
global IC FI FIExact XC Xmin Xmax N P
L=zeros(5,12);
for IC=0:4
    Main1DConvDif1stOrder
    L(IC+1,:)=FI;
end

plot(XC,FIExact,'-go',XC,L(1,:),'^:b',XC,L(2,:),'s-.g',XC,L(3,:),'d--r',...
   XC,L(5,:),'d:c','LineWidth',2)
xlabel('X')
ylabel('\phi')
xlim([Xmin Xmax])
ylim([min(FIExact) max(FIExact)])
title('Comparison Between Convection Term Discritisation Schemes')
legend('Exact','Centeral Difference','First Order Upwind'...  
,'Hybrid','Power Law','Location','NorthWest')
text(0.03,0.6,strcat('N=',num2str(N)),'Units','Normalized','Edge','blue') 
text(0.03,0.8,strcat('Peclet=',num2str(P)),'Units','Normalized','Edge','red') 
