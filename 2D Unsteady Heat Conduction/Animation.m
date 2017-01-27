function  Animation(IT,hT,hQ,Dt)
%Creates Animations During the Solution
%Animation my be slow the solution
global animateT animateQ intAnimT intAnimQ intSaveT intSaveQ movT movQ qx qy q_l q_r  X Y T mCountT mCountQ
%Animation
%----Plot  T Countours
if animateT
    if mod(IT,intAnimT)==0 || IT==1
        figure(hT);
        %----Plot  T Contour
        contourf(X,Y,T,'LineWidth',2);
        xlabel('X Coordinate');
        ylabel('Y Coordinate');
        title(strcat('Temperature Contour  ','  (','Time Step=',num2str(IT),' ,Time= ',num2str(Dt*(IT-1)),')'));
        drawnow
        pause(0.5);
        movT(mCountT)=getframe(gcf);
        mCountT=mCountT+1;
        if intSaveT, print(hT,'-dtiff',strcat('Tconntour',num2str(IT))) ,end
    end
end
%----Plot  Q Profile
if animateQ % Use with precuation, Slow down the code very much
    if mod(IT,intAnimQ)==0 || IT==1
        %---Compute Fluxes @ Boundary Points
        %Left & Right Boundaries
        for j=1:n;
            q_l(j)=ke(1,j)*(T(1,j)-T(2,j))/abs(X(1,j)-X(2,j))*Dy(1,j);   %Flux @ the Left Boundary
            q_r(j)=kw(m,j)*(T(m-1,j)-T(m,j))/abs(X(m-1)-X(m,j))*Dy(m,j);   %Flux @ the Left Boundary
        end
        %Top & Bottom Boundaries
        for i=1:m;
            q_t(i)=kn(i,n)*(T(i,n-1)-T(i,n))/abs(Y(i,n-1)-Y(i,n))*Dx(i,n);   %Flux @ the Bottom Boundary
            q_b(i)=ks(i,1)*(T(i,1)-T(i,2))/abs(Y(i,1)-Y(i,2))*Dx(i,2);   %Flux @ the Top Boundary
        end
        %-----Compute and Plot Fluxes @ East & South CV Faces
        for j=1:n;
            for i=2:m-1;
                qx(i,j)=ke(i,j)*(T(i,j)-T(i+1,j))/abs(X(i+1,j)-X(i,j))*Dy(i,j);
            end
            qx(1,j)=q_l(j);
            qx(m,j)=q_r(j);
        end
        for i=1:m;
            for j=2:n-1;
                qy(i,j)=ks(i,j)*(T(i,j)-T(i,j+1))/abs(Y(i,j+1)-Y(i,j))*Dx(i,j);
            end
            qy(i,1)=q_t(i);
            qy(i,n)=q_b(i);
        end
        figure(hQ);
        %Plot Flux vectors
        contour(X,Y,T)
        hold on
        quiver(X,Y,qx,qy)
        colormap hsv
        hold off
        title(strcat('Profile of Heat Flux  ','  (','Time Step=',num2str(IT-1),' ,Time= ',num2str(Dt*(IT-1)),')'));
        drawnow
        pause(0.5);
        movQ(mCountQ)=getframe(gcf);
        mCountQ=mCountQ+1;
        if intSaveQ, print(hQ,'-dtiff',strcat('Qprofile',num2str(IT))) ,end
    end
end




