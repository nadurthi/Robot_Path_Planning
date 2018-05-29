function xk1=dummymodel(tk,xk,dt)

xk1(1)=xk(1)+dt*xk(2);
xk1(2)=xk(2)+dt*(-xk(1)-xk(1)^3);



end
