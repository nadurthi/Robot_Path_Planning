function xk1=UM_motion(T,xk)
% global T;
% T=para(1);
% T=1;
% keyboard

xk=xk(:);

xk1=[1,0,T,0;
    0,1,0,T;
    0,0,1,0;
    0,0,0,1]*xk;
end
