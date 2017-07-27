% 2 body problem pure propagation
% function Sims_2bodyProb
clear
clc
close
global kappa

time.t0 = 0;
time.dt = 29.173785213265145;
time.tf = 21600;
T=time.t0:time.dt:time.tf;
Re=6378.1;

model.fn = 6;               % state space dimensionality
model.fx = 'twoBody';

P0 =blkdiag(0.01,0.01,0.01,0.000001,0.000001,0.000001);
r0 = [7000 0 0 0 -1.0374090357 7.4771288355]';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=10000;
wmc=ones(N,1)./N;
Xmc0=mvnrnd(r0,P0,N);
opt = odeset('reltol',1e-10,'abstol',1e-10);
Xmc=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,Xmc0(i,:)',opt);
   Xmc(:,:,i)=X;
end
save('body2sims','Xmc','wmc')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kappa=1;
[Xut0,wut]=UT_sigmapoints(r0,P0,2);

N=length(wut);
Xut=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,Xut0(i,:)',opt);
   Xut(:,:,i)=X;
end
save('body2sims','Xmc','wmc','Xut','wut')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xcut40,wcut4]=conjugate_dir_gausspts(r0,P0);
N=length(wcut4);
Xcut4=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,Xcut40(i,:)',opt);
   Xcut4(:,:,i)=X;
end
save('body2sims','Xmc','wmc','Xut','wut','Xcut4','wcut4')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xcut60,wcut6]=conjugate_dir_gausspts_till_6moment_scheme2(r0,P0);
N=length(wcut6);
Xcut6=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,Xcut60(i,:)',opt);
   Xcut6(:,:,i)=X;
end
save('body2sims','Xmc','wmc','Xut','wut','Xcut4','wcut4','Xcut6','wcut6')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xcut80,wcut8]=conjugate_dir_gausspts_till_8moment(r0,P0);
N=length(wcut8);
Xcut8=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,Xcut80(i,:)',opt);
   Xcut8(:,:,i)=X;
end
save('body2sims','Xmc','wmc','Xut','wut','Xcut4','wcut4','Xcut6','wcut6','Xcut8','wcut8')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xgh30,wgh3]=GH_points(r0,P0,3);
N=length(wgh3);
Xgh3=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,Xgh30(i,:)',opt);
   Xgh3(:,:,i)=X;
end
save('body2sims','Xmc','wmc','Xut','wut','Xcut4','wcut4','Xcut6','wcut6','Xcut8','wcut8','Xgh3','wgh3')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xgh40,wgh4]=GH_points(r0,P0,4);
N=length(wgh4);
Xgh4=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,Xgh40(i,:)',opt);
   Xgh4(:,:,i)=X;
end
save('body2sims','Xmc','wmc','Xut','wut','Xcut4','wcut4','Xcut6','wcut6','Xcut8','wcut8','Xgh3','wgh3','Xgh4','wgh4')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xgh50,wgh5]=GH_points(r0,P0,5);
N=length(wgh5);
Xgh5=zeros(length(T),6,N);
parfor i=1:N
    i
   [t,X]= ode45(@twoBody,T,Xgh50(i,:)',opt);
   Xgh5(:,:,i)=X;
end

save('body2sims','Xmc','wmc','Xut','wut','Xcut4','wcut4','Xcut6','wcut6','Xcut8','wcut8','Xgh3','wgh3','Xgh4','wgh4','Xgh5','wgh5')


X=Xmc./Re;
X1=X(:,:,1:50000);
% X(:,4:6,:)=[];
w1=ones(50000,1)./50000
[y1,Mmc1]=Evol_moments_samples(X1,w1,1,'central');
[y2,Mmc2]=Evol_moments_samples(X1,w1,2,'central');
[y3,Mmc3]=Evol_moments_samples(X1,w1,3,'central');
[y4,Mmc4]=Evol_moments_samples(X1,w1,4,'central');

X=Xut./Re;
X(:,4:6,:)=[];
[y1,Mut1]=Evol_moments_samples(X,wut,1,'central');
[y2,Mut2]=Evol_moments_samples(X,wut,2,'central');
[y3,Mut3]=Evol_moments_samples(X,wut,3,'central');
[y4,Mut4]=Evol_moments_samples(X,wut,4,'central');

X=Xcut4./Re;
X(:,4:6,:)=[];
[y1,Mcut41]=Evol_moments_samples(X,wcut4,1,'central');
[y2,Mcut42]=Evol_moments_samples(X,wcut4,2,'central');
[y3,Mcut43]=Evol_moments_samples(X,wcut4,3,'central');
[y4,Mcut44]=Evol_moments_samples(X,wcut4,4,'central');

X=Xcut6./Re;
X(:,4:6,:)=[];
[y1,Mcut61]=Evol_moments_samples(X,wcut6,1,'central');
[y2,Mcut62]=Evol_moments_samples(X,wcut6,2,'central');
[y3,Mcut63]=Evol_moments_samples(X,wcut6,3,'central');
[y4,Mcut64]=Evol_moments_samples(X,wcut6,4,'central');

X=Xcut8./Re;
X(:,4:6,:)=[];
[y1,Mcut81]=Evol_moments_samples(X,wcut8,1,'central');
[y2,Mcut82]=Evol_moments_samples(X,wcut8,2,'central');
[y3,Mcut83]=Evol_moments_samples(X,wcut8,3,'central');
[y4,Mcut84]=Evol_moments_samples(X,wcut8,4,'central');

X=Xgh3;
X(:,4:6,:)=[];
[y1,Mgh31]=Evol_moments_samples(X,wgh3,1,'central');
[y2,Mgh32]=Evol_moments_samples(X,wgh3,2,'central');
[y3,Mgh33]=Evol_moments_samples(X,wgh3,3,'central');
[y4,Mgh34]=Evol_moments_samples(X,wgh3,4,'central');

X=Xgh4;
X(:,4:6,:)=[];
[y1,Mgh41]=Evol_moments_samples(X,wgh4,1,'central');
[y2,Mgh42]=Evol_moments_samples(X,wgh4,2,'central');
[y3,Mgh43]=Evol_moments_samples(X,wgh4,3,'central');
[y4,Mgh44]=Evol_moments_samples(X,wgh4,4,'central');

X=Xgh5./Re;
X(:,4:6,:)=[];
[y1,Mgh51]=Evol_moments_samples(X,wgh5,1,'central');
[y2,Mgh52]=Evol_moments_samples(X,wgh5,2,'central');
[y3,Mgh53]=Evol_moments_samples(X,wgh5,3,'central');
[y4,Mgh54]=Evol_moments_samples(X,wgh5,4,'central');

save('body2sims_final','Xmc','wmc','Xut','wut','Xcut4','wcut4','Xcut6','wcut6','Xcut8','wcut8','Xgh3','wgh3','Xgh4','wgh4','Xgh5','wgh5','y1','y2','y3','y4','Mmc1','Mmc2','Mmc3','Mmc4','Mut1','Mut2','Mut3','Mut4','Mcut41','Mcut42','Mcut43','Mcut44','Mcut61','Mcut62','Mcut63','Mcut64','Mcut81','Mcut82','Mcut83','Mcut84','Mgh31','Mgh32','Mgh33','Mgh34','Mgh41','Mgh42','Mgh43','Mgh44','Mgh51','Mgh52','Mgh53','Mgh54')




err1_mc1e4=sqrt(sum(((Mmc1_1e4-Mmc1)).^2,2));
err1_ut   =sqrt(sum(((Mut1-Mmc1)).^2,2));
err1_cut4 =sqrt(sum(((Mcut41-Mmc1)).^2,2));
err1_cut6 =sqrt(sum(((Mcut61-Mmc1)).^2,2));
err1_cut8 =sqrt(sum(((Mcut81-Mmc1)).^2,2));
 err1_gh5 =sqrt(sum(((Mgh51-Mmc1)).^2,2));

err2_mc1e4=(sum(abs((Mmc2_1e4-Mmc2)).^(1/2),2));
err2_ut   =(sum(abs((Mut2-Mmc2)).^(1/2),2));
err2_cut4 =(sum(abs((Mcut42-Mmc2)).^(1/2),2));
err2_cut6 =(sum(abs((Mcut62-Mmc2)).^(1/2),2));
err2_cut8 =(sum(abs((Mcut82-Mmc2)).^(1/2),2));
err2_gh5 =(sum(abs((Mgh52-Mmc2)).^(1/2),2));

err3_mc1e4=(sum(abs((Mmc3_1e4-Mmc3)).^(1/3),2));
err3_ut   =(sum(abs((Mut3-Mmc3)).^(1/3),2));
err3_cut4 =(sum(abs((Mcut43-Mmc3)).^(1/3),2));
err3_cut6 =(sum(abs((Mcut63-Mmc3)).^(1/3),2));
err3_cut8 =(sum(abs((Mcut83-Mmc3)).^(1/3),2));
 err3_gh5 =(sum(abs((Mgh53-Mmc3)).^(1/3),2));

err4_mc1e4=(sum(abs((Mmc4_1e4-Mmc4)).^(1/4),2));
err4_ut   =(sum(abs((Mut4-Mmc4)).^(1/4),2));
err4_cut4 =(sum(abs((Mcut44-Mmc4)).^(1/4),2));
err4_cut6 =(sum(abs((Mcut64-Mmc4)).^(1/4),2));
err4_cut8 =(sum(abs((Mcut84-Mmc4)).^(1/4),2));
err4_gh5 =(sum(abs((Mgh54-Mmc4)).^(1/4),2));

[max(err1_mc1e4),max(err1_ut),max(err1_cut4),max(err1_cut6),max(err1_cut8),max(err1_gh5);...
 max(err2_mc1e4),max(err2_ut),max(err2_cut4),max(err2_cut6),max(err2_cut8),max(err2_gh5);...
 max(err3_mc1e4),max(err3_ut),max(err3_cut4),max(err3_cut6),max(err3_cut8),max(err3_gh5);...
 max(err4_mc1e4),max(err4_ut),max(err4_cut4),max(err4_cut6),max(err4_cut8),max(err4_gh5)]./10

% [norm(err1_ut),norm(err1_cut4),norm(err1_cut6),norm(err1_cut8);...
%  norm(err2_ut),norm(err2_cut4),norm(err2_cut6),norm(err2_cut8);...
%  norm(err3_ut),norm(err3_cut4),norm(err3_cut6),norm(err3_cut8);...
%  norm(err4_ut),norm(err4_cut4),norm(err4_cut6),norm(err4_cut8)]

t=T;
figure
plot(t,sqrt(sum((Mmc1_1e4-Mmc1).^2,2)),'--',t,sqrt(sum((Mut1-Mmc1).^2,2)),t,sqrt(sum((Mcut41-Mmc1).^2,2)),t,sqrt(sum((Mcut61-Mmc1).^2,2)),t,sqrt(sum((Mcut81-Mmc1).^2,2)))
legend('MC-10000','UT','CUT4','CUT6','CUT8')

figure
plot(t,sqrt(sum((Mmc1_1e4-Mmc1).^2,2)),'--',t,sqrt(sum((Mut1-Mmc1).^2,2)),t,sqrt(sum((Mcut41-Mmc1).^2,2)),t,sqrt(sum((Mcut61-Mmc1).^2,2)),t,sqrt(sum((Mcut81-Mmc1).^2,2)))
legend('MC-10000','UT','CUT4','CUT6','CUT8')

figure
plot(t,sqrt(sum((Mmc1_1e4-Mmc1).^2,2)),'--',t,sqrt(sum((Mut1-Mmc1).^2,2)),t,sqrt(sum((Mcut41-Mmc1).^2,2)),t,sqrt(sum((Mcut61-Mmc1).^2,2)),t,sqrt(sum((Mcut81-Mmc1).^2,2)))
legend('MC-10000','UT','CUT4','CUT6','CUT8')

figure
plot(t,sqrt(sum((Mmc1_1e4-Mmc1).^2,2)),'--',t,sqrt(sum((Mut1-Mmc1).^2,2)),t,sqrt(sum((Mcut41-Mmc1).^2,2)),t,sqrt(sum((Mcut61-Mmc1).^2,2)),t,sqrt(sum((Mcut81-Mmc1).^2,2)))
legend('MC-10000','UT','CUT4','CUT6','CUT8')

% figure
% plot3(Mmc1(:,1),Mmc1(:,2),Mmc1(:,3),Xmc0(:,1),Xmc0(:,2),Xmc0(:,3),'ro')
% view(130,38)
% grid
% xlabel('x')
% ylabel('y')
% zlabel('z')
% axis([-1e4,1e4,-1e3,1e3,-1e4,1e4])

% figure
% plot(T,Mmc2(:,6),T,Mmc2(:,19))
% % view(130,38)
% % grid
% xlabel('x')
% ylabel('y')
% % zlabel('z')
% % axis([-1e4,1e4,-1e3,1e3,-1e4,1e4])

% for Tp=[601,699]
% Y=zeros(N,3);
% Y(:,1)=Xmc(Tp,1,:);
% Y(:,2)=Xmc(Tp,2,:);
% Y(:,3)=Xmc(Tp,3,:);

% figure
% 
% subplot(1,3,1)
%  hist3([Y(:,1),Y(:,2)])
%  view(122,16)
%  xlabel('x')
%  ylabel('y')
%  
%  subplot(1,3,2)
%  hist3([Y(:,2),Y(:,3)])
%  view(122,16)
%  xlabel('y')
%  ylabel('z')
%  
%  subplot(1,3,3)
%  hist3([Y(:,1),Y(:,3)])
%  view(122,16)
%  xlabel('x')
%  ylabel('z')
%  
% figure
% subplot(1,3,1)
%  plot(Y(:,1),Y(:,2),'ro')
% %  view(122,16)
%  xlabel('x')
%  ylabel('y')
%  
%  subplot(1,3,2)
%  plot(Y(:,2),Y(:,3),'ro')
% %  view(122,16)
%  xlabel('y')
%  ylabel('z')
%  
%  subplot(1,3,3)
%  plot(Y(:,1),Y(:,3),'ro')
% %  view(122,16)
%  xlabel('x')
%  ylabel('z')
%  title(num2str(Tp))
%  pause(1)
% end
 
% end
% 
% y=[zeros(1,3);eye(3);y2;y3]%;y4];
% M=[1;zeros(3,1);Mmc2(Tp,:)';Mmc3(Tp,:)']%;Mmc4(Tp,:)'];
% [y,lam,xl,xu]=MaxEntPdf(y,M,10);
% 
% [X,W] = GLeg_pts(15, xl(3), xu(3));
% 
% [xx,zz]=meshgrid(xl(1):1:xu(1),xl(3):1:xu(3));
% p=zeros(size(xx));
% for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%         for k=1:1:length(W)
% p(i,j)=p(i,j)+W(k)*pdf_MaxEnt([xx(i,j),X(k),zz(i,j)],lam,y);
%         end
%     end
% end
% [xx,zz]=meshgrid(xl(1)+Mmc1(Tp,1):1:xu(1)+Mmc1(Tp,1),xl(3)+Mmc1(Tp,3):1:xu(3)+Mmc1(Tp,3));
% 
%  plot(Y(:,1),Y(:,3),'ro')
% %  view(122,16)
%  xlabel('x')
%  ylabel('z')
%  hold on
%  
%  contour(xx,zz,p,27)
% 
% options = statset('MaxIter',5000); 
%  obj = gmdistribution.fit([Y(:,1),Y(:,3)],5,'Regularize',0.1,'Options',options);
%  
%  [xx,zz]=meshgrid(xl(1)+Mmc1(Tp,1):1:xu(1)+Mmc1(Tp,1),xl(3)+Mmc1(Tp,3):1:xu(3)+Mmc1(Tp,3));
%  p=zeros(size(xx));
%  for i=1:1:size(xx,1)
%     for j=1:1:size(xx,2)
%  p(i,j)=pdf(obj,[xx(i,j),zz(i,j)]);
%     end
%  end
%  
%   plot(Y(:,1),Y(:,3),'ro')
% %  view(122,16)
%  xlabel('x')
%  ylabel('z')
%  hold on
%   contour(xx,zz,p,5)