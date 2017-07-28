function Xsat0=getInitialrv_3D_3(n)


% global mue re;

re = 6372.797;   % Mean earth radius (Km)
mue = 398601.2;  % Gravitational parameter of earth (Km^3/Sec^2)


r0 = [7000; 0; 0];
rdot0 = [0; -1.03741; 7.977129];

OE = XYZ2OE([r0(:);rdot0(:)]);
a0=OE(1);
e0=OE(2);
E0=OE(3);
w0=OE(4);
i0=OE(5);
Om0=OE(6);

M0 = E0 - e0*sin(E0);

Nsat = n;

% uniform distribution parameters
M0a = 0; M0b = pi/5;

M0sat = M0 + M0a + (M0b-M0a)*rand(Nsat,1);
inc=i0-pi/8+(i0+pi/8-(i0-pi/8))*rand(Nsat,1);
% inc=i0*ones(Nsat,1);
for j = 1:Nsat
%     E0sat = kepler(M0sat(j), e0, 0, 1e-12);
    E0sat = kepler(e0,M0sat(j));
    XX= OE2XYZ([a0, e0, E0sat, w0, inc(j), Om0, mue]);
    X1=XX(1:3);
    Xdot1 =XX(4:6);
    Xhist(:,j) = X1;
    Xdothist(:,j) = Xdot1;
end
Xsat0=[Xhist',Xdothist'];
end
