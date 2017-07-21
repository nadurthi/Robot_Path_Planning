% Orbit groundtrack plot Latitude longitude lat long
% Richard Rieber
% December 18, 2006
% rrieber@gmail.com
%
% Revision 8/21/07: Fixed small typo in help comments
%                   Added H1 line for lookfor functionality
%          9/27/09: Added inputs for:
%                    - plot string, s, to customize plot
%                    - mu, gravitational constant
%                   Also fixed a bug drawing a horizontal line across the
%                   plot when the orbit propogates from east to west
%                   Removed Re, since it is not used.
%
% function Groundtrack(Kepler,GMSTo,Tf,fig,dT,s,mu)
%
% Purpose: This function plots the groundtrack of a given orbit.
% 
% Inputs:  o Kepler - A vector of length 6 containing all of the keplerian
%                    orbital elements [a,ecc,inc,Omega,w,nu] in km and radians
%          o GMSTo  - The Greenwich Mean Siderial Time at the given initial position
%                     in radians.
%          o Tf     - The length of time to plot the groundtrack in seconds
%          o fig    - The figure number on which to plot the groundtrack [OPTIONAL]
%          o dT     - Timesteps for groundtrack, defaults to 60 seconds [OPTIONAL]
%          o s      - String for customizing the plot (example: '--b*').
%                     See, help plot, for more information [OPTIONAL]
%          o mu     - Gravitational constant of Earth. Defaults to
%                     3998600.4415 km^3/s^2 [OPTIONAL]
%
% Outputs: o H      - The figure handle

function [h] = Groundtrack(Kepler,GMSTo,Tf,fig,dT,s,mu)

h = 0;

r2d = 180/pi;  %Radians to degrees conversion
w_earth  = 7.2921158553e-5; %rad/s  %Rotation rate of Earth

if nargin < 3 || nargin > 7
    error('Incorrect number of inputs, see help groundtrack.m')
elseif nargin == 3
	h = figure;
    dT  = 60;
    s   = 'b';
    mu  = 398600.4415;     %km^3/s^2  Gravitational constant of Earth
elseif nargin == 4
    s = 'b';
    dT  = 60;
    s   = 'b';
    mu  = 398600.4415;     %km^3/s^2  Gravitational constant of Earth
elseif nargin == 5
    s   = 'b';
    mu  = 398600.4415;     %km^3/s^2  Gravitational constant of Earth
elseif nargin == 6
    mu  = 398600.4415;     %km^3/s^2  Gravitational constant of Earth
end

if h == 0
    h = figure(fig);
end

a   = Kepler(1);
ecc = Kepler(2);
inc = Kepler(3);
O   = Kepler(4);
w   = Kepler(5);
nuo = Kepler(6);

n = (mu/a^3)^.5;  %Mean motion of satellite

E  = atan2(sin(nuo)*(1-ecc^2)^.5,ecc+cos(nuo));  %Eccentric anomaly
MA = E-ecc*sin(E); %Initial Mean anomaly

time = [0:dT:Tf]; %time vector

bp(1)    = 0;
bp(2)    = length(time);
k        = 2;
Lat_rad  = zeros(1,length(time));
Long_rad = zeros(1,length(time));

for j = 1:length(time)
    GMST = zeroTo360(GMSTo + w_earth*time(j),1);  %GMST in radians
    
    M  = zeroTo360(MA + n*time(j),1);  %Mean anomaly in rad
    nu = nuFromM(M,ecc,10^-12);        %True anomaly in rad
    
    [ECI,Veci] = randv(a,ecc,inc,O,w,nu);
    clear veci
    ECEF = eci2ecef(ECI,GMST);
    
    [Lat_rad(j),Long_rad(j)] = RVtoLatLong(ECEF);
    
    if j > 1 && ((Long_rad(j)-Long_rad(j-1)) < -pi || (Long_rad(j)-Long_rad(j-1)) > pi)
        bp(k) = j-1;
        k     = k+1;
    end
end

bp(length(bp)+1) = length(time);
Lat              = Lat_rad.*r2d;
Long             = Long_rad.*r2d;

x = load('Coastline.dat');
plot(x(:,1),x(:,2),'g')
clear x

hold on
for j = 2:length(bp)
	plot(Long(bp(j-1)+1:bp(j)),Lat(bp(j-1)+1:bp(j)),s)
end

axis([-180,180,-90,90])

xlabel('Longitude (\circ)')
ylabel('Latitude (\circ)')