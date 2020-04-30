%% Yo-yo De-spin Simulation and Analysis
% Daniel Huang
% AE 140, Team/Individual Project
%% Setup 
clc;clear;close all
global m me r h Ix Iz T G Re mw mr l

% Tunable variables
w = 15.5; % Mass of both weights [kg]
l_max = 31.2; % Length of string [m]

init_psi = deg2rad(-120.5724); % Initial longitude
init_the = deg2rad(90-34.7420); % Initial latitude
init_spin = pi; % Initial spin of rocket [rad/s]

r = 3.7/2; % Radius of rocket [m]
m1 = 3.5E5; % Mass of rocket [kg]
h1 = 70; % Height of rocket [m]
T1 = 7.607E6; % First stage thrust [N]
m2 = 5E4; % Mass of second stage [kg]
h2 = 30; % Height of second stage [m]
T2 = 9.34E5; % Second stage thrust [N]

al = 1E-5; % Fixed offset of thrust from nominal axis [rad]

% Other variables
G = 6.67408E-11; % Gravitational constant
Re = 6.378E6; % Radius of Earth [m]
me = 5.972E24; % Mass of earth [kg]
rot_e = 2*pi/(24*60*60); % Earths rotation rate [rad/s]


% Rotation Matrices
nRa = @(psi)[cos(psi) -sin(psi) 0;sin(psi) cos(psi) 0;0 0 1];
aRb = @(the)[1 0 0;0 cos(the) -sin(the);0 sin(the) cos(the)];
bRr = @(phi)[cos(phi) -sin(phi) 0;sin(phi) cos(phi) 0;0 0 1];

%% Setup Initial Launch
m = m1; % Mass of rocket [kg]
h = h1; % Height of rocket [m]
Ix = 0.5*m*r^2; % Ixx/Iyy of rocket about CM [kg*m^2]
Iz = m/12*(3*r^2+h^2); % Izz of rocket about CM [kg*m^2]
T = T1; % Thrust [N]
tspan = [0 150];

% Initial Conditions
x0 = Re*sin(init_the)*sin(init_psi);
xd0 = Re*rot_e*sin(init_the)*cos(init_psi);
y0 = -Re*sin(init_the)*cos(init_psi);
yd0 = Re*rot_e*sin(init_the)*sin(init_psi);
z0 = Re*cos(init_the);
zd0 = 0;
psi0 = init_psi; % Precession angle [rad]
psid0 = 0; % Precession rate [rad/s]
the0 = init_the; % Pitch angle [rad]
thed0 = 0; % Pitch rate [rad/s]
phi0 = 0; % Spin angle [rad/s]
phid0 = init_spin; % Spin rate [rad/s]
Initials = [x0 xd0 y0 yd0 z0 zd0 psi0 psid0 the0 thed0 phi0 phid0]; 

%% Simulate Launch
fprintf('Initial launch solving...\n')
burn1 = ode45(@(t,x)eqn(t,x,al),tspan,Initials,[]);
fprintf('Initial launch solved.\n')

%% Setup Delay 
T = 0; % Turn off thrust
% Detach first stage
m = m2; % Mass of Second stage [kg]
h = h2; % Height of second stage [m]
Initials_t = burn1.y(:,end);
tspan_t = [0 30]; 

%% Simulate Delay
fprintf('First to second stage transition solving...\n')
trans = ode45(@(t,x)eqn(t,x,al),tspan_t,Initials_t,[]);
fprintf('First to second stage transition solved.\n')

%% Setup Second Stage burn
Initials2 = trans.y(:,end); % Setup initials after launch
Initials2(7) = Initials2(7)+120*pi/180; % Pitch over maneuver
Initials2(9) = Initials2(9); 
Ix = 0.5*m*r^2; % Ixx/Iyy of rocket about CM [kg*m^2]
Iz = m/12*(3*r^2+h^2); % Izz of rocket about CM [kg*m^2]
T = T2; % Thrust [N]
tspan2 = [0 450]; % Time for second stage burn

%% Simulate Second Stage burn
fprintf('Second stage burn solving...\n')
burn2 = ode45(@(t,x)eqn(t,x,al),tspan2,Initials2,[]); 
fprintf('Second stage burn solved.\n')

%% Setup despin stage 1
mr = 0.2*m; % Loss of weight due to no more fuel
mw = w; % Mass of both weights [kg]
l = l_max; % Maximum length of string [m]

% Initial Conditions
beta0 = burn2.y(11,end)+0.00001;
betad0 = burn2.y(12,end);
phi0 = burn2.y(11,end);
phid0 = burn2.y(12,end);
Initials3 = [burn2.y(1:6,end);beta0;betad0;phi0;phid0];
tspan3 = [0 10];

options = odeset('Events',@(t,x)Stage1Halt(t,x));

%% Simulate despin stage 1
fprintf('Despin stage 1 solving...\n')
despin1 = ode113(@(t,x)EoM1(t,x),tspan3,Initials3,options);
fprintf('Despin stage 1 solved.\n')

%% Setup despin stage 2
Init = despin1.y(:,end);
Initials4 = Init;
Initials4(7:8) = [0 Init(8)*(Init(7)-Init(9))]; % Change beta to gamma initials
tspan4 = [0 10];

options2 = odeset('Events',@(t,x)Stage2Halt(t,x));

%% Simulate despin stage 2
fprintf('Despin stage 2 solving...\n')
despin2 = ode113(@(t,x)EoM2(t,x),tspan4,Initials4,options2);
fprintf('Despin stage 2 solved.\n')

%% Setup final trajectory
offset = despin1.y(7,end)-despin1.y(9,end); % Offset angle btw cx and rx 
gamd = despin2.y(8,end);
phid = despin2.y(10,end);
% Position and Velocity vectors of weights 
posr = (r+l)*[sin(offset);cos(offset);0]; % R frame position
velr = (l*(gamd+phid)+phid*r)*[cos(offset);-sin(offset);0]; % R frame vel
% Position and Velocity vectors rotated to N frame
pos = nRa(burn2.y(7,end))*aRb(burn2.y(9,end))*bRr(despin2.y(9,end))*posr;
vel = nRa(burn2.y(7,end))*aRb(burn2.y(9,end))*bRr(despin2.y(9,end))*velr;
% Upper stage Initials
Initials5{1} = despin2.y(1:6,end); % No translational modifications
% Weight 1 Initials
Initials5{2}(1) = Initials5{1}(1)+pos(1);
Initials5{2}(3) = Initials5{1}(3)+pos(2);
Initials5{2}(5) = Initials5{1}(5)+pos(3);
Initials5{2}(2) = Initials5{1}(2)+vel(1);
Initials5{2}(4) = Initials5{1}(4)+vel(2);
Initials5{2}(6) = Initials5{1}(6)+vel(3);
% For weight 2
Initials5{3}(1) = Initials5{1}(1)+pos(1);
Initials5{3}(3) = Initials5{1}(3)+pos(2);
Initials5{3}(5) = Initials5{1}(5)+pos(3);
Initials5{3}(2) = Initials5{1}(2)-vel(1);
Initials5{3}(4) = Initials5{1}(4)-vel(2);
Initials5{3}(6) = Initials5{1}(6)-vel(3);
% First stage
Initials5{4} = burn1.y(1:6,end);

tspan5 = [0 100000];

options3 = odeset('Events',@(t,x)HitEarth(t,x));

%% Simulate final trajectory
names = {'upper stage' 'weight 1' 'weight 2' 'first stage'};
for ii = 1:length(names)
    fprintf('Free fall trajectory for %s solving...\n',names{ii})
    final{ii} = ode113(@(t,x)EoMFree(t,x),tspan5,Initials5{ii},options3);
    fprintf('Free fall trajectory for %s solved.\n',names{ii})
end

%% Process Results
% Total trajectory up to end of despin
nx = [burn1.y(1,:) trans.y(1,:) burn2.y(1,:) despin1.y(1,:) despin2.y(1,:)];
ny = [burn1.y(3,:) trans.y(3,:) burn2.y(3,:) despin1.y(3,:) despin2.y(3,:)];
nz = [burn1.y(5,:) trans.y(5,:) burn2.y(5,:) despin1.y(5,:) despin2.y(5,:)];
% Trajectory of Upper stage
posx{1} = [nx final{1}.y(1,:)];
posy{1} = [ny final{1}.y(3,:)];
posz{1} = [nz final{1}.y(5,:)];
% Trajectory of Weight 1
posx{2} = [nx final{2}.y(1,:)];
posy{2} = [ny final{2}.y(3,:)];
posz{2} = [nz final{2}.y(5,:)];
% Trajectory of Weight 2
posx{3} = [nx final{3}.y(1,:)];
posy{3} = [ny final{3}.y(3,:)];
posz{3} = [nz final{3}.y(5,:)];
% Trajectory of First stage
posx{4} = [burn1.y(1,:) final{4}.y(1,:)];
posy{4} = [burn1.y(3,:) final{4}.y(3,:)];
posz{4} = [burn1.y(5,:) final{4}.y(5,:)];

% Obtain time vectors
%       Common time
t = [burn1.x burn1.x(end)+trans.x burn1.x(end)+trans.x(end)+burn2.x ...
    burn1.x(end)+trans.x(end)+burn2.x(end)+despin1.x ...
    burn1.x(end)+trans.x(end)+burn2.x(end)+despin1.x(end)+despin2.x]; 

t_total{1} = [t t(end)+final{1}.x]; % Second stage of rocket
t_total{2} = [t t(end)+final{2}.x]; % First weight 
t_total{3} = [t t(end)+final{3}.x]; % Second weight
t_total{4} = [burn1.x burn1.x(end)+final{4}.x]; % First stage of rocket

% String position in B frame 
l_v = [r*(despin1.y(7,:)-despin1.y(9,:)) l*ones(1,length(despin2.x))];
gam = [zeros(1,length(despin1.x)) despin2.y(7,:)];
beta = [despin1.y(7,:) despin2.y(9,:)+offset];
sBx = (r+l_v.*sin(gam)).*cos(beta)+l_v.*cos(gam).*sin(beta);
sBy = (r+l_v.*sin(gam)).*sin(beta)-l_v.*cos(gam).*cos(beta);
pivot = [r*cos(beta);r*sin(beta)];
t_local = [despin1.x despin1.x(end)+despin2.x];
phi = [despin1.y(9,:) despin2.y(9,:)];

% Earth Sphere
[a,b,c] = sphere;
a = Re*a;
b = Re*b;
c = Re*c;

% Rocket Circle
cx = r*cos(0:0.01:2*pi); % x coords
cy = r*sin(0:0.01:2*pi); % y coords
events = zeros(4,1);
event_loc = zeros(4,1);

%% Static Plot Results
close all
% Trajectory plot
figure
plot3(nx(1),ny(1),nz(1),'O','DisplayName','Launch location','MarkerSize',...
    6,'LineWidth',3)
hold on
plot3(posx{1},posy{1},posz{1},'g','DisplayName','Upper stage') % Upper Stage Trajectory
plot3(posx{2},posy{2},posz{2},'b','DisplayName','Weight 1') % Weight 1 trajectory
plot3(posx{3},posy{3},posz{3},'m','DisplayName','Weight 2') % Weight 2 trajectory
plot3(posx{4},posy{4},posz{4},':k','DisplayName','First stage') % First stage trajectory
% Impact positions
for ii = 1:4
    if final{ii}.ie
        name = sprintf('Impact position of %s\nt = %.2f s',names{ii},t_total{ii}(end));
        plot3(final{ii}.ye(1),final{ii}.ye(3),final{ii}.ye(5),'x',...
            'DisplayName',name,'MarkerSize',9,'Linewidth',3);
        % Both weights impact Earth, upper stage maintains proper orbit.  
    end
end
n = mesh(a,b,c);
set(get(get(n,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
set(gca,'Color','k')
%set(gcf,'Color','k')
axis equal
legend('Location','Best','Color','w')
title('Trajectory Plot')
xlabel('x [m]')
ylabel('y [m]')
zlabel('z [m]')

% Spin rate of upper stage 
figure(2)
plot(despin1.x,despin1.y(10,:),despin1.x(end)+despin2.x,despin2.y(10,:),...
    t_local(end),despin2.y(10,end),'kx','MarkerSize',8)
title('Spin Rate of Upper Stage during Yo-Yo De-spin Maneuver')
xlabel('time [s]')
ylabel('\phi'' [rad/s]')
endrot = sprintf('Ending spin rate: %f rad/s',despin2.y(10,end));
legend('Stage 1','Stage 2',endrot,'Location','Best')

%% Despin animation plot (In B x-y plane)
% figure(3)
% plot(cx,cy,'k') % Uncomment all if saving GIF file
% title('De-spin Mechanism Animation')
% xlabel('Bx [m]')
% ylabel('By [m]')
% ylim([-40 40])
% axis equal
% legend('Rocket','Location','East')
% gif('Despin_animation.gif','DelayTime',1/10,'LoopCount',1) % Save animation as GIF file

for ii = 1:length(l_v)
    figure(3)
    if gam(ii)
        stage = 'Stage 2';
    else
        stage = 'Stage 1';
    end
    plot(cx,cy,'k',[pivot(1,ii) sBx(ii)],[pivot(2,ii) sBy(ii)],...
        'r',-[pivot(1,ii) sBx(ii)],-[pivot(2,ii) sBy(ii)],'b',...
        [0 r*cos(phi(ii))],[0 r*sin(phi(ii))],'.-k')
    t_str = sprintf('De-spin Mechanism Animation: %s',stage);
    title(t_str)
    xlabel('Bx [m]')
    ylabel('By [m]')
    ylim([-40 40])
    axis equal
    legstr = sprintf('Rocket\nt = %f s',t_local(ii));
    legend(legstr,'Weight 1','Weight 2','Rocket Reference','Location','East')
    drawnow
    % pause(0.2)
    %gif % Uncomment if saving GIF file
end
