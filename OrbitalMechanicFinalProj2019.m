% Joseph Poncini
% Final Project Orbit propogation


clc
close all
clear all

mu = 398600;
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

%1 Titan
rtitan0 = [-34679.94608950792, 23714.16356235581, -983.9725080343926];
vtitan0 = [-1.7308940382298856, -2.527619702112526, -0.20633135605217293];

%2 Pegasus
rpegasus0 =  [5655.710092112706, -4112.6760591427455, 193.52641090744842];
vpegasus0 =  [4.429404439223712, 6.062499824675882, -0.76520997255162];

%3 Falcon
rfalcon0 = [-2361.2244281594963, -6512.964230350284, -1082.1106081104983];
vfalcon0 =  [7.109745254222879, -2.5029912865911577, -0.22095098382149542];

%4 Ayame
rayame0 = [522.0863977348827, -20832.59952640893, -2129.5033134102673];
vayame0 = [4.752905144825563, 1.347547151794975, -0.48726994774098664];

%coes
[a1, e1, i1, W1, w1, theta1, h_mag1,rp1,ra1] = find_coes_degrees(rtitan0, vtitan0);
[a2, e2, i2, W2, w2, theta2, h_mag2,rp2,ra2] = find_coes_degrees(rpegasus0, vpegasus0);
[a3, e3, i3, W3, w3, theta3, h_mag3,rp3,ra3] = find_coes_degrees(rfalcon0, vfalcon0);
[a4, e4, i4, W4, w4, theta4, h_mag4,rp4,ra4] = find_coes_degrees(rayame0, vayame0);

[T1] = Tfinder(mu,a1);
[T2] = Tfinder(mu,a2);
[T3] = Tfinder(mu,a3);
[T4] = Tfinder(mu,a4);

%propgating

T = 24*60*60;
initstate1 = [rtitan0 vtitan0];
initstate2 = [rpegasus0 vpegasus0];
initstate3 = [rfalcon0 vfalcon0];
initstate4 = [rayame0 vayame0];

initstates = {initstate1,initstate2,initstate3,initstate4};
t1 = 0;
t2 = T;

[~,States,Times]=propogate(initstates,t1,t2,'on','n','n',1);
%plot3(ra4(1),ra4(2),ra4(3),'go')
%% mission1
t1 = 0 ;
t2 = 5*T2;
[endstates1,states1,times1] = propogate(initstates,t1,t2,'of','n','n',1); %five loops around peg
%%  INC and Phase
kk = 0;

state2 = States{2};
state3 = States{3};
for ii = 1:length(state2)
    for jj = 1:length(state3)
        diff = norm(state2(ii,1:3) - state3(jj,1:3));
       if diff < 13.1626
          index = [ii,jj];
%           disp(state2(ii,1:3))
%             disp(state3(jj,1:3))
            kk = 1;
            break
       end
    end
    if kk == 1
        break
    end
end

r = state2(ii,1:3);
v1 = state2(ii,4:6);
v2 = state3(jj,4:6);
%deltav1 = norm(v2-v1);
theta = acos((dot(rpegasus0,r))/(norm(rpegasus0)*norm(r)));
ratio = theta/(2*pi);
dt =  T2*ratio;
t3 = (5*T2) + dt;


endstate2 = endstates1{2};
rpeg = endstate2(1:3);
vpeg = endstate2(4:6);
rfal = state3(jj,1:3);
vfal = state3(jj,4:6);

[vd,va] = lamberts(mu,rpeg,rfal,dt,1.5,'n');
dv1 = norm(vd - vpeg);
dv2 = norm(vfal - va);
deltavlam1 = dv1+dv2;


[endstates2,states2,times2] = propogate(endstates1,t2,t3,'of','n','n',1); %pegasus to intersections
[ttranslam, stranslam] = ode45(@odefun,[t2 t3],[rpeg vd],options, mu);
plot3(stranslam(:,1),stranslam(:,2),stranslam(:,3),'k')
%% mission

endstate2 = endstates2{2};
rpeg0 = endstate2(1:3);
vpeg0 = endstate2(4:6);

endstate3 = endstates2{3};
rfal0 = endstate3(1:3);
vfal0 = endstate3(4:6);

%finding dt
theta = acos((dot(rfal0,rpeg0))/(norm(rfal0)*norm(rpeg0)));
ratio = theta/(2*pi);
Tphase =  T3 + T3*ratio;
t4 = t3 + Tphase;

[endstates3,states3,times3] = propogate(endstates2,t3,t4,'of','n','n',1); %after phase shift

endstate2 = endstates3{2};
rpeg1 = endstate2(1:3);
vpeg1 = endstate2(4:6);

endstate3 = endstates3{3};
rfal1 = endstate3(1:3);
vfal1 = endstate3(4:6);

%finding dv
aphase = afinder(Tphase,mu);
vmagphase = vfinder1(mu,norm(rfal1),aphase);
vphase = (vmagphase/norm(vfal1))*vfal1;
deltavphase1 = 2*norm(vphase - vfal1);

[ttrans1, strans1] = ode45(@odefun,[t3 t4],[rfal1 vphase],options, mu); %phase shift orbit
plot3(strans1(:,1),strans1(:,2),strans1(:,3),'k')

%% lamberts
t5 = t4 + 5*T3;
[endstates4,states4,times4] = propogate(endstates3,t4,t5,'of','n','n',1); % going around falcon orbit


r1 = rfal1;
% dt = (16)*(60*60) + (20)*60; % 12 and 16 hrs work nice
dt = T4 - (t5-t3) - 8.5*(60) ;
t6 = t5 + dt;
endstate4 = endstates4{4};
[tayame, stateayame] = ode45(@odefun,[t5 t6],endstate4,options, mu);
r2a = stateayame(end,1:3);
v2a = stateayame(end,4:6);
[vd,va] = lamberts(mu,r1,r2a,dt,1.5,'n');
dv1 = norm(vd - vfal1);
dv2 = norm(v2a - va);
deltav2 = dv1+dv2;

[endstates5,states5,times5] = propogate(endstates4,t5,t6,'of','n','n',1);
[ttrans2, strans2] = ode45(@odefun,[t5 t6],[r1 vd],options, mu);
plot3(strans2(:,1),strans2(:,2),strans2(:,3),'k')

%%
t7 = t6 + 5*T4;
[endstates6,states6,times6] = propogate(endstates5,t6,t7,'of','n','n',1); %go around ayame
state1 = states6{1};


alt1 = norm(ra4);
% alt2 = norm(rp1);
alt2 = 4.175126761968658e+04;
a = (alt1 + alt2)/2;
h1 = hfinder3(alt1,e4,pi,mu);
V1 = vfinder1(mu,alt1,a4);
Vd = vfinder1(mu,alt1,a);
Va = vfinder1(mu,alt2,a);
V2 = sqrt(mu/alt2);
Tcir = ((2*pi)/sqrt(mu))*alt2^(3/2);

dV1 = Vd - V1;
dV2 = V2 - Va;


deltav3 = dV1 + dV2;

endstate4 = endstates6{4};


r4 = endstate4(:,1:3);
vayame = endstate4(1,4:6);
vdvector = (Vd/norm(vayame))*vayame;
Tell = Tfinder(mu,a);
Ttransit = Tell/2;
t8 = t7 + Ttransit;

[endstates7,states7,times7] = propogate(endstates6,t7,t8,'of','n','n',1); %during hohmann transfer
[ttrans3, strans3] = ode45(@odefun,[t7 t8],[r4 vdvector],options, mu);
plot3(strans3(:,1),strans3(:,2),strans3(:,3),'k')

%hohmann done

r5 = strans3(end,1:3);
vcir = (V2/norm(strans3(end,4:6)))*strans3(end,4:6);
[~, strans41] = ode45(@odefun,[0 T1],[r5 vcir],options, mu); %this doesn't happen, used to find orbit intersection

ll = 1;
kk = 0;
% diffe = zeros(1,(length(statecir)*length(state1)));
for ii = 1:length(strans41)
    for jj = 1:length(state1)
        diff = norm(strans41(ii,1:3) - state1(jj,1:3));
%         diffe(ll) = diff;
%         ll =ll + 1;
       if diff < 25
          index = [ii,jj];
%           disp(strans41(ii,1:3))
%             disp(state1(jj,1:3))
            kk = 1;
            break
       end
    end
    if kk == 1
        break
    end
end

% [min,index] = min(diffe);

r = state1(jj,1:3);
v1 = strans41(ii,4:6);
v2 = state1(jj,4:6);
deltav4 = norm(v2-v1);

angle = acos(dot(r,r5)/(norm(r)*norm(r5)));
% Tcir = Tfinder(mu,norm(r)); %approximately
Tquartercir =  Tcir*((angle/(2*pi)));


t9 = t8 + Tquartercir;
[endstates8,states8,times8] = propogate(endstates7,t8,t9,'of','n','n',1); %cirular orbit transfer
[ttrans4, strans4] = ode45(@odefun,[t8 t9],[r5 vcir],options, mu);
plot3(strans4(:,1),strans4(:,2),strans4(:,3),'k')

endstate1 = endstates8{1};
rtitan = endstate1(1:3);
rsat = strans4(end,1:3);

angle = acos(dot(rtitan,rsat)/(norm(rtitan)*norm(rsat)));
ratio = angle/(2*pi);
Tphase = T1 + T1*ratio;
t10 = t9 + Tphase;

aphase = afinder(Tphase,mu);
vmagphase = vfinder1(mu,norm(rsat),aphase);
vphase = (vmagphase/norm(v2))*v2;
deltavphase2 = 2*norm(vphase - v2);

[endstates9,states9,times9] = propogate(endstates8,t9,t10,'of','n','n',1);
[ttrans5, strans5] = ode45(@odefun,[t9 t10],[rsat vphase],options, mu);
plot3(strans5(:,1),strans5(:,2),strans5(:,3),'k')

t11 = t10 + 5*T1;
[endstates10,states10,times10] = propogate(endstates9,t10,t11,'of','n','n',1);

legend ('Titan','Pegasus','Falcon','Ayame','Mission Trajectory')

TotalDeltaV = deltavlam1+deltav2+deltav3+deltav4+deltavphase1+deltavphase2
TotalTime = t11/24/60/60
%


state1 = [states1{1};states2{1};states3{1};states4{1};states5{1};states6{1};states7{1};states8{1};states9{1};states10{1}];
state2 = [states1{2};states2{2};states3{2};states4{2};states5{2};states6{2};states7{2};states8{2};states9{2};states10{2}];
state3 = [states1{3};states2{3};states3{3};states4{3};states5{3};states6{3};states7{3};states8{3};states9{3};states10{3}];
state4 = [states1{4};states2{4};states3{4};states4{4};states5{4};states6{4};states7{4};states8{4};states9{4};states10{4}];

missionstate = [states1{2};states2{2};strans1;states4{3};strans2;states6{4};strans3;strans4;strans5;states10{1}];

time1 = [times1{1};times2{1};times3{1};times4{1};times5{1};times6{1};times7{1};times8{1};times9{1};times10{1}];
time2 = [times1{2};times2{2};times3{2};times4{2};times5{2};times6{2};times7{2};times8{2};times9{2};times10{2}];
time3 = [times1{3};times2{3};times3{3};times4{3};times5{3};times6{3};times7{3};times8{3};times9{3};times10{3}];
time4 = [times1{4};times2{4};times3{4};times4{4};times5{4};times6{4};times7{4};times8{4};times9{4};times10{4}];

missiontime = [times1{2};times2{2};ttrans1;times4{3};ttrans2;times6{4};ttrans3;ttrans4;ttrans5;times10{1}];

figure
plot3(missionstate(:,1),missionstate(:,2),missionstate(:,3))
hold on
plot3(state1(:,1),state1(:,2),state1(:,3))
plot3(state2(:,1),state2(:,2),state2(:,3))
plot3(state3(:,1),state3(:,2),state3(:,3))
plot3(state4(:,1),state4(:,2),state4(:,3))
axis equal

x0 = missionstate(1,1); y0 = missionstate(1,2); z0 = missionstate(1,3);
x1 = state1(1,1); y1 = state1(1,2); z1 = state1(1,3);
x2 = state2(1,1); y2 = state2(1,2); z2 = state2(1,3);
x3 = state3(1,1); y3 = state3(1,2); z3 = state3(1,3);
x4 = state4(1,1); y4 = state4(1,2); z4 = state4(1,3);


m0 = plot3(x0,y0,z0,'bs');
m1 = plot3(x1,y1,z1,'k*');
m2 = plot3(x2,y2,z2,'k*');
m3 = plot3(x3,y3,z3,'k*');
m4 = plot3(x4,y4,z4,'k*');
legend ('Mission Trajectory','Titan','Pegasus','Falcon','Ayame')
zz=1;aa=1;bb=1;cc=1;dd=1;

delta = 1;
    nm = length(missiontime);
    n1 = length(time1);
    n2 = length(time2);
    n3 = length(time3);
    n4 = length(time4);

input('Hit Enter to start aninmation: ', 's');
for ii = 1:delta:t11
    

    if missiontime(zz) < ii
        x0 = missionstate(zz,1); y0 = missionstate(zz,2); z0 = missionstate(zz,3);
        m0.XData = x0; m0.YData = y0; m0.ZData = z0;
        zz = zz + delta;
        
    end

    

    if time1(aa) < ii
        x1 = state1(aa,1); y1 = state1(aa,2); z1 = state1(aa,3);
        m1.XData = x1; m1.YData = y1; m1.ZData = z1;
        aa = aa + delta;
    end

    
    
    if time2(bb) < ii
        x2 = state2(bb,1); y2 = state2(bb,2); z2 = state2(bb,3);
        m2.XData = x2; m2.YData = y2; m2.ZData = z2;
        bb = bb + delta;
    end
    
    
    
    if time3(cc) < ii
        x3 = state3(cc,1); y3 = state3(cc,2); z3 = state3(cc,3);
        m3.XData = x3; m3.YData = y3; m3.ZData = z3;
        cc = cc + 1; 
    end
    
    
    
    if time4(dd) < ii
        x4 = state4(dd,1); y4 = state4(dd,2); z4 = state4(dd,3);
        m4.XData = x4; m4.YData = y4; m4.ZData = z4;
        dd = dd + delta;
    end
    
pause(0.000001)
   
end




function [h] = hfinder3(r, ecc,theta, mu)
h2 = r*(1 + (ecc*cos(theta)))*mu;
h = sqrt(h2);
end

function [a] = afinder(T,mu)
top = sqrt(mu)*T;
bot = 2*pi;
a = (top/bot)^(2/3);

end

function [v] = vfinder1(mu,r,a)
v = sqrt(2*((mu/r) - (mu/(2*a))));
end

function dydt = odefun(t,state, mu)
%credit to Dr. A's Thursday lecture
x = state(1);
y = state(2);
z = state(3);
dx = state(4);
dy = state(5);
dz = state(6);

rmag = norm([ x y z]); %km

ddx = (-mu*x/(rmag^3));
ddy = (-mu*y/(rmag^3));
ddz = (-mu*z/(rmag^3));

dydt = [dx; dy; dz; ddx; ddy; ddz];
end

function [v1,v2] = lamberts(mu,r1,r2,dt,z0,prograde)

r1_ = norm(r1);
r2_ = norm(r2);

[~,y,A] = zyAfinder(r1,r2,dt,z0,prograde);

f = 1 - (y/r1_);
g = A*sqrt(y/mu);
dg = 1 - (y/r2_);

v1 = (1/g)*(r2 - (f*r1));
v2 = (1/g)*((dg*r2) - r1);
end


function [z,y,A] = zyAfinder(r1,r2,dt,z0,prograde)
r1_ = norm(r1);
r2_ = norm(r2);
mu = 398600;

dtheta = dthetafinder(r1,r2,prograde);
top = r1_*r2_;
bot = 1 - cos(dtheta);
A = sin(dtheta) * sqrt(top/bot);

z = z0;

[S,C] = SCtrig(z);

y = r1_ + r2_ + A*((z*S)-1)/(sqrt(C));

F = ( ( (y/C)^(3/2) )* S) + (A*sqrt(y)) - (dt*sqrt(mu));

if z ~= 0
    one = (y/C)^(3/2);
    two = (1/(2*z))*(C - (3*S)/(2*C));
    three = ((3/4)* ((S^2)/C));
    four = (A/8);
    five = 3*(S/C)*sqrt(y);
    six = A*sqrt(C/y);
    
    dF = one*(two+three)+four*(five+six);
else
    one = (sqrt(2)/40)*(y^(3/2));
    two = A/8;
    three = sqrt(y);
    four = A*sqrt(1/(2*y));
    
    dF = one + (two*(three+four));
end
z0 = z;
z = z - (F/dF);

err = z -z0;

while abs(err) > 10^-8

[S,C] = SCtrig(z);

y = r1_ + r2_ + A*((z*S)-1)/(sqrt(C));

F = ( (y/C)^(3/2) )* S + (A*sqrt(y)) - (sqrt(mu)*dt);

if z ~= 0
    one = (y/C)^(3/2);
    two = (1/(2*z))*(C - (3*S)/(2*C));
    three = ((3/4)* ((S^2)/C));
    four = (A/8);
    five = 3*(S/C)*sqrt(y);
    six = A*sqrt(C/y);
    
    dF = one*(two+three)+four*(five+six);
else
    one = (sqrt(2)/40)*(y^(3/2));
    two = A/8;
    three = sqrt(y);
    four = A*sqrt(1/(2*y));
    
    dF = one + (two*(three+four));
end
z0 = z;
z = z - (F/dF);

err = z -z0;
    
end

y = r1_ + r2_ + A*((z*S)-1)/(sqrt(C));

end

function [dtheta] = dthetafinder(r1,r2,prograde)
r1_ = norm(r1);
r2_ = norm(r2);
rcr= cross(r1,r2);

angle = dot(r1,r2)/(r1_*r2_);

if prograde == 'n'
if rcr(3) >= 0
    dtheta = acos(angle);
else
    dtheta = (2*pi) - acos(angle);
end

else
    
if rcr(3) <= 0
    dtheta = acos(angle);
else
    dtheta = (2*pi) - acos(angle);
end

end

end

function [T] = Tfinder(mu,a)
T = (2*pi*(a^(3/2)))/sqrt(mu);
end

function [a, e, i, RA_Omega, AP_omega, nu, h_mag,rp,ra] = find_coes_degrees(R, V)

%set the relevant reference vectors
K = [0,0,1];
I = [1,0,0];
%to make things cleaner, calculate the magnitudes and dot product
V_mag = norm(V);
R_mag = norm(R);
RoV = dot(R,V);
%set mu
mu = 3.986 * 10^5;
%find epsilon
epsi = (0.5)*((V_mag)^2) - (mu/R_mag);
%find a
a = (-1)*mu/(2*epsi);
%find e vector and e
e_vector = (1/mu)*(((((V_mag)^2) -(mu/R_mag))*(R)) - ((RoV)*(V)));
e = norm(e_vector);
e_hat = e_vector/e;
%find h and its magnitude
h = cross(R,V);
h_mag = norm(h);
%create this important dot product
Koh = dot(K,h);
%find i
i = acosd ((Koh)/(h_mag));
%fine line of nodes and n's magnitude and the dot product of n and i hat
n = cross(K,h);
n_mag = norm(n);
Ion = dot(I,n);
%calculate the big omega
RA_Omega = acosd((Ion)/(n_mag));
if n(2) < 0 %quadrant check
RA_Omega = 360 - RA_Omega;
end
%dot product
noe = dot(n,e_vector);
%calculate little omega
AP_omega = acosd((noe)/((n_mag)*(e)));
if e_vector(3) < 0 % quadrant check
AP_omega = 360 - AP_omega;
end
%useful dot product
eoR = dot(e_vector,R);
%calcult v
nu = acosd((eoR)/((e)*(R_mag)));
if RoV < 0 %quadrant check
nu = 360 - nu;
end

rp_mag = ((h_mag^2)/mu)*(1/(1+e));
rp = rp_mag*e_hat;

ra_mag = ((h_mag^2)/mu)*(1/(1-e));
ra = -ra_mag*e_hat;
end

function [S,C] = SCtrig(z)

if z > 0
    S = (sqrt(z) - sin(sqrt(z)))/(sqrt(z^3));
elseif z < 0 
    S = (sinh(sqrt(-z)) - sqrt(-z))/(sqrt(-z)^3);
elseif z == 0 
    S = (1/6);
end

if z > 0
    C = (1 - cos(sqrt(z)))/z;
elseif z < 0 
    C = (cosh(sqrt(-z)) - 1)/(-z);
elseif z == 0 
    C = 0.5;
end
end

function [endstates,states,times] = propogate(initstates,t1,t2,plot,initpoints,finpoints,n)
mu =398600;

initstate1 = initstates{1};
initstate2 = initstates{2};
initstate3 = initstates{3};
initstate4 = initstates{4};
% ODE CAll
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

[tnew1, state1] = ode45(@odefun,[t1 t2],initstate1,options, mu);
[tnew2, state2] = ode45(@odefun,[t1 t2],initstate2,options, mu);
[tnew3, state3] = ode45(@odefun,[t1 t2],initstate3,options, mu);
[tnew4, state4] = ode45(@odefun,[t1 t2],initstate4,options, mu);

states = {state1,state2,state3,state4};
times = {tnew1,tnew2,tnew3,tnew4};

if plot == 'on'

figure(n)
plot3(state1(:,1),state1(:,2),state1(:,3))
hold on
plot3(state2(:,1),state2(:,2),state2(:,3))
plot3(state3(:,1),state3(:,2),state3(:,3))
plot3(state4(:,1),state4(:,2),state4(:,3))

end

x1 = state1(1,1); y1 = state1(1,2); z1 = state1(1,3);
x2 = state2(1,1); y2 = state2(1,2); z2 = state2(1,3);
x3 = state3(1,1); y3 = state3(1,2); z3 = state3(1,3);
x4 = state4(1,1); y4 = state4(1,2); z4 = state4(1,3);

x10= state1(end,1); y10= state1(end,2); z10= state1(end,3);
x20= state2(end,1); y20= state2(end,2); z20= state2(end,3);
x30= state3(end,1); y30= state3(end,2); z30= state3(end,3);
x40= state4(end,1); y40= state4(end,2); z40= state4(end,3);

if initpoints == 'y'
plot3(x1,y1,z1,'k*');
plot3(x2,y2,z2,'k*');
plot3(x3,y3,z3,'k*');
plot3(x4,y4,z4,'k*');
end

if finpoints == 'y'
plot3(x10,y10,z10,'b*');
plot3(x20,y20,z20,'b*');
plot3(x30,y30,z30,'b*');
plot3(x40,y40,z40,'b*');
end

xlabel('x [km]')
ylabel('y [km]')
zlabel('z [km]')
axis equal
legend ('Titan','Pegasus','Falcon','Ayame')

endstates = {state1(end,1:6),state2(end,1:6),state3(end,1:6),state4(end,1:6)};
end
