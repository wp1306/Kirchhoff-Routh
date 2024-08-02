%% conduct Newton method 
% for theta=pi/2
clear;
close;

rng(33);

%addpath('./../useful_function');

plotf1 = 0;
plotf2 = 0;
calc1 = 0;

Zpp = exp(1i*pi*(0:0.01:2));
alpha_val = 0.5;

%% figure property
set(0,'defaultAxesFontSize',15);
set(0,'defaultAxesFontName','Arial')
set(0,'defaultlegendFontName','Arial')
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',1.2);

%zlist = 0.85:0.01:0.95;
z0 = 0.45+0.4i; z1 = -0.3i; z2 = -0.2; z3 = -0.3 - 0.3i;
zn = [z0, z1,z2]; 
qn = [1/2, 3/2,-1];
M = length(zn);
%f = @(z) fd(z,qn,zn);
dfdz = @(z) dfdz_func(z,qn,zn);
logn = @(z,a) log((z - a).*exp(-1i*angle(-a))) + 1i*angle(-a);
epp = 1e-6;

%% Region: Quadrature domain
a = 1.5; A = (a^6 - 1)/(a^3-a^2 + 1);
gz = @(z) A*z./(a.^3 - z.^3);
dgdz = @(z) A*(a^3 + 2*z.^3)./(a.^3 - z.^3).^2;
dg2dg = @(z) 6*z.^2./(a^3+2*z.^3) + 6*z.^2./(a.^3 - z.^3);
wn = gz(zn);
hz  = @(z) - 1i*log(dgdz(z)); % 
dhz = @(z) -1i*dg2dg(z);

tw = @(w) w.*(sqrt(3)*sqrt(27*a^6+4*A^3./w.^3) + 9*a^3).^(1/3);
ginvw1 = @(w) tw(w)./(2.^(1/3)*3^(2/3)*w) - (2/3)^(1/3)*A./tw(w);
ginvw2 = @(w) (1+1i*sqrt(3))*A./(2.^(2/3)*3.^(1/3)*tw(w)) - (1-1i*sqrt(3))*tw(w)./(2*2^(1/3)*3^(2/3)*w);
ginvw3 = @(w) (1-1i*sqrt(3))*A./(2.^(2/3)*3.^(1/3)*tw(w)) - (1+1i*sqrt(3))*tw(w)./(2*2^(1/3)*3^(2/3)*w);
ginvw = @(w) ginvw1(w).*(abs(angle(w))<pi/3) + ginvw2(w).*(abs(angle(w)-2*pi/3)<pi/3) + ginvw3(w).*(abs(angle(w)+2*pi/3)<pi/3);

f = @(z) fd(z,qn,zn);
fDw = @(w) f(ginvw(w));
hw = @(w) hz(ginvw(w));

fmodz = @(z) f(z) + hz(z);
fmodw = @(w) f(ginvw(w)) + hz(ginvw(w));
dfmoddz = @(z) dfdz(z) + dhz(z);
rbmod = @(z,n) fmodz(z) + 1i*qn(n)*logn(z,zn(n));
rb = @(z,n) f(z) + 1i*qn(n)*logn(z,zn(n));
dfmoddw = @(w) dfdz(ginvw(w)) + dhz(ginvw(w));


phi = @(w) real(fDw(w)) + real(hw(w));
phih = @(w) phi(w) - pi*floor((phi(w)+pi/2)/pi);
phiz = @(z) real(f(z)) + real(hz(z));
phihz = @(z) phiz(z) - pi*floor((phiz(z)+pi/2)/pi);

ep = 1e-3;
dth = pi/100000;
th = dth/2:dth:2*pi-dth/2; Zp = exp(1i*th);
% calculate energy

Fd = 0; Fep = 0;
for j = 1:length(zn)
    Fep = Fep - 2*pi*qn(j)^2*log(ep);
    Fd = Fd + (qn(j).^2 - 2*qn(j))*log(abs(dgdz(zn(j)))) + qn(j).*imag(rb(zn(j)+1i*1e-10,j));
end

addS = 1/2i * sum(1i*dhz(Zp).*conj(hz(Zp)).*Zp) * dth; 
Fg = addS + 4*pi*log(abs(dgdz(0))); 
total_energy = Fep + 2*pi*Fd + Fg;

%% line integral around zk;
fint = @(z) dfmoddz(z).*conj(fmodz(z));
epr = ep./abs(dgdz(zn));  % radius;

% around the circle
for j = 1:length(zn)
    lep(j) = -sum(1/2i * fint(zn(j)+epr(j)*Zp).*Zp*epr(j)*1i)*dth;
end
lc0 = sum(1/2i * fint(Zp).*Zp.*1i)*dth;
% line integral; from zn -> 1
dx = 1e-8; smallep = 1e-8;
for j = 1:M
edx = epr(j);
thg(j) = mod(angle(zn(j)),2*pi);
l1 = flip(exp(1i*thg(j))*(edx:dx:1-abs(zn(j))) + zn(j));
lt{j} = l1;
chk(j,1) = trapz(l1,dfmoddz(l1+1i*smallep))*2*pi*qn(j)/(2i); % 2*pi*q1 がぷらす．
end





% dd = 0.001;
% [Xc,Yc] = meshgrid(-1:dd:1);
% maskc = (Xc.^2 + Yc.^2 < 1) + 0;
% Fg2 = sum(abs(dg2dg(Xc+1i*Yc)).^2.*maskc,"all")*dd.^2; % this is equal to
% addS













plotf = 1;
if plotf == 1
boundary = gz(Zpp);
dd = 0.02; Dx = 2.2;
[X,Y] = meshgrid(-Dx+dd/2:dd:Dx+dd/2);
mask = inpolygon(X,Y,real(boundary),imag(boundary)) + 0;

h = figure();
mask(mask==0) = NaN;
phihXY = real(phih(X+1i*Y)).*mask;
p = pcolor(X,Y,phihXY);
colormap(h,'hsv');
clim([-pi/2,pi/2]);
colorbar('Ticks',[-pi/2+1e-3,-pi/4,0,pi/4,pi/2-1e-3],'TickLabels',{'$-\pi/2$','$-\pi/4$','$0$','$\pi/4$','$\pi/2$'},'TickLabelInterpreter','latex','FontSize',20);
alpha(p,alpha_val);
hold on;
axis equal;
shading flat;
plot(gz(Zpp),"k-","LineWidth",3.0);
ll = 0.2; di = 0.2;
[Xp,Yp] = meshgrid(-Dx+di/2:di:Dx-di/2);
for xi = 1:length(Xp)
    for yi = 1:length(Yp)
        wp = Xp(xi,yi)+1i*Yp(xi,yi);
        if inpolygon(real(wp),imag(wp),real(boundary),imag(boundary))
            angp = phih(wp);
            vec = ll*exp(1i*angp)+1i*0.00001;
            plot([wp - vec/2, wp + vec/2],"k-","LineWidth",0.5);
        end
    end
end
plot(gz(zn)+1i*1e-10,"ro","MarkerSize",10.0);
plot(gz(zn(1)) + ep*Zpp,"b-");


boundary = Zpp;
dd = 0.02; Dx = 1.0;
[X,Y] = meshgrid(-Dx+dd/2:dd:Dx+dd/2);
mask = inpolygon(X,Y,real(boundary),imag(boundary)) + 0;

h = figure();
mask(mask==0) = NaN;
phihXY = real(phiz(X+1i*Y)).*mask;
p = pcolor(X,Y,phihXY);
colormap(h,'turbo');
%clim([-pi/2,pi/2]);
colorbar('Ticks',[-pi/2+1e-3,-pi/4,0,pi/4,pi/2-1e-3],'TickLabels',{'$-\pi/2$','$-\pi/4$','$0$','$\pi/4$','$\pi/2$'},'TickLabelInterpreter','latex','FontSize',20);
alpha(p,alpha_val);
hold on;
axis equal;
shading flat;
plot(Zp,"k-","LineWidth",3.0);
ll = 0.2; di = 0.2;
[Xp,Yp] = meshgrid(-Dx+di/2:di:Dx-di/2);
for xi = 1:length(Xp)
    for yi = 1:length(Yp)
        wp = Xp(xi,yi)+1i*Yp(xi,yi);
        if inpolygon(real(wp),imag(wp),real(boundary),imag(boundary))
            angp = phiz(wp);
            vec = ll*exp(1i*angp)+1i*0.00001;
            plot([wp - vec/2, wp + vec/2],"k-","LineWidth",0.5);
        end
    end
end
plot(zn+1i*1e-10,"ro","MarkerSize",10.0);
plot(ginvw(gz(zn(1)) + ep*Zp),"b-");  % the radius goes to ep/abs(dgdz(zn(1)));
for j = 1:length(zn)
plot(lt{j},"k-");
end

end


function f = fd(z,qn,zn)
M = length(qn);
f = - pi/2;
logn = @(z,a) log((z - a).*exp(-1i*angle(-a))) + 1i*angle(-a);
for j = 1:M
    %f = f - 1i*qn(j)*(log(z-zn(j)) + log(1 - conj(zn(j))*z));
    %f = f - 1i*qn(j)*(log(z) + log(1 - zn(j)./z) + log(1 - z.*conj(zn(j))));
    f = f - 1i*qn(j)*(logn(z,zn(j)) + log(1 - conj(zn(j))*z));
    %f = f - 1i*qn(j)*(log(z) + logn(1,zn(j)./z) + log(1-conj(zn(j)).*z));
end
end

function df_dz = dfdz_func(z,qn,zn)
M = length(qn);
df_dz = 0;
for j = 1:M
    df_dz = df_dz - 1i*qn(j)./(z - zn(j)) + 1i*qn(j).*conj(zn(j))./(1 - conj(zn(j))*z);
end
end




