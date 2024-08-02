%%%% conduct Newton method 
% for theta=pi/2
clear;
close;

rng(33);

%addpath('./../../useful_function');
%addpath('./../../newton_optimization');
%addpath('./comp_function');

plotf1 = 0;
    

%% figure property
set(0,'defaultAxesFontSize',17);
set(0,'defaultAxesFontName','Arial')
set(0,'defaultlegendFontName','Arial')
set(0,'defaulttextinterpreter','latex');
set(0,'defaultLineLineWidth',1.5);

thp = pi*(0:0.00025:2);
Zp = exp(1i*thp);
ZpL = exp(1i*pi*(1:0.01:2));
ZpU = exp(1i*pi*(0:0.005:1));
dth = pi/1000;
Zpp = exp(1i*(dth/2:dth:2*pi-dth/2));

z1 = 5^(-1/4); z2 = -z1;
f = @(z) -1i/2*(log(z-z1)+log(1-z.*conj(z1))) - 1i/2*(log(z-z2) + log(1-z.*conj(z2))) -pi/2;
ft = @(z) -1i/2*(log(z-z1) + log(1-z.*conj(z1))) - 1i/2*(log(1-z.*conj(z)));
phi = @(z) real(f(z));
phih = @(z) phi(z) - pi*floor((phi(z)+pi/2)/pi);
psi = @(z) imag(f(z));

boundary = Zp;
Dx = 1.0;
dd = 1e-2;
alpha_val = 0.4;
[X,Y] = meshgrid(-Dx+dd/2:dd:Dx+dd/2);
mask = inpolygon(X,Y,real(boundary),imag(boundary)) + 0;

h = figure();
mask(mask==0) = NaN;
phihXY = real(phih(X+1i*Y)).*mask;
psihXY = psi(X+1i*Y);
p = pcolor(X,Y,phihXY);
colormap(h,'hsv');
clim([-pi/2,pi/2]);
colorbar('Ticks',[-pi/2+1e-3,-pi/4,0,pi/4,pi/2-1e-3],'TickLabels',{'$-\pi/2$','$-\pi/4$','$0$','$\pi/4$','$\pi/2$'},'TickLabelInterpreter','latex','FontSize',20);
alpha(p,alpha_val);
hold on;
axis equal;
shading flat;
plot(Zp,"k-","LineWidth",2.0);
plot(z1+1i*1e-10,'ro','MarkerSize',10);
plot(z2+1i*1e-10,'ro','MarkerSize',10);
ll = 0.05; di = 0.05;
[Xp,Yp] = meshgrid(-Dx+di/2:di:Dx-di/2);
for xi = 1:length(Xp)
    for yi = 1:length(Yp)
        wp = Xp(xi,yi)+1i*Yp(xi,yi);
        %if inpolygon(real(wp),imag(wp),real(boundary),imag(boundary))
        if abs(wp) < 1
            angp = phih(wp);
            vec = ll*exp(1i*angp)+1i*0.00001;
            plot([wp - vec/2, wp + vec/2],"k-","LineWidth",1.0);
        end
    end
end
p3 = contour(X,Y,psihXY,30,"LineWidth",1.0);
%p4 = contour(X,Y,phihXY,30,"LineWidth",1.0);
alpha_val = 0.1;
alpha(p,alpha_val);

h = figure();
mask(mask==0) = NaN;
psiht = imag(ft(X+1i*Y)).*mask;
phiht = real(ft(X+1i*Y)).*mask;
%p = pcolor(X,Y,phiht);
p = contour(X,Y,psiht,300,'LineWidth',1.0);
hold on;
p2 = contour(X,Y,phiht,20,'LineWidth',0.5);
axis equal;
colormap(h,'parula');

