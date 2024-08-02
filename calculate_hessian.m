%% conduct Newton method 
% for theta=pi/2
clear;
close;

rng(33);

addpath('./../../useful_function');
addpath('./../../newton_optimization');
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

%% ------------------- 

% check the quadrature domain with different 


%M = length(zn);
%A = (a^6-1)/(a^3-a^2 + 1);%R = (a^4 - 1)/sqrt((1+a^4)/2);
%gz = @(z) A*z./(a.^3 - z.^3);
%f = @(z) fd(z,qn,zn);
%dgdz = @(z) A*(a^3 + 2*z.^3)./((a^3 - z.^3).^2);
%dg2dg_p = @(z) 12.*z.*(a^3 - z.^3)./(a^3+2*z.^3).^2 + 6*z.*(2*a^3+z.^3)./(a^3-z.^3).^2;
% dfdz = @(z) dfdz_func(z,qn,zn);
% hz  = @(z) - 1i*log(dgdz(z));
% wn = gz(zn);
% d = gz(1/a); 
% tw = @(w) w.*(sqrt(3)*sqrt(27*a^6+4*A^3./w.^3) + 9*a^3).^(1/3);
% ginvw1 = @(w) tw(w)./(2.^(1/3)*3^(2/3)*w) - (2/3)^(1/3)*A./tw(w);
% ginvw2 = @(w) (1+1i*sqrt(3))*A./(2.^(2/3)*3.^(1/3)*tw(w)) - (1-1i*sqrt(3))*tw(w)./(2*2^(1/3)*3^(2/3)*w);
% ginvw3 = @(w) (1-1i*sqrt(3))*A./(2.^(2/3)*3.^(1/3)*tw(w)) - (1+1i*sqrt(3))*tw(w)./(2*2^(1/3)*3^(2/3)*w);
% ginvw = @(w) ginvw1(w).*(abs(angle(w))<pi/3) + ginvw2(w).*(abs(angle(w)-2*pi/3)<pi/3) + ginvw3(w).*(abs(angle(w)+2*pi/3)<pi/3);
% 
% 
% %ginvw = @(w) -(1-1i*sqrt(3))*tw(w).^(1/3)./(2*2.^(1/3)*3^(2/3)*w) + (1-1i*sqrt(3))*A./(2^(2/3)*3^(1/3)*tw(w).^(1/3));
% 
% %fDw = @(w) f(ginvw(w));
% %hw = @(w) hz(ginvw(w));
% fmodz = @(z) f(z) + hz(z);
% fmodw = @(w) f(ginvw(w)) + hz(ginvw(w));
% phi = @(w) real(fmodw(w));
% phih = @(w) phi(w) - pi*floor((phi(w)+pi/2)/pi);


alist = [2^(1/3):0.01:4.5];
eig_two = zeros(length(alist),1);
eig_four = zeros(length(alist),1);
Delta_list = zeros(length(alist),1);

for j = 1:length(alist)
%% calculate Hessian.
a = alist(j);
A = (a^6-1)/(a^3-a^2 + 1);%R = (a^4 - 1)/sqrt((1+a^4)/2);
gz = @(z) A*z./(a.^3 - z.^3);
f = @(z) fd(z,qn,zn);
dgdz = @(z) A*(a^3 + 2*z.^3)./((a^3 - z.^3).^2);
dg2dz = @(z) 6*A*z.^2.*(z.^3 + 2*a^3)./(a^3 - z.^3).^3;
dg2dg = @(z) dg2dz(z)./dgdz(z);
hz = @(z) -1i*log(dg2dg(z));

dg2dg_p = @(z) 12.*z.*(a^3 - z.^3)./(a^3+2*z.^3).^2 + 6*z.*(2*a^3+z.^3)./(a^3-z.^3).^2;

J = [1, 1; 1i, -1i];
Delta = sqrt(3)*real(gz(1/a));
Delta_list(j,1) = Delta;

% triple
fg = @(s) [s(1)./(2*(1-s(1).^2)) - 1./(2*(s(1)-s(2))) + s(2)./(2*(1-s(1).*s(2))) - 3/4.*(6*s(1).^2./(a^3+2*s(1).^3) - 6*s(1).^2./(s(1).^3 - a^3)),...
s(2)./(2*(1-s(2).^2)) - 1./(2*(s(2)-s(1))) + s(1)./(2*(1-s(2).*s(1))) - 3/4.*(6*s(2).^2./(a^3+2*s(2).^3) - 6*s(2).^2./(s(2).^3 - a^3))];
options = optimset('TolFun',1e-7,'TolX',1e-15,'MaxFunEvals',3000,'MaxIter',3000,'Display','off');
sr = fminsearch(@(s) sum(abs(fg(s))),[0.3,0.8],options);
zn = [sr(1), sr(2)]; 
qn = [1/2, 1/2];
z1 = zn(1); z2 = zn(2);
H11 = [z1.^2/4/(1-z1^2)^2 + 1/4/(z1-z2)^2 + z2.^2/4/(1-z2*z1)^2 - 3/8*dg2dg_p(z1), 1/4/(1-z1^2)^2;
           1/4/(1-z1^2)^2, conj(z1.^2/4/(1-z1^2)^2 + 1/4/(z1-z2)^2 + z2.^2/4/(1-z2*z1)^2 - 3/8*dg2dg_p(z1))];
H22 = [z2.^2/4/(1-z2^2)^2 + 1/4/(z2-z1)^2 + z1.^2/4/(1-z1*z2)^2 - 3/8*dg2dg_p(z2), 1/4/(1-z2^2)^2;
           1/4/(1-z2^2)^2, conj(z2.^2/4/(1-z2^2)^2 + 1/4/(z2 - z1)^2 + z1.^2/4/(1-z1*z2)^2 - 3/8*dg2dg_p(z2))];
H12 = [1./4/(z1-z2)^2,1./4/(1-z1*z2)^2;1/4/(1-z1*z2)^2,1/4/(z1-z2)^2];
H21 = [1./4/(z2-z1)^2,1./4/(1-z2*z1)^2;1/4/(1-z2*z1)^2,1/4/(z2-z1)^2];
BigJ = [J,zeros(2,2);zeros(2,2),J];
Hessian = BigJ*[H11,H12;H21,H22]*transpose(BigJ);
eig_H1r = eig(Hessian);
eig_two(j,1) = min(eig_H1r);

%% eig two again
options = optimset('TolFun',1e-7,'TolX',1e-15,'MaxFunEvals',3000,'MaxIter',3000,'Display','off');
fg2 = @(rth) [-conj(rth(1)*exp(1i*rth(2)))./(1-rth(1).^2) + 1./(rth(1)*exp(1i*rth(2)) - rth(1)*exp(-1i*rth(2))) - rth(1)*exp(1i*rth(2))./(1-rth(1).^2*exp(2i*rth(2))) + 3/2*dg2dg(rth(1)*exp(1i*rth(2))), ...
          -conj(rth(1)*exp(-1i*rth(2)))./(1-rth(1).^2) + 1./(rth(1)*exp(-1i*rth(2)) - rth(1)*exp(1i*rth(2))) - rth(1)*exp(-1i*rth(2))./(1-rth(1).^2*exp(-2i*rth(2))) + 3/2*dg2dg(rth(1)*exp(-1i*rth(2)))];
rth = fminsearch(@(s) sum(abs(fg2(s))),[0.8930    2.0407],options);
zn = [rth(1)*exp(1i*rth(2)),rth(1)*exp(-1i*rth(2))];
qn = [1/2, 1/2];
z1 = zn(1); z2 = zn(2);
H11 = [conj(z1).^2/4/(1-z1*conj(z1))^2 + 1/4/(z1 - z2)^2 + conj(z2).^2/4/(1-conj(z2)*z1)^2 - 3/8*dg2dg_p(z1), 1/4/(1-z1*conj(z1))^2;
       1/4/(1-z1*conj(z1))^2, conj(conj(z1).^2/4/(1-z1*conj(z1))^2 + 1/4/(z1 - z2)^2 + conj(z2).^2/4/(1-conj(z2)*z1)^2 - 3/8*dg2dg_p(z1))];
H22 = [conj(z2).^2/4/(1-z2*conj(z2))^2 + 1/4/(z2 - z1)^2 + conj(z1).^2/4/(1-conj(z1)*z2)^2 - 3/8*dg2dg_p(z2), 1/4/(1-z2*conj(z2))^2;
       1/4/(1-z2*conj(z2))^2, conj(conj(z2).^2/4/(1-z2*conj(z2))^2 + 1/4/(z2 - z1)^2 + conj(z1).^2/4/(1-conj(z1)*z2)^2 - 3/8*dg2dg_p(z2))];

%H22 = [z2.^2/4/(1-z2^2)^2 + 1/4/(z2-z1)^2 + z1.^2/4/(1-z1*z2)^2 - 3/8*dg2dg_p(z2), 1/4/(1-z2^2)^2;
%       1/4/(1-z2^2)^2, conj(z2.^2/4/(1-z2^2)^2 + 1/4/(z2 - z1)^2 + z1.^2/4/(1-z1*z2)^2 - 3/8*dg2dg_p(z2))];
H12 = [1./4/(z1-z2)^2,1./4/(1-z1*conj(z2))^2;1/4/(1-conj(z1)*z2)^2,conj(1/4/(z1-z2)^2)];
%H21 = [1./4/(z2-z1)^2,1./4/(1-z2*z1)^2;1/4/(1-z2*z1)^2,1/4/(z2-z1)^2];
H21 = [1./4/(z2-z1)^2,1./4/(1-z2*conj(z1))^2;1/4/(1-conj(z2)*z1)^2,conj(1/4/(z2-z1)^2)];
BigJ = [J,zeros(2,2);zeros(2,2),J];
Hessian = BigJ*[H11,H12;H21,H22]*transpose(BigJ);
eig_H1 = eig(Hessian);
%disp(eig_H1);
eig_twor(j,1) = min(eig_H1);



%% calculate Hessian numerically??

% 
fg2 = @(s) s.^9 + 7*a^3*s.^6 - (3 - a.^6)*s.^3 - 6*a^3;
options = optimset('TolFun',1e-7,'TolX',1e-15,'MaxFunEvals',3000,'MaxIter',3000,'Display','off');
s = fminsearch(@(s) abs(fg2(s)),0.8,options);
omega = exp(2i*pi/3);
zn = [0,s, s*omega,s*omega^2]; 
qn = [-1/2,1/2,1/2, 1/2];
%  calculate Hessian for 4 topological defects
z0 = zn(1); z1 = zn(2); z2 = zn(3); z3 = zn(4); zk = [z1,z2,z3];
E = zeros(2,2);
BigJ = [J,E,E,E;E,J,E,E;E,E,J,E;E,E,E,J];
H00_1 = 5/8*dg2dg_p(z0) - 1/4 *sum(1./(z0-zk).^2 + (conj(zk)).^2); 
H00 = [H00_1,1/4;1/4, conj(H00_1)];
for m = 1:3
    zm = zk(m);
    H0m{m} = [1./4/(z0-zm).^2,-1/4/(1-conj(zm)*z0);conj(-1/4/(1-conj(zm)*z0)), conj(1./4/(z0-zm).^2)];
end
for m = 1:3
    zm = zk(m); zmt = zk(zk ~= zk(m));
    Hmm{m} = [-3/8*(dg2dg_p(zm)) + conj(zm).^2./4./(1-zm*conj(zm)).^2 - 1/4*(1./(zm-z0).^2) + 1/4*sum(1./(zmt-zm).^2 + conj(zmt).^2./(1 - zm*conj(zmt)).^2), ...
             1/4/(1-zm.*conj(zm)).^2;1/4/(1-zm.*conj(zm)).^2,  ...
             conj(-3/8*(dg2dg_p(zm)) + conj(zm).^2./4./(1-zm*conj(zm)).^2 - 1/4*(1./(zm-z0).^2) + 1/4*sum(1./(zmt-zm).^2 + conj(zmt).^2./(1 - zm*conj(zmt)).^2))];
end
H12 = [-1/4/(zk(1)-zk(2)).^2,1/4/(1-zk(1)*conj(zk(2))).^2;conj(1/4/(1-zk(1)*conj(zk(2))).^2),conj(-1/4/(zk(1)-zk(2)).^2)];
H13 = [-1/4/(zk(1)-zk(3)).^2,1/4/(1-zk(1)*conj(zk(3))).^2;conj(1/4/(1-zk(1)*conj(zk(3))).^2),conj(-1/4/(zk(1)-zk(3)).^2)];
H23 = [-1/4/(zk(2)-zk(3)).^2,1/4/(1-zk(2)*conj(zk(3))).^2;conj(1/4/(1-zk(2)*conj(zk(3))).^2),conj(-1/4/(zk(2)-zk(3)).^2)];
BigH = [H00,H0m{1},H0m{2},H0m{3};...
        transpose(H0m{1}),Hmm{1},H12,H13;...
        transpose(H0m{2}),transpose(H12),Hmm{2},H23;
        transpose(H0m{3}),transpose(H13),transpose(H23),Hmm{3}];
Hessian = BigJ*BigH*transpose(BigJ);
eig_H2 = eig(Hessian);

eig_four(j,1) = min(eig_H2);

end

figure()
c1 = plot(Delta_list,eig_two,"k-",'LineWidth',2.5);
hold on;
c2 = plot(Delta_list,eig_twor,"r-.",'MarkerSize',10.0,'LineWidth',2.5);
c3 = plot(Delta_list,eig_four,"b--",'LineWidth',2.5);
b1 = find(abs(eig_two)<1e-3);
b2 = find(abs(eig_four)<1e-3);
plot(Delta_list(b1(end)),0,"k.","MarkerSize",25);
plot(Delta_list(b2(end)),0,"b.","MarkerSize",25);
plot(Delta_list(1),eig_two(1),"ko","MarkerSize",10);
plot(Delta_list(1),eig_four(1),"bo","MarkerSize",10);
plot(Delta_list(1),eig_twor(1),"ro","MarkerSize",10);
plot(0.5:0.01:2,0*[0.5:0.01:2],"k:","LineWidth",1.5);
legend([c1,c2,c3],["two defects (i)","two defects (ii)","four defects"],"FontSize",20,"Location","northwest","Interpreter","latex");
title("The smallest eigenvalue of Hessian matrix","FontSize",20);
xlabel("$\Delta$","FontSize",20);
ylabel("$\lambda_{\rm min}$","FontSize",20);




function f = fd(z,qn,zn)
M = length(qn);
f = -pi/2;
for j = 1:M
    f = f - 1i*qn(j)*(log(z - zn(j)) + log(1 - conj(zn(j))*z));
end
end
