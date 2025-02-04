clear;
close;

rng(33);

addpath('sk_function')

Zpp = exp(1i*pi*(0:0.01:2));
dthp = pi/200; thp = dthp/2:dthp:2*pi-dthp/2;
Zp = exp(1i*thp);

% confmap


%% domainの設定
logn = @(z,a) log((z - a).*exp(-1i*angle(-a))) + 1i*angle(-a);
rho = 0.3;
Pz = @(z) P(z,rho); Phz = @(z) Ph(z,rho);
Kz = @(z) K(z,rho); Khz = @(z) Kh(z,rho);
Lz = @(z) L(z,rho); Lhz = @(z) Lh(z,rho);
a = 1.3; A = exp(1i*pi/4);
g = @(z) A*z.*Pz(z/a*sqrt(rho)).*Pz(1i*z/a*sqrt(rho)).*Pz(-z/a*sqrt(rho)).*Pz(-1i*z/a*sqrt(rho)) ...
       ./Pz(z/a)./Pz(-z/a)./Pz(1i*z/a)./Pz(-1i*z/a);
dgdz = @(z) g(z).*(1 + Kz(z/a*sqrt(rho)) + Kz(z/a*sqrt(rho)*1i) + Kz(-z/a*sqrt(rho)) + Kz(-1i*z/a*sqrt(rho)) ...
            -  Kz(z/a) - Kz(-z/a) -  Kz(z/a*1i) - Kz(-z/a*1i))./z;
dg2dz = @(z) (dgdz(z) - g(z)./z).*((1 + Kz(z/a*sqrt(rho)) + Kz(-z/a*sqrt(rho)) + Kz(z/a*sqrt(rho)*1i) + Kz(-z/a*sqrt(rho)*1i) ...
        -  Kz(z/a) - Kz(-z/a) -  Kz(z/a*1i) - Kz(-z/a*1i))./z) + ...
        g(z).*(Lz(z/a*sqrt(rho)) + Lz(-z/a*sqrt(rho)) + Lz(z/a*sqrt(rho)*1i) + Lz(-z/a*sqrt(rho)*1i) - Lz(z/a) - Lz(-z/a) - Lz(z/a*1i) - Lz(-z/a*1i))./z.^2;

hz = @(z) - log(dgdz(z));
dhdz = @(z) - dg2dz(z)./dgdz(z);

logp = @(z,a) -log(-a) + logn(z,a) + log(Phz(z/a));
dlogpdz = @(z,a) 1./(z - a) + 1./z.*Khz(z./a);

dh = 1e-3;
%zn = [+0.8, -0.8  - 0.1i, 0.8i, -0.8i];
zn = [0.8+0.2i, -0.5+0.1i];
qn = [-0.5, 0.5];
%c0 = sum(qn.*angle(zn)); have to modify because we are considering the map
c0 = imag(sum(qn.*log(zn)));
c1 = - c0;
Period = 0;
Qn = c1 + pi*Period;
fhat = @(z) fhat_doubly(z,zn,qn,rho);
f = @(z) fhat(z) - (1i*Qn./log(rho) + 1)*log(z) - 1i*c0 + pi*1i/2 + hz(z); %
fmod = @(z) fhat(z) - (1i*Qn./log(rho) + 1)*log(z) - 1i*c0 + pi*1i/2;
dfhatdz = @(z) dfhatdz_doubly(z,zn,qn,rho);
dfdz = @(z) dfhatdz(z) - (1i*Qn/log(rho) + 1)./z + dhdz(z); % OK

ep = 1e-3; dx = 1e-5;
dth = pi/8000;
Zp = exp(1i*(dth/2:dth:2*pi-dth/2));

calc_en = 1;
if calc_en == 1
%% integral function
intf = @(z) conj(f(z)).*dfdz(z);
lint_C0 = sum(intf(Zp).*Zp*1i)*dth/2i;
lint_Crho = -sum(intf(rho*Zp).*Zp*1i*rho)*dth/2i;

% integral from outer boundary to lz2
for k = 1:length(zn)
    epn(k) = ep./abs(dgdz(zn(k)));
    lzk = exp(1i*angle(zn(k))) * (abs(zn(k))+epn(k):dx:1);
    lzlist{k} = lzk;
    %lint_lk(k) = trapz(lzk,(intf(lzk*exp(1i*dx))))/2i - trapz(lzk,(intf(lzk./exp(1i*dx))))/2i;
    lint_lk(k) =  - 1i*2*pi*qn(k)*trapz(lzk,(dfdz(lzk)))/2i;
    lint_zk(k) =  - sum(intf(zn(k)+epn(k)*Zp).*epn(k).*Zp*1i)*dth/2i;
end

l1rho = -1:dx:-rho;
lint_1rho = trapz(l1rho,(dfdz(l1rho)))/2i* (2*pi*1i*(-1i*Qn/log(rho)+1));
Fint = sum(lint_lk) + sum(lint_zk) + lint_C0 + lint_Crho + lint_1rho;
disp("Integral");
disp(Fint);
end

%% robin function and ,,,?

Fep = - sum(qn.^2) * log(ep); 
Fg = sum((qn.^2 - 2*qn).*log(abs(dgdz(zn))));
Fcap = -(Qn.^2/log(rho).^2 + 1)*log(rho);
%rbk = @(z,k) fmod(z) + qn(k)*logn(z,zn(k));
for k = 1:length(zn)
    rbkt = rbk_doubly(zn,qn,rho,k,Period); 
    %rbkt = rbk(zn(k)+1i*1e-10,k);
    Fr_mod(k) = qn(k)*real(rbkt);
    %disp(rbkt)
    %disp(rbk(zn(k)+1i*1e-8,k));
end
addS = 1/2i * sum(1i*dhdz(Zp).*conj(hz(Zp)).*Zp) * dth - 1/2i * sum(1i*dhdz(rho*Zp).*conj(hz(rho*Zp)).*Zp*rho) * dth; 
F3 = - 2*sum(qn.*real(log(zn)));
total_energy = 2*pi*(Fep + Fg + Fcap + sum(Fr_mod) + F3) + addS;
disp(total_energy)
energy_list = total_energy;
disp(Fint - total_energy);




figure()
bv0 = exp(1i*pi/4)*bv0;
bv1 = exp(1i*pi/4)*bv1;
fill(real(bv0),imag(bv0),[0.9,0.9,0.9]);
hold on;
fill(real(bv1),imag(bv1),[1,1,1]);
plot(bv1,'k-','LineWidth',2.0);
plot(bv0,'k-','LineWidth',2.0);
%shading interp;
axis square;
axis off;


function fhat = fhat_doubly(z,zn,qn,rho)
logn = @(z,a) log((z - a).*exp(-1i*angle(-a))) + 1i*angle(-a);
Pz = @(z) P(z,rho); Phz = @(z) Ph(z,rho);
logp = @(z,a) -log(-a) + logn(z,a) + log(Phz(z/a));
N = length(zn);
fhat = 0;
for k = 1:N
    fhat = fhat - qn(k)*(logp(z,zn(k)) + log(Pz(z.*conj(zn(k)))));
end
end

function rbkt = rbk_doubly(zn,qn,rho,k,Period)
c0 = imag(sum(qn.*log(zn)));
c1 = - c0;
Qn = c1 + pi*Period;
Pz = @(z) P(z,rho); Phz = @(z) Ph(z,rho);
N = length(zn);
rbkt = -qn(k)*(log(Phz(1)) + log(Pz(zn(k).*conj(zn(k)))) - log(-zn(k))) ;
rbkt2 = -1i*c0 - (1+1i*Qn/log(rho))*log(zn(k)) + 1i*pi/2;
for l = 1:N
    if l ~= k 
        rbkt = rbkt - qn(l)*(log(Pz(zn(k)/zn(l))) + log(Pz(zn(k).*conj(zn(l)))));
    end
end
rbkt = rbkt;
end


function dfhatdz = dfhatdz_doubly(z,zn,qn,rho)

Kz = @(z) K(z,rho); Khz = @(z) Kh(z,rho);
dlogpdz = @(z,a) 1./(z - a) + 1./z.*Khz(z./a);
N = length(zn);
%dfhatdz = @(z) -1i*qn(1)*(dlogpdz(z,zn(1)) + dlogpdz(z,1./conj(zn(1)))) + ...
%               -1i*qn(2)*(dlogpdz(z,zn(2)) + dlogpdz(z,1./conj(zn(2))));
dfhatdz = 0;
for k = 1:N
    dfhatdz = dfhatdz - qn(k)*(dlogpdz(z,zn(k)) + dlogpdz(z,1./conj(zn(k))));
end
end


