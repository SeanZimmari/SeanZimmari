%%Plot delle funzioni di corrente e potenziale di un cilindro circolare

%Definzione dello stato iniziale


mi = 10; %Intensità della doppietta ( m^3/s )
U = 1; %Velocità della corrente indisturbata ( m/s )
a = sqrt(mi/(2*pi*U)); %in m
visc_din = 1.81e-5; %Viscosità dinamica dell'aria ( Pa s )
ro = 1.225; %Densità dell'aria ( kg/m^3 )
L = 2*a;


%Definizione del campo di lavoro

theta = 0:.02:2*pi; %Angolo ( rad ) 
r = 1:.02:3; %Distanza dalla sorgente ( m )
[R, Theta] = meshgrid(r,theta);
x = (R.*cos(Theta))./a;
y = (R.*sin(Theta))./a;

%Plot delle isolinee psi

%( 1 ) gamma = 0

gamma1 = 0; %Intensità del vortice irrotazionale ( m^2/s )
psi1 = U.*R.*sin(Theta).*((a^2)./(R.^2) - 1) - (gamma1/(2*pi)).*(log(R)-log(a));
phi1 = -(mi./(2*pi.*R)).*cos(Theta) - U.*x + ((gamma1.*Theta)./(2*pi));

% ( 2 ) 0 < gamma < 4*pi*a
gamma2 = 10;
psi2 = U.*R.*sin(Theta).*((a^2)./(R.^2) - 1) - (gamma2/(2*pi)).*(log(R)-log(a));
phi2 = -(mi./(2*pi.*R)).*cos(Theta) - U.*x + ((gamma2.*Theta)./(2*pi));

% ( 3 ) gamma = 4*pi*a
gamma3 = 4*pi*a;
psi3 = U.*R.*sin(Theta).*((a^2)./(R.^2) - 1) - (gamma3/(2*pi)).*(log(R)-log(a));
phi3 = -(mi./(2*pi.*R)).*cos(Theta) - U.*x + ((gamma3.*Theta)./(2*pi));

%( 4 ) gamma > 4*pi*a
gamma4 = 20;
psi4 = U.*R.*sin(Theta).*((a^2)./(R.^2) - 1) - (gamma4/(2*pi)).*(log(R)-log(a));
phi4 = -(mi./(2*pi.*R)).*cos(Theta) - U.*x + ((gamma4.*Theta)./(2*pi));

figure(1)
t = tiledlayout(2,2);
nexttile
[~,h] = contour(x,y,psi1,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,psi1,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \phi ( \gamma = 0 )','FontSize',18)
plot(0,0,'ko','LineWidth',2)

nexttile

[~,h] = contour(x,y,psi2,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,psi2,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \phi ( 0 < \gamma < 4\pia )','FontSize',18)
plot(0,0,'ko','LineWidth',2)

nexttile

[~,h] = contour(x,y,psi3,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,psi3,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \phi ( \gamma = 4\pia )','FontSize',18)
plot(0,0,'ko','LineWidth',2)

nexttile

[~,h] = contour(x,y,psi4,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,psi4,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \phi ( \gamma > 4\pia )','FontSize',18)
plot(0,0,'ko','LineWidth',2)

figure(2)
t = tiledlayout(2,2);
nexttile
[~,h] = contour(x,y,phi1,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,phi1,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \phi ( \gamma = 0 )','FontSize',18)
plot(0,0,'ko','LineWidth',2)

nexttile

[~,h] = contour(x,y,phi2,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,phi2,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \phi ( 0 < \gamma < 4\pia )','FontSize',18)
plot(0,0,'ko','LineWidth',2)

nexttile

[~,h] = contour(x,y,phi3,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,phi3,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \phi ( \gamma = 4\pia )','FontSize',18)
plot(0,0,'ko','LineWidth',2)

nexttile

[~,h] = contour(x,y,phi4,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,phi4,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \phi ( \gamma > 4\pia )','FontSize',18)
plot(0,0,'ko','LineWidth',2)

%Calcolo Cp
theta = 0:.02:pi;
r = a;
[R Theta] = meshgrid(r,theta);

%Considerando r = a, U = 1 e gamma = 0:
V_r = U.*cos(Theta).*((a^2)./(R.^2) - 1);
V_theta = U.*sin(Theta).*((a^2)./(R.^2) + 1) + gamma2/(2*pi.*R);
V = sqrt(V_r.^2 + V_theta.^2);
Cp0 = 1- (V./U0).^2;

%Linee di corrente del campo di velocità
clear x y
[x y] = meshgrid(-1:.01:1, -1:.01:1);
Theta = atan(y./x);
R = x./cos(Theta);
V_r = U.*cos(Theta).*((a^2)./(R.^2) - 1);
V_theta = U.*sin(Theta).*((a^2)./(R.^2) + 1) + gamma2./(2*pi.*R);
figure(3)
streamslice(x./a,y./a,V_r,V_theta,1.5)
hold on
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',15)
title('Linee di corrente del campo di velocità (V_r,V_\theta)','FontSize',18)

clear R Theta
R = a;
Theta = 0:.02:pi;
Angle = Theta*(180/3.14);
figure(4)
plot(Angle,Cp0,'m-','LineWidth',2)
xlabel('Angle ( deg )','FontSize',10)
ylabel('Cp','FontSize',10)
title('Cp per \gamma > 4\pia','FontSize',15)

%Calcolo la forza agente sul cilindro

%Considerando una superficie infinitesima ds = a*dtheta, la forza sarà pari
%al prodotto della pressione agente su ds per la superficie stessa
%Considerando il caso in cui gamma < 4*pi*a, la forza risultante si può
%scrivere come:
%dF = -(p - p_inf)*a*dtheta = -((ro*U*gamma)/pi)*sin(theta)
syms theta
dD = ((ro*U*gamma2)/pi)*sin(theta)*cos(theta);
dL = ((ro*U*gamma2)/pi)*sin(theta)*sin(theta);

D = int(dD,theta,0,2*pi)
L = int(dL,theta,0,2*pi)

%Si è dimostrato quindi numericamente il paradosso di D'Alembert, ovvero
%che un vortice irrotazionale induce la distribuzione di pressione ad
%essere non simmetrica tra la parte superiore e quella inferiore del cilindro, 
%pertanto questo porta ad un aumento di velocità nella parte superiore del
%cilindro ed una diminuzione di velocità nella parte inferiore.
%La forza che deriva da questa differenza di pressione viene denominata
%portanza e la presenza di una portanza nel caso di un cilindro circolare è
%detto effetto magnus
