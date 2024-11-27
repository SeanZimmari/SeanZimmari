%%Plot delle funzioni di corrente e potenziale di una Prua di Fuhrmann

%Definizione dello stato iniziale

Q = 5; %Intensità della sorgente puntiforme
U = 1; %Intensità della velocità della corrente indisturbata V_inf = -Ui

%Definizione del campo di lavoro

r = -10:.02:10;
theta = 0:.02:2*pi; %Semipiano positivo delle y

[R, Theta] = meshgrid(r, theta);
x = R.*cos(Theta);
y = R.*sin(Theta);

%Integrazione delle derivate del potenziale phi
u = (Q/(2*pi))*(x./(x.^2 + y.^2)) - U;
v = (Q/(2*pi))*(y./(x.^2 + y.^2));
fx = (5*pi*log(x.^2 + y.^2))/4 - x;
fy = (1791925356007081*log(x.^2 + y.^2))/4503599627370496;

%Definizione della funzione di corrente psi

psi = ((Q*Theta)/(2*pi)) - U.*y;

%Plot delle isolinee psi
figure
[~,h] = contour(x,y,psi,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,psi,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \psi','FontSize',18)
plot(0,0,'ko','LineWidth',2)

%Definizione della funzione potenziale psi
phi = fx + fy;

%Plot delle isolinee phi

figure
[~,h] = contour(x,y,phi,50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
hold on
contour(x,y,phi,[0 0],'LineWidth',2);
colormap jet
colorbar
axis equal;
xlabel('x','FontSize',15)
ylabel('y','Rotation',0,'FontSize',18)
title('Isolinee \phi','FontSize',18)
plot(0,0,'ko','LineWidth',2)
