%%Plot del campo di velocità risultante da una sorgente puntiforme
%Definisco il piano di lavoro in coordinate polari
r = -6:.02:6;
theta = 0:.02:2*pi;
[R,Theta] = meshgrid(r,theta)
x = R.*cos(Theta);
y = R.*sin(Theta);
x0 = 6;
y0 = 0;
%Definisco l'intensità della sorgente Q
Q = 5;

%Calcolo la funzione di corrente
psi = (Q.*Theta)/(2*pi);

%Plotto le isolinee psi
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

%Calcolo la funzione potenziale phi
phi = (Q/(2*pi)).*log((x.^2 + y.^2)./(x0.^2 + y0.^2));

%Plotto le isolinee phi
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

hold on
