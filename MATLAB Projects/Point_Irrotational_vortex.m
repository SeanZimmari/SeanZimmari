%%Plot del campo di velocità di un vortice irrotazionale puntiforme
%Definisco il piano su cui andrò a lavorare
r = -6:.02:6;
theta = 0:.02:2*pi;
[R, Theta]= meshgrid(r,theta);
x = R.*cos(Theta);
y = R.*sin(Theta);
x0 = 10;
y0 = 0;
r0 = sqrt(x0^2 + y0^2);
circ = 5;

%Calcolo la funzione psi
psi = -(circ/(2*pi)).*log((x.^2 + y.^2)./(x0.^2 + y0.^2));

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

%Calcolo la funzione phi
phi = (circ.*Theta)./(2*pi);

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