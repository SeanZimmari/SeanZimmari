%%Traformazione del flusso attorno ad un cerchio in flusso attorno a
%%profilo alare tramite Traformazione di Kutta-Joukowsky

%Dati in input
clear all
close all
clc
disp('-------------------------------------------')
disp('Kutta-Joukovsky transformation')
disp('-------------------------------------------')

v_inf = input(' Speed downstream of the cylinder [m/s]: ');

disp('-------------------------------------------')

x_center = input('Insert the x-coordinate of the center of the circle [m]: ');
y_center = input('Insert the y-coordinate of the center of the circle [m]: ');

z0 = x_center + 1i*y_center;

disp('-------------------------------------------')

alpha = input('Insert the angle of incidence [deg]: ');
alpha = alpha*pi/180;

disp('-------------------------------------------')

a = input("Insert cylinder's radius [m]: ");

disp('-------------------------------------------')

beta = asin(y_center/a);
b = a*cos(beta) - x_center;
e = x_center/b;


%Vt = V_inf*sin(alpha+beta) = gamma/(2*pi*r)------>Gamma =
%2*pi*r*V_inf*sin(alpha+beta)
theta = linspace(0,2*pi,200);
Gamma2pi = 2*v_inf*a*sin(alpha+beta);
Gamma = Gamma2pi/(2*pi);
rho = 1.225;


%Dimensioni del piano complesso
d = 8;
ll = -d-d*1i;
ur = d+d*1i;

%Definisco il piano complesso x-iy
x = linspace(real(ll),real(ur),100);
y = linspace(imag(ll),imag(ur),100);
[X,Y] = meshgrid(x,y);
Z = complex(X,Y);

%Escludo i punti all'interno del cilindro
for i=1:length(X)
    for j=1:length(Y)
        if(abs(Z(i,j) - z0) < a)
            Z(i,j) = NaN;
        end
    end
end

%Definisco il potenziale complesso per un caso generico
wf = @(z) (abs(z) > 1).*((-v_inf.*((z.*exp(1i*alpha)) + ((a^2)./z).*exp(-1i*alpha))) - 1i*(Gamma2pi/(2*pi)).*log(z./a));
w = wf(Z - z0);

%Plotto le isolinee w nel piano x-iy
L = rho*v_inf*Gamma;
z_circle = z0 + a*exp(1i*theta);
figure(1)
plot(0,0,'Marker','o','Color','k')
hold on
fill(real(z_circle),imag(z_circle),'y')
hold on
colormap jet
colorbar
[~,h] = contour(real(Z),imag(Z),imag(w),50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
axis([real(ll),real(ur),imag(ll),imag(ur)],'on','equal')
xlabel('x','FontSize',18)
ylabel('iy','Rotation',0,'FontSize',18)
title(strcat('Lift around cylinder = ',num2str(L),'N/m'))

%Trasformo il piano complesso x-iy nel piano epsilon-i(eta) tramite la
%trasformazione di Kutta-Joukowsky

ZETA = Z + (b^2./Z);
z_airfoil = z_circle + b^2./z_circle;
figure(2)
plot(0,0,'Marker','o','Color','k')
hold on
fill(real(z_airfoil),imag(z_airfoil),'k')
hold on
colormap jet
colorbar
[~,h] = contour(real(ZETA),imag(ZETA),imag(w),50,'LineWidth',1);
set(h,'ShowText','off','TextStep',get(h,'LevelStep'));
axis([real(ll),real(ur),imag(ll),imag(ur)],'on','equal')
xlabel('\xi','FontSize',18)
ylabel('i\eta','Rotation',0,'FontSize',18)
title(strcat('Lift around Airfoil = ',num2str(L),'N/m'))