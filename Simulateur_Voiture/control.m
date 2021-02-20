
function [output]=control(inp)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction qui calcule la loi de commande
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% On renomme les variables
  pr= inp(1:2); % position de la référence 
  thetar= inp(3); % orientation de la référence 
  Xr=[pr;thetar]; % Etat complet de la référence
  p= inp(4:5); % position de la voiture
  theta = inp(6) % orientation de la voiture
  phi = inp(7) % angle de braquage de la voiture
  X= [p;theta;phi];
  vr= inp(8); % vitesse linéaire de la référence
  wr= inp(9); % vitesse angulaire de la référence 
  
  %on fixe phi reference nul comme pour la commande LT 
  phir = 0;
  % Distance entre les essieus avant et arrière
  L=2;
  
  %gain commande lin exact de retour de sortie
  k1 = 9;
  k21 = 800;
  k22 = 850;
  k23 = 180;
  
  %calcul de p tilde pour x tilde
  p_tild = rot(-thetar)*(p-pr);
  
  %calcul xsi pour changement de variable
  xi1 = p_tild(1);
  xi2 = p_tild(2);
  xi3 = tan(theta-thetar);
  xi4 = tan(phi)/(L*cos(theta-thetar)^3);
  
  
  a1 = vr* tan(phir)/L;
  %^lin exact pour w1 sur p1 découplé du reste
  w1 = -k1*xi1 - a1*xi2 + vr;

  alpha = -2*a1*w1*xi3*(-3*a1*xi3^2 + 3*w1*xi4 - 2.5*a1)...
          +3*(-a1 + w1*xi4/(1+xi3^2))*xi3*xi4*w1^2;
      
  xi34 = ((1 + xi3^2)^1.5)/L + L*xi4^2/((1 + xi3^2)^1.5);
  beta = xi34*w1^2;
  
  dxi2 = xi3*w1 - a1 * xi1;
  ddxi2 = -2*a1*w1* xi3^2 + w1*(w1*xi4 - a1);
  
  %lin exacte par sortie pour p2
  w2 = (-alpha - k21*xi2 -k22*dxi2 -k23*ddxi2)/beta;

  w = [w1, w2];
  
  %changement de variable inverse pour avoir les u tilde
  M = [1/cos(theta-thetar) 0;
       0                   1];
  u = M*w';

% %Calcul de la commande du linearise tangent autour de l'equilibre (x,u) =
% %(0, 0) et phir nul ne marche que pour la 1ere partie de la trajectoire
% %(explications en 6)
%   u1 = -20* p_tild(1);
%   u2 = -0.1* p_tild(2)- 10 * (theta-thetar) - 2 * phi;
%   u= [u1;u2];
output= [u;X;Xr];

function out= rot(theta)
out= [cos(theta) -sin(theta); sin(theta) cos(theta)];
  