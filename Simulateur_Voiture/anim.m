function [szOut]=anim

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% anim.m: Fonction qui fait l'animation
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% on récupère les data
load voit;
voit= ans';
save voit.data voit -ascii;
load ref;
ref= ans';
save ref.data ref -ascii;

% on dessine le premier
xytp=[voit(1,2),voit(1,3),voit(1,4),voit(1,5)];
xyt_r= [ref(1,2),ref(1,3),ref(1,4)];
draw(xytp,xyt_r,1);
refresh;

% on fait l'animation
cl0=clock;
  for index = 1:1:size(voit,1)
    while  etime(clock,cl0) < voit(index,1),
     etime(clock,cl0);
    end;
    xytp=[voit(index,2),voit(index,3),voit(index,4),voit(index,5)];
    xyt_r= [ref(index,2),ref(index,3),ref(index,4)];
    draw(xytp,xyt_r,index);
  end;
 disp('Taper sur une touche pour fermer la fenetre');
 pause;
 close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [szOut] = draw(Xv,Xr,index)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% draw.m: Fonction qui fait les tracés
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global hc;

% Définition de la géométrie du robot et du repère de référence
L=2;
p= Xv(1:2)';
theta= Xv(3);
phi= Xv(4);
Rv= [cos(theta) -sin(theta); sin(theta) cos(theta)];
Rphi= [cos(theta+phi) -sin(theta+phi); sin(theta+phi) cos(theta+phi)];
XYT	= [ 0 L nan 0  0 nan   -0.5*L  0.5*L nan -0.5*L -0.5*L nan...
          -0.5*L  0.5*L nan 0.5*L L nan L 0.5*L nan -0.3*L 0.3*L nan...
          -0.3*L 0.3*L; 
          0 0 nan -0.7*L 0.7*L nan 0.5*L 0.5*L nan -0.5*L 0.5*L nan...
          -0.5*L -0.5*L nan -0.5*L 0 nan 0 0.5*L nan -0.7*L -0.7*L nan...
          0.7*L 0.7*L];
RAv= [-0.3 0.3; 0 0];

pr= Xr(1:2)';
thetar= Xr(3);
Rr= [cos(thetar) -sin(thetar); sin(thetar) cos(thetar)];
XYTr= 2*[0 1 nan 0 0 nan 1 1-0.2 nan 1 1-0.2 nan 0 0.2 nan 0 -0.2;...
       0 0 nan 0 1 nan 0 0.2 nan 0 -0.2 nan 1 1-0.2 nan 1 1-0.2];

for i=1:26
  XYT(:,i) = Rv*XYT(:,i) + p;
end;

for i=1:2
  RAv(:,i) = Rphi*RAv(:,i) + p+Rv*[L;0];
end;

for i=1:17
  XYTr(:,i) = Rr*XYTr(:,i) + pr;
end;

% On dessine les objets
if index==1   
% Initialisation des propriétés de la figure
  figure(1);
  clf;
  set(gcf, 'Name', 'Animation Window');
  set(gcf,'Position',[ 10 50 950 600]);
  set(1,'BackingStore','on');
  hold on;
  set(gca, 'UserData', hc,'NextPlot', 'add','Visible', 'on','DataAspectRatio', [1 1 1], ...
	  'Color', 'k', 'SortMethod', 'childorder','Xlim',[-2 40],'Ylim',[-10 15]);

  hc(1)=plot(XYT(1,:),XYT(2,:),'g-','Linewidth',2);
  hc(2)=plot(RAv(1,:),RAv(2,:),'r-','Linewidth',2);
  hc(3)= plot(XYTr(1,:),XYTr(2,:),'w','Linewidth',2);
  drawnow;
 
else
  set(hc(1),'XData',XYT(1,:),'YData',XYT(2,:));
  set(hc(2),'XData',RAv(1,:),'YData',RAv(2,:));
  set(hc(3),'XData',XYTr(1,:),'YData',XYTr(2,:));
  drawnow;

end
