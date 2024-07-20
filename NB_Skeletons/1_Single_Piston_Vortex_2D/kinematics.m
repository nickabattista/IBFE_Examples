function kinematics()

tVec=0:0.001:1;
maxSpeed=0.375;
d1=-0.4+4*maxSpeed*(0.5*tVec.^2-1/3*tVec.^3);
d2=-0.4+4*maxSpeed*tVec.*(1-tVec);
plot(tVec,d1,'b-','LineWidth',5); hold on;
plot(tVec,d2,'r-','LineWidth',5); hold on;
pause();

yL = linspace(-0.05,0.05,50);
xL = -0.4*ones(size(yL));

xT=linspace(-0.4,-0.4+0.25,10);
yT=0.05*ones(size(xT));
yB=-yT;

dt=5e-4;
tVec=0:dt:1;

X=xL;
maxSpeed=0.375;

for n=1:length(tVec)
   
    tTilde=mod(tVec(n),2);
    
    X = X + 4.0*(maxSpeed)*(tTilde)*(1.0-tTilde)*dt;
    
    if mod(n,10)==0
        plot(xT,yT,'k-','LineWidth',3); hold on;
        plot(xT,yB,'k-','LineWidth',3); hold on;
        plot(X,yL,'.','MarkerSize',20); hold on;
        axis([-1 1 -1 1]);
        pause(0.1)
    end
    
end