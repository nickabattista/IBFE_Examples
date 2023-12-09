function test_Spinning()

% Loads: "XY" 
load('XY_Data.mat');
s0 = XY(:,2);
s1 = XY(:,3);
    
XY_Aux = zeros( length(XY(:,1)), 2);

tVec=0:0.01:2;
for n=1:length(tVec)
    
    time=tVec(n);
    

    
    if time<1
       xAux=s0*cos(2*3.14*0.25*time) - s1*sin(2*3.14*0.25*time); 
       yAux=s0*sin(2*3.14*0.25*time) + s1*cos(2*3.14*0.25*time)+.5*time;
    end
    
    
    plot(XY(:,2),XY(:,3),'b.','MarkerSize',10); hold on;
    plot(xAux,yAux,'r.','MarkerSize',10); hold on;
    axis([-0.25 0.5 -0.25 0.5])
    pause(0.1);
    clf;
    
    
end