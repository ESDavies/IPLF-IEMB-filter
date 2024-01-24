%This function runs scenario 2 for all the Gaussian multi-target filters in

%E. S. Davies and Á. F. García-Fernández, "Information Exchange Track-Before-Detect Multi-Bernoulli Filter for Superpositional Sensors,"
% in IEEE Transactions on Signal Processing, vol. 72, pp. 607-621, 2024, doi: 10.1109/TSP.2024.3349769

%The proposed filter is the Iterated Posterior Linearisation Filter (IPLF)
%implementation of the IEMB filter (IEMB-IPLF).

%This scenario considers a multi-Bernoulli birth model and each Bernoulli with small spatial
%uncertainty

% Author: Elinor S. Davies

clear
% close all hidden
rng(2,'twister')
%number of targets
numtruth = 5;

%create Sensor Network
Area=[120 120];
Cells=[12 12];
cell_width = Area(1)/Cells(1);

%initalise measurement vector
[sensor_pos , sensor_index]  = createSensorNetwork(Cells,Area);
sensor_num = length(sensor_pos);

%set target state dimension
n=4;

%Nearly constant velocity model
T=1;
F=[1 T 0 0;0 1 0 0; 0 0 1 T;0 0 0 1];
sigmaU=0.5;
Q=(sigmaU)^2*[T^3/3 T^2/2 0 0; T^2/2 T 0 0;0 0 T^3/3 T^2/2; 0 0 T^2/2 T];

%Set the measurement covariance
R=1*eye(Cells(1)*Cells(2));
p_s_constant=0.99;

%Number of steps
Nsteps=81;

nmc=100; % Number of Monte Carlo Loops

%We generate ground truth according to the model
[X_truth] = Trajectory_cross(Nsteps, F, numtruth, Q, Area);

% Makes sure the target trajectories are contained within the surveillence region
[X_truth,m_b,t_birth,t_death,targetLeg] = Trajectory_region(X_truth,numtruth,Nsteps);
figure(1)
hold on
plot(sensor_pos(1,:),sensor_pos(2,:),'.k','Markersize',4)
targetLeg=cat(1,targetLeg,"Target Birth","Sensor Position");
legend(targetLeg);
% legend(targetLeg,'Location','southeast');

% saveas(figure(1),'Trajectories\Trajectories_with_birth.fig')
% saveas(figure(1),'Trajectories\Trajectories_with_birth.png')
% saveas(figure(1),'Trajectories\Trajectories_with_birth.pdf','pdf')
% saveas(figure(1),'Trajectories\Trajectories_with_birth.esp','epsc')

%Birth model

m_b=[75,0,100,0;...
    20,0,105,0;...
    110,0,50,0;...
    110,0,85,0]';

P_b=diag([10 10 10 10]);


P_ini=repmat(P_b,1,1,numtruth);

%Set the number of iterations for IPLF
it=20;

KLD=true; %If true IPLF stops the iteration once KLD becomes less than kld
kld=10^(-1);

%Define GOSPA variables
iplf_squared_gospa_t_tot=zeros(1,Nsteps); %Total error
iplf_squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
iplf_squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
iplf_squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error

%Define GOSPA variables for UKF
ukf_squared_gospa_t_tot=zeros(1,Nsteps); %Total error
ukf_squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
ukf_squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
ukf_squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error

%Define GOSPA variables for IEKF
iekf_squared_gospa_t_tot=zeros(1,Nsteps); %Total error
iekf_squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
iekf_squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
iekf_squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error

%Define GOSPA variables for EKF
ekf_squared_gospa_t_tot=zeros(1,Nsteps); %Total error
ekf_squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
ekf_squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
ekf_squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error

%Define GOSPA variables IMB-IPLF
imb_iplf_squared_gospa_t_tot=zeros(1,Nsteps); %Total error
imb_iplf_squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
imb_iplf_squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
imb_iplf_squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error

%Define GOSPA variables for IMB-UKF
imb_ukf_squared_gospa_t_tot=zeros(1,Nsteps); %Total error
imb_ukf_squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
imb_ukf_squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
imb_ukf_squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error


for q=1 :nmc

    tic

    z=measure(X_truth,t_birth,t_death,R,Nsteps,n,sensor_num,sensor_pos,sensor_index,Cells);

    %IPLF-IEMB
    [iplf_squared_gospa,iplf_gospa_loc,iplf_gospa_mis,iplf_gospa_fal] = TBDkFil(X_truth,z,m_b,t_birth,t_death,P_ini,p_s_constant,it,KLD,kld,Cells,cell_width,n,F,Q,R,Nsteps,sensor_pos,sensor_index,sensor_num);
    %UKF-IEMB
    [ukf_squared_gospa,ukf_gospa_loc,ukf_gospa_mis,ukf_gospa_fal] = UKFil(X_truth,z,m_b,t_birth,t_death,P_ini,p_s_constant,Cells,cell_width,n,F,Q,R,Nsteps,sensor_pos,sensor_index,sensor_num);
    %IEKF-IEMB
    [iekf_squared_gospa,iekf_gospa_loc,iekf_gospa_mis,iekf_gospa_fal]= IEKFil(X_truth,z,m_b,t_birth,t_death,P_ini,p_s_constant,it,KLD,kld,Cells,cell_width,n,F,Q,R,Nsteps,sensor_pos,sensor_index,sensor_num);
    %EKF-IEMB
    [ekf_squared_gospa,ekf_gospa_loc,ekf_gospa_mis,ekf_gospa_fal]= EKFil(X_truth,z,m_b,t_birth,t_death,P_ini,p_s_constant,Cells,cell_width,n,F,Q,R,Nsteps,sensor_pos,sensor_index,sensor_num);
    %IPLF-IMB
    [imb_iplf_squared_gospa,imb_iplf_gospa_loc,imb_iplf_gospa_mis,imb_iplf_gospa_fal] = IMB_IPLF(X_truth,z,m_b,t_birth,t_death,P_ini,p_s_constant,it,KLD,kld,Cells,cell_width,n,F,Q,R,Nsteps,sensor_pos,sensor_index,sensor_num);
    %UKF-IMB
    [imb_ukf_squared_gospa,imb_ukf_gospa_loc,imb_ukf_gospa_mis,imb_ukf_gospa_fal] = IMB_IPLF(X_truth,z,m_b,t_birth,t_death,P_ini,p_s_constant,1,false,kld,Cells,cell_width,n,F,Q,R,Nsteps,sensor_pos,sensor_index,sensor_num);

    %IPLF
    iplf_squared_gospa_t_tot=iplf_squared_gospa_t_tot+iplf_squared_gospa;
    iplf_squared_gospa_loc_t_tot=iplf_squared_gospa_loc_t_tot+iplf_gospa_loc;
    iplf_squared_gospa_false_t_tot=iplf_squared_gospa_false_t_tot+iplf_gospa_fal;
    iplf_squared_gospa_mis_t_tot=iplf_squared_gospa_mis_t_tot+iplf_gospa_mis;
    %UKF
    ukf_squared_gospa_t_tot=ukf_squared_gospa_t_tot+ukf_squared_gospa;
    ukf_squared_gospa_loc_t_tot=ukf_squared_gospa_loc_t_tot+ukf_gospa_loc;
    ukf_squared_gospa_false_t_tot=ukf_squared_gospa_false_t_tot+ukf_gospa_fal;
    ukf_squared_gospa_mis_t_tot=ukf_squared_gospa_mis_t_tot+ukf_gospa_mis;
    %IEKF
    iekf_squared_gospa_t_tot=iekf_squared_gospa_t_tot+iekf_squared_gospa;
    iekf_squared_gospa_loc_t_tot=iekf_squared_gospa_loc_t_tot+iekf_gospa_loc;
    iekf_squared_gospa_false_t_tot=iekf_squared_gospa_false_t_tot+iekf_gospa_fal;
    iekf_squared_gospa_mis_t_tot=iekf_squared_gospa_mis_t_tot+iekf_gospa_mis;
    %EKF
    ekf_squared_gospa_t_tot=ekf_squared_gospa_t_tot+ekf_squared_gospa;
    ekf_squared_gospa_loc_t_tot=ekf_squared_gospa_loc_t_tot+ekf_gospa_loc;
    ekf_squared_gospa_false_t_tot=ekf_squared_gospa_false_t_tot+ekf_gospa_fal;
    ekf_squared_gospa_mis_t_tot=ekf_squared_gospa_mis_t_tot+ekf_gospa_mis;
    %IMB-IPLF
    imb_iplf_squared_gospa_t_tot=imb_iplf_squared_gospa_t_tot+imb_iplf_squared_gospa;
    imb_iplf_squared_gospa_loc_t_tot=imb_iplf_squared_gospa_loc_t_tot+imb_iplf_gospa_loc;
    imb_iplf_squared_gospa_false_t_tot=imb_iplf_squared_gospa_false_t_tot+imb_iplf_gospa_fal;
    imb_iplf_squared_gospa_mis_t_tot=imb_iplf_squared_gospa_mis_t_tot+imb_iplf_gospa_mis;
    %IMB-UKF
    imb_ukf_squared_gospa_t_tot=imb_ukf_squared_gospa_t_tot+imb_ukf_squared_gospa;
    imb_ukf_squared_gospa_loc_t_tot=imb_ukf_squared_gospa_loc_t_tot+imb_ukf_gospa_loc;
    imb_ukf_squared_gospa_false_t_tot=imb_ukf_squared_gospa_false_t_tot+imb_ukf_gospa_fal;
    imb_ukf_squared_gospa_mis_t_tot=imb_ukf_squared_gospa_mis_t_tot+imb_ukf_gospa_mis;



    t=toc;
    display(['Completed iteration number ', num2str(q),' time ', num2str(t), ' sec'])

end

%Root mean square GOSPA errors at each time step for IPLF
iplf_rms_gospa_t=sqrt(iplf_squared_gospa_t_tot/nmc);
iplf_rms_gospa_loc_t=sqrt(iplf_squared_gospa_loc_t_tot/nmc);
iplf_rms_gospa_false_t=sqrt(iplf_squared_gospa_false_t_tot/nmc);
iplf_rms_gospa_mis_t=sqrt(iplf_squared_gospa_mis_t_tot/nmc);

%Root mean square GOSPA errors at each time step for UKF
ukf_rms_gospa_t=sqrt(ukf_squared_gospa_t_tot/nmc);
ukf_rms_gospa_loc_t=sqrt(ukf_squared_gospa_loc_t_tot/nmc);
ukf_rms_gospa_false_t=sqrt(ukf_squared_gospa_false_t_tot/nmc);
ukf_rms_gospa_mis_t=sqrt(ukf_squared_gospa_mis_t_tot/nmc);

%Root mean square GOSPA errors at each time step for IEKF
iekf_rms_gospa_t=sqrt(iekf_squared_gospa_t_tot/nmc);
iekf_rms_gospa_loc_t=sqrt(iekf_squared_gospa_loc_t_tot/nmc);
iekf_rms_gospa_false_t=sqrt(iekf_squared_gospa_false_t_tot/nmc);
iekf_rms_gospa_mis_t=sqrt(iekf_squared_gospa_mis_t_tot/nmc);

%Root mean square GOSPA errors at each time step for EKF
ekf_rms_gospa_t=sqrt(ekf_squared_gospa_t_tot/nmc);
ekf_rms_gospa_loc_t=sqrt(ekf_squared_gospa_loc_t_tot/nmc);
ekf_rms_gospa_false_t=sqrt(ekf_squared_gospa_false_t_tot/nmc);
ekf_rms_gospa_mis_t=sqrt(ekf_squared_gospa_mis_t_tot/nmc);

%Root mean square GOSPA errors at each time step for IMB-IPLF
imb_iplf_rms_gospa_t=sqrt(imb_iplf_squared_gospa_t_tot/nmc);
imb_iplf_rms_gospa_loc_t=sqrt(imb_iplf_squared_gospa_loc_t_tot/nmc);
imb_iplf_rms_gospa_false_t=sqrt(imb_iplf_squared_gospa_false_t_tot/nmc);
imb_iplf_rms_gospa_mis_t=sqrt(imb_iplf_squared_gospa_mis_t_tot/nmc);

%Root mean square GOSPA errors at each time step for IMB-UKF
imb_ukf_rms_gospa_t=sqrt(imb_ukf_squared_gospa_t_tot/nmc);
imb_ukf_rms_gospa_loc_t=sqrt(imb_ukf_squared_gospa_loc_t_tot/nmc);
imb_ukf_rms_gospa_false_t=sqrt(imb_ukf_squared_gospa_false_t_tot/nmc);
imb_ukf_rms_gospa_mis_t=sqrt(imb_ukf_squared_gospa_mis_t_tot/nmc);



%Root mean square GOSPA errors across all time steps for IPLF
iplf_rms_gospa_tot=sqrt(sum(iplf_squared_gospa_t_tot)/(nmc*Nsteps));
iplf_rms_gospa_loc_tot=sqrt(sum(iplf_squared_gospa_loc_t_tot)/(nmc*Nsteps));
iplf_rms_gospa_false_tot=sqrt(sum(iplf_squared_gospa_false_t_tot)/(nmc*Nsteps));
iplf_rms_gospa_mis_tot=sqrt(sum(iplf_squared_gospa_mis_t_tot)/(nmc*Nsteps));

%Root mean square GOSPA errors across all time steps for UKF
ukf_rms_gospa_tot=sqrt(sum(ukf_squared_gospa_t_tot)/(nmc*Nsteps));
ukf_rms_gospa_loc_tot=sqrt(sum(ukf_squared_gospa_loc_t_tot)/(nmc*Nsteps));
ukf_rms_gospa_false_tot=sqrt(sum(ukf_squared_gospa_false_t_tot)/(nmc*Nsteps));
ukf_rms_gospa_mis_tot=sqrt(sum(ukf_squared_gospa_mis_t_tot)/(nmc*Nsteps));

%Root mean square GOSPA errors across all time steps for IEKF
iekf_rms_gospa_tot=sqrt(sum(iekf_squared_gospa_t_tot)/(nmc*Nsteps));
iekf_rms_gospa_loc_tot=sqrt(sum(iekf_squared_gospa_loc_t_tot)/(nmc*Nsteps));
iekf_rms_gospa_false_tot=sqrt(sum(iekf_squared_gospa_false_t_tot)/(nmc*Nsteps));
iekf_rms_gospa_mis_tot=sqrt(sum(iekf_squared_gospa_mis_t_tot)/(nmc*Nsteps));

%Root mean square GOSPA errors across all time steps for EKF
ekf_rms_gospa_tot=sqrt(sum(ekf_squared_gospa_t_tot)/(nmc*Nsteps));
ekf_rms_gospa_loc_tot=sqrt(sum(ekf_squared_gospa_loc_t_tot)/(nmc*Nsteps));
ekf_rms_gospa_false_tot=sqrt(sum(ekf_squared_gospa_false_t_tot)/(nmc*Nsteps));
ekf_rms_gospa_mis_tot=sqrt(sum(ekf_squared_gospa_mis_t_tot)/(nmc*Nsteps));

%Root mean square GOSPA errors across all time steps for IMB-IPLF
imb_iplf_rms_gospa_tot=sqrt(sum(imb_iplf_squared_gospa_t_tot)/(nmc*Nsteps));
imb_iplf_rms_gospa_loc_tot=sqrt(sum(imb_iplf_squared_gospa_loc_t_tot)/(nmc*Nsteps));
imb_iplf_rms_gospa_false_tot=sqrt(sum(imb_iplf_squared_gospa_false_t_tot)/(nmc*Nsteps));
imb_iplf_rms_gospa_mis_tot=sqrt(sum(imb_iplf_squared_gospa_mis_t_tot)/(nmc*Nsteps));

%Root mean square GOSPA errors across all time steps for IMB_UKF
imb_ukf_rms_gospa_tot=sqrt(sum(imb_ukf_squared_gospa_t_tot)/(nmc*Nsteps));
imb_ukf_rms_gospa_loc_tot=sqrt(sum(imb_ukf_squared_gospa_loc_t_tot)/(nmc*Nsteps));
imb_ukf_rms_gospa_false_tot=sqrt(sum(imb_ukf_squared_gospa_false_t_tot)/(nmc*Nsteps));
imb_ukf_rms_gospa_mis_tot=sqrt(sum(imb_ukf_squared_gospa_mis_t_tot)/(nmc*Nsteps));

leg={'IEMB-IPLF','IEMB-UKF','IEMB-EKF','IEMB-IEKF','IMB-IPLF','IMB-UKF'};


figure(1)
clf
plot(1:Nsteps,iplf_rms_gospa_t,'blue','Linewidth',1.3)
hold on
plot(1:Nsteps,ukf_rms_gospa_t,'--','Linewidth',1.3)
plot(1:Nsteps,ekf_rms_gospa_t,'--','Linewidth',1.3)
plot(1:Nsteps,iekf_rms_gospa_t,'--','Linewidth',1.3)
plot(1:Nsteps,imb_iplf_rms_gospa_t,'--','Linewidth',1.3)
plot(1:Nsteps,imb_ukf_rms_gospa_t,'--','Linewidth',1.3)
lim =axis;
ax=[0 Nsteps lim(3:4)];
axis(ax);
grid on
xlabel('Time step')
ylabel('RMS GOSPA error')
legend(leg,'Location','northwest')
set(gca,'FontSize',15)


figure(2)
clf
plot(1:Nsteps,iplf_rms_gospa_loc_t,'blue','Linewidth',1.3)
hold on
plot(1:Nsteps,ukf_rms_gospa_loc_t,'--','Linewidth',1.3)
plot(1:Nsteps,ekf_rms_gospa_loc_t,'--','Linewidth',1.3)
plot(1:Nsteps,iekf_rms_gospa_loc_t,'--','Linewidth',1.3)
plot(1:Nsteps,imb_iplf_rms_gospa_loc_t,'--','Linewidth',1.3)
plot(1:Nsteps,imb_ukf_rms_gospa_loc_t,'--','Linewidth',1.3)
axis(ax);
grid on
xlabel('Time step')
ylabel('RMS GOSPA localisation error')
legend(leg,'Location','northwest')
set(gca,'FontSize',15)

figure(3)
clf
plot(1:Nsteps,iplf_rms_gospa_false_t,'blue','Linewidth',1.3)
hold on
plot(1:Nsteps,ukf_rms_gospa_false_t,'--','Linewidth',1.3)
plot(1:Nsteps,ekf_rms_gospa_false_t,'--','Linewidth',1.3)
plot(1:Nsteps,iekf_rms_gospa_false_t,'--','Linewidth',1.3)
plot(1:Nsteps,imb_iplf_rms_gospa_false_t,'--','Linewidth',1.3)
plot(1:Nsteps,imb_ukf_rms_gospa_false_t,'--','Linewidth',1.3)
axis(ax);
grid on
xlabel('Time step')
ylabel('RMS GOSPA false target error')
legend(leg,'Location','northwest')
set(gca,'FontSize',15)

figure(4)
clf
plot(1:Nsteps,iplf_rms_gospa_mis_t,'blue','Linewidth',1.3)
hold on
plot(1:Nsteps,ukf_rms_gospa_mis_t,'--','Linewidth',1.3)
plot(1:Nsteps,ekf_rms_gospa_mis_t,'--','Linewidth',1.3)
plot(1:Nsteps,iekf_rms_gospa_mis_t,'--','Linewidth',1.3)
plot(1:Nsteps,imb_iplf_rms_gospa_mis_t,'--','Linewidth',1.3)
plot(1:Nsteps,imb_ukf_rms_gospa_mis_t,'--','Linewidth',1.3)
axis(ax);
grid on
xlabel('Time step')
ylabel('RMS GOSPA missed target error')
legend(leg,'Location','northwest')
set(gca,'FontSize',15)

save('Workspace')

writematrix([...
    "Filter",    "Error tot",            "Error loc.",               "Error Fal.",                 "Error Mis.";...
    "IEMB-IPLF", iplf_rms_gospa_tot,     iplf_rms_gospa_loc_tot,     iplf_rms_gospa_false_tot,     iplf_rms_gospa_mis_tot;...
    "IEMB-UKF",  ukf_rms_gospa_tot,      ukf_rms_gospa_loc_tot,      ukf_rms_gospa_false_tot,      ukf_rms_gospa_mis_tot;...
    "IEMB-IEKF", iekf_rms_gospa_tot,     iekf_rms_gospa_loc_tot,     iekf_rms_gospa_false_tot,     iekf_rms_gospa_mis_tot;...
    "IEMB-EKF",  ekf_rms_gospa_tot,      ekf_rms_gospa_loc_tot,      ekf_rms_gospa_false_tot,      ekf_rms_gospa_mis_tot;...
    "IMB-IPLF",  imb_iplf_rms_gospa_tot, imb_iplf_rms_gospa_loc_tot, imb_iplf_rms_gospa_false_tot, imb_iplf_rms_gospa_mis_tot;...
    "IMB-UKF",   imb_ukf_rms_gospa_tot,  imb_ukf_rms_gospa_loc_tot,  imb_ukf_rms_gospa_false_tot,  imb_ukf_rms_gospa_mis_tot;...
    ],'Results Table.csv',"WriteMode","append");

% b=["total" "localisation" "false" "missed"];
% for i=1:4
%     saveas(figure(i),'Graphs\'+b(i)+'.fig')
%     saveas(figure(i),'Graphs\'+b(i)+'.png')
%     saveas(figure(i),'Graphs\'+b(i)+'.pdf','pdf')
%     saveas(figure(i),'Graphs\'+b(i)+'.esp','epsc')
% end
