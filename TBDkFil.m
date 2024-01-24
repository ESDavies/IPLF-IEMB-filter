function [squared_gospa_t_tot,squared_gospa_loc_t_tot,squared_gospa_mis_t_tot,squared_gospa_false_t_tot] = TBDkFil(X_truth,z,m_b,t_birth,t_death,P_ini,p_s_constant,it,KLD,kld,Cells,cell_width,n,F,Q,R,Nsteps,sensor_pos,sensor_index,sensor_num)

%Bernoulli birth probability
p_b= 0.0001;

%Birth covariance matrix
P_b =P_ini;

%Define GOSPA variables
squared_gospa_t_tot=zeros(1,Nsteps);
squared_gospa_loc_t_tot=zeros(1,Nsteps); %Localisation error
squared_gospa_false_t_tot=zeros(1,Nsteps); %False target error
squared_gospa_mis_t_tot=zeros(1,Nsteps); %Misdetection error
p=2;
alpha=2;

%Set intial values for variables
P_pred = []; %Predicted Covariance
mean_pred = []; %Predicted Mean
r=[]; % Probability of existance
L=0; % Number of potential targets
for k=1 : Nsteps

    detected_targets=and(t_birth<=k,t_death>=k); %get the index positions for detected targets
    N_detected_targets = sum(detected_targets);

    %Prediction        
    for i=1 : size(mean_pred,2) % Predict the state of targets that survived previous time step
        mean_pred(:,i) = F*mean_pred(:,i);

        % Remove targets outside of surveillence region by setting
        % probability of survial to zero
        if mean_pred(1,i)<0 ||  mean_pred(1,i)>120 ||  mean_pred(3,i)<0 ||  mean_pred(3,i)>120
            p_s=0;
        else
            p_s=p_s_constant;
        end
        P_pred(:,:,i) = F*P_pred(:,:,i)*F' + Q;            
        r(i) = p_s*r(i);            
    end
    % Birth model
    mean_pred = [mean_pred m_b];
    P_pred = cat(3,P_pred,P_b);
    r = [r ones(1,size(m_b,2))*p_b];
    L = size(mean_pred,2); 

    % Initalise covariance matrices 
    S = zeros(sensor_num,sensor_num,L);
    S_corr = zeros(sensor_num,sensor_num,L); 
    C = zeros(n,sensor_num,L);
    A = zeros(sensor_num,n,L);
    b = zeros(sensor_num,L);
    omega = zeros(sensor_num,sensor_num,L);

    % Initalise matrices for the predicted measurment mean
    mean_y = zeros(sensor_num,L);
    mean_yy = zeros(sensor_num,sensor_num,L);
    mean_z = zeros(sensor_num,L);            
    z_corr=zeros(sensor_num,L);

    for i=1:L % Iterate over each potential target
        %Use sigma points to find predicted measurments mean and
        %the mean product of the predicted measurments and its transpose
        [~,~,mean_y(:,i),mean_yy(:,:,i)]= sigma_approximation(mean_pred(:,i),P_pred(:,:,i),n,sensor_pos,sensor_index,sensor_num,Cells);

    end
    for i=1:L % Iterate over each potential target
        for u=1:L % Iterate over each potential target again                    
            if u~=i % find the correction for the measurements and the measurement covariance
                z_corr(:,i)=z_corr(:,i)+r(u)*mean_y(:,u);
                S_corr(:,:,i)=S_corr(:,:,i)+r(u)*mean_yy(:,:,u)-r(u)*r(u)*mean_y(:,u)*mean_y(:,u)';
            end
        end                
    end
    %IPLF
   u_p=mean_pred;
   W=P_pred;      


    for i=1 : L % Iterate over each potential target
        for j=1:it

            [S(:,:,i),C(:,:,i),mean_y(:,i),~]= sigma_approximation(u_p(:,i),W(:,:,i),n,sensor_pos,sensor_index,sensor_num,Cells);

            %Find A, b, omega
            A=C(:,:,i)'/W(:,:,i);
            b=mean_y(:,i)-A*u_p(:,i);
            omega=S(:,:,i)-A*W(:,:,i)*A';      

            %Calculate the predicted measurement mean and covariance                                                
            S(:,:,i) = A*P_pred(:,:,i)*A'+ omega+ R + S_corr(:,:,i);
            mean_z(:,i)=A*mean_pred(:,i)+ b + z_corr(:,i);
            
            %Set mu0 and var0 to equal the mean and covariance before update 
            mu0= u_p(:,i);
            var0= W(:,:,i);

            %update the mean and covariance
            u_p(:,i)=mean_pred(:,i)+P_pred(:,:,i)*A'/(S(:,:,i))*(z(:,k)-mean_z(:,i));
            W(:,:,i)=P_pred(:,:,i)-P_pred(:,:,i)*A'/(S(:,:,i))*A*P_pred(:,:,i);

            %Set mu1 and var1 to equal the mean and covariance after update
            mu1= u_p(:,i);
            var1= W(:,:,i);
            %Find  the KLD
            d=dist_kullback(mu0,var0,mu1,var1);
            if KLD %If KLD is set to true then we stop the iteration once d<kld
                if d<kld %Currently kld=10^(-5)
                    break
                end
            end
        end

        % Calulate the log likelihoods
        Log_L1 = -0.5*log(det(S(:,:,i)))-0.5*((z(:,k)- mean_z(:,i))'/S(:,:,i))*(z(:,k)- mean_z(:,i)); % If the target exists
        Log_L0 = -0.5*log(det(R+S_corr(:,:,i)))-0.5*((z(:,k)-z_corr(:,i))'/(R+S_corr(:,:,i)))*(z(:,k)-z_corr(:,i)); % If the target does not exists

        %Calculate the probabilities
        % Reduce rounding error by dividing both 
        % probabilities by the largest
        max_L = max(Log_L0,Log_L1);
        L1_prob= exp(Log_L1 - max_L);
        L0_prob= exp(Log_L0 - max_L);

        %Update the existance probability
        numerater=(r(i)*L1_prob);
        denominater=((1-r(i))*L0_prob + r(i)*L1_prob);            

        r(i) = numerater/denominater;

    end

    % Remove targets that don't exist
    mean_pred=u_p(:,r>=0.01);
    P_pred=W(:,:,r>=0.01);
    r=r(r>=0.01);

    %Start GOSPA

    X_truth_pos = reshape([detected_targets;and(false,detected_targets);detected_targets;and(false,detected_targets)],[],1);
    X_truth_pos = reshape(X_truth(X_truth_pos,k),2,[]);
    mean_pred_pos = reshape(mean_pred([r>=0.5;and(false,r>=0.5);r>=0.5;and(false,r>=0.5)]),2,[]);

    [d_gospa, ~, decomp_cost] = GOSPA(X_truth_pos, mean_pred_pos, p, cell_width/2, alpha);

    squared_gospa=d_gospa^2;
    gospa_loc=decomp_cost.localisation;
    gospa_mis=decomp_cost.missed;
    gospa_fal=decomp_cost.false;

    squared_gospa_t_tot(k)=squared_gospa;
    squared_gospa_loc_t_tot(k)=gospa_loc;
    squared_gospa_false_t_tot(k)=gospa_fal;
    squared_gospa_mis_t_tot(k)=gospa_mis;
    %End GOSPA 

end    
end
