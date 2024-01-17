function [S,C,mean_y,mean_yy,sigma_points,W,y]= sigma_approximation(mean_pred,P,n,sensor_pos,sensor_index,sensor_num,Cells)

    % define variables
    % number of sigma point
    sigmanum = 2*n+1;
    sigma_points= zeros(4,sigmanum);
    
    %Find the weights for the sigma points
    W = zeros(1,sigmanum);
    W(1,1) = 1/3;
    W(1,2:2*n+1) = (1-W(1,1))/(2*n);

        %Generate Sigma points

    [chol_var_mult]=chol((n/(1-W(1,1))*P));
    sigma_points=[zeros(n,1),chol_var_mult',-chol_var_mult'];
    sigma_points=repmat(mean_pred,1,length(W))+sigma_points;
    
    % set covariance matrices to zero
    S = zeros(sensor_num,sensor_num);
    C = zeros(n,sensor_num);
    mean_y = zeros(sensor_num,1);
    mean_yy = zeros(sensor_num,sensor_num);

    %Generate measurements for the sigma points
    for j=1 : sigmanum
        y(:,j) = sensorMeasurement(sigma_points(:,j),sensor_num, sensor_pos, sensor_index,Cells);
    end

    %calculate the mean measurement of the sigma points           
    for j = 1:size(sigma_points,2)
        Wy=W(j) * y(:,j);
        mean_y = mean_y + Wy;
        mean_yy=mean_yy + Wy*y(:,j)';
    end

    % Find the covariance for each sigma point
    for j=1 : sigmanum
        S = S + W(1,j)*(y(:,j)- mean_y)*(y(:,j)- mean_y)';
        C = C + W(1,j)*(sigma_points(:,j)- mean_pred)*(y(:,j)- mean_y)';
    end

end