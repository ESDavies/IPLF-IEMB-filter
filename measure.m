function [z]=measure(X_truth,t_birth,t_death,R,Nsteps,n,sensor_num,sensor_pos,sensor_index,Cells)
    z=zeros(sensor_num,Nsteps);

    for k=1 : Nsteps
        index_targets=find(and(t_birth<=k,t_death>=k)); %get the index positions for detected targets
        N_detected_targets=length(index_targets); %find number of targets
        if N_detected_targets > 0
            for i= 1 : N_detected_targets %Find measurements for each target and sum them together
                z(:,k) = z(:,k)+ sensorMeasurement(X_truth((n*(index_targets(i)-1)+1):index_targets(i)*n,k),sensor_num,sensor_pos,sensor_index,Cells);
            end
        end
        z(:,k)=z(:,k)+chol(R)*randn(size(z(:,k))); %Add noise to the measurements
    end
end