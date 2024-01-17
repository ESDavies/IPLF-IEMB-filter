function [y,Jac] = sensorMeasurement(X_truth,sensor_num,sensor_pos,sensor_index,Cells)

    a=2;
    epsilon= 25;
    phi = 500;
   
        
    heat= zeros(Cells);
    H = zeros(sensor_num,1);
    Jac= zeros(sensor_num,4);

%loop through the sensor positions

   for i = 1 : sensor_num

            coord = [(sensor_pos(1,i) - X_truth(1)).^(2);(sensor_pos(2,i) - X_truth(3)).^(2)];  %find x distance and y distance of object from each sensor at time k
            d = sqrt(coord(1)+coord(2)); %distance of object from each sensor at time k
            H(i) = phi/(d^a + epsilon); %Calculate sensor measurement using received signal strength indicator
            
            f= (d^a + epsilon)^(-2); %Calculate partial differential dh/df(d) 
            f_1=d^(a-2); %Calculate partial differential df(d)/dd
            Jac(i,1)=a*phi*(sensor_pos(1,i) - X_truth(1))*f*f_1; %Calculate x jacobian 
            Jac(i,3)=a*phi*(sensor_pos(2,i) - X_truth(3))*f*f_1; %Calculate y jacobian


         heat(sensor_index(1,i),sensor_index(2,i))= H(i); %Measurements arranged in matrix

   end

    y = reshape(heat(:,:),(size(heat,2)*size(heat,1)),1); %Measurements arranged in vector
    
end