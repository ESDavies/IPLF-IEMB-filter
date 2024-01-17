function [sensor_pos , sensor_index] = createSensorNetwork(Cells,Area)
    
    cell_width = Area(1)/Cells(1);
    cell_height = Area(2)/Cells(2);

    x_index = (1: Cells(1));
    y_index = (1: Cells(2));
    sensor_pos = zeros(2,Cells(1)*Cells(2));
    [X_index, Y_index] = meshgrid(x_index , y_index);

    sensor_index(1,:) = X_index(:);
    sensor_index(2,:) = Y_index(:);
    sensor_pos(1,:) = sensor_index(1,:)*cell_width - cell_width/2;
    sensor_pos(2,:) = sensor_index(2,:)*cell_height - cell_height/2;

end
