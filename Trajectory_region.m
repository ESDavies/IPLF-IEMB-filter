function [X_truth,m_b,t_birth,t_death,targetLeg] = Trajectory_region(X_truth,numtruth,Nsteps)

t_birth= [1:numtruth];
t_death= [Nsteps-1:-1:Nsteps-numtruth];

m_bK=zeros(4,numtruth);
midpoint = (Nsteps+1)/2;


x_pos=4*([1:numtruth]-1)+1;
y_pos=4*([1:numtruth]-1)+3;

targetLeg=[];

for j=1:numtruth
%targets must be born in first half of simulation    
    for k=1: midpoint-1
        %If the target is in the surveillence region by 10m then trajectory
        %is unchanged
        if X_truth(x_pos(j),k)>10 && X_truth(x_pos(j),k)<110 && X_truth(y_pos(j),k)>10 && X_truth(y_pos(j),k)<110
        else
            %If the target is outside of the surveillence region or less
            %that 10m from the edge of the surveillence region then we set
            %set the target birth to the next time step
            if ismember(k+1,t_birth)
                %If a target is already being born in this time step then
                %we skip to the next time-step because two targets can't be
                %born in the same time step
                continue
            end


            t_birth(j)=k+1;
        end
        %
        X_truth(x_pos(j):y_pos(j)+1,1:t_birth(j)-1)= zeros(4,length(1:t_birth(j)-1));
        m_bK(:,j)=X_truth(x_pos(j):y_pos(j)+1,t_birth(j));
        
    end
    for k= Nsteps:-1:midpoint+1
        if X_truth(x_pos(j),k)>10 && X_truth(x_pos(j),k)<110 && X_truth(y_pos(j),k)>10 && X_truth(y_pos(j),k)<110
        else
            if ismember(k-1,t_death)
                continue
            end
            t_death(j)=k-1;
        end
        X_truth(x_pos(j):y_pos(j)+1,t_death(j)+1:Nsteps)= zeros(4,length(t_death(j)+1:Nsteps));
    end
    
     figure(1)
     hold on
    axis([0 120 0 120])
    plot(X_truth(x_pos(j),t_birth(j):t_death(j)),X_truth(y_pos(j),t_birth(j):t_death(j)),'Linewidth',1.3)    
    target=['Target '+string(j)];
    targetLeg=cat(1,targetLeg,target);
    xlabel('x axis (m)')
    ylabel('y axis (m)')


end
   plot(X_truth(sub2ind(size(X_truth), x_pos,t_birth)),X_truth(sub2ind(size(X_truth),y_pos,t_birth)),'*k','Markersize',9)

% for a known_birth
    m_b=round(m_bK.*[1 0 1 0]');
    m_b=m_bK;
        
end