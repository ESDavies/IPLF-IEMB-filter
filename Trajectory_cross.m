function [X_truth,t_birth,t_death]=Trajectory_cross(Nsteps,F,numtruth,Q,Area)

%Code modified from the source files in https://arxiv.org/abs/1203.2995


birthtime = zeros(1,numtruth);
Pmid = 0.1*eye(4);

simlen = Nsteps; % must be odd
midpoint = (simlen+1)/2;
%redundent but kept to keep the trajectories consistant
randperm((Nsteps-1)/2,numtruth);


%If one wants to change birth/death times, do it here
numfb=repelem(midpoint-1,numtruth);

chol_Q=chol(Q)';


% Initialise at time midpoint and propagate forward and backwards (I sum
% 100, 100 to center it in my region of interest)

X_truth=zeros(4*numtruth,Nsteps);
for j= 1 : numtruth
%%disp(midpoint);
    x = chol(Pmid)'*randn(size(F,1),1)+repmat([Area(1)/2;0;Area(2)/2;0],1);
    xlog{midpoint} = x;
    X_truth(4*j-3:4*j,midpoint)=x(:);
    xf = x; xb = x;

    for t = 1:numfb(j)
        % Run forward and backward simulation process
        xf=F*xf + chol_Q*randn(size(F,1),size(x,2));
        xb = F\(xb + chol_Q*randn(size(F,1),size(x,2)));  
        xlog{midpoint-t} = xb(:,midpoint-t>birthtime(j));
        xlog{midpoint+t} = xf; % note that all targets exist after midpoint(j)   
        X_truth(4*j-3:4*j,midpoint-t)=xb(:);
        X_truth(4*j-3:4*j,midpoint+t)=xf(:);
    end

end
