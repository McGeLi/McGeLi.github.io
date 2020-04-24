clc;clear all;
%% Parameter setting
Vel = 9;%Homogeneous
delta= 10;
depthEvent = 10;
error_tolerance = 10^-7;
hypoCenter = [10 3.0 10.0]';%True hypocenter
hypoCenter0 =[5 4 20]';%Initial guess of hypocenter
evolution_path(:,1) = hypoCenter0;
t0_actual = 0; % Real origin time
t0 = 2; % Initial origin time
sta = [1 0 0;0 1 0; 2.5 0 0; 0 6.5 0;
      10 7 0;0 2.5 0; 20 5 0; 13 2 0;
      0.5 10 0; 14 10 0]';%10 Station locations
staT = sta';
dist = sqrt(sum(abs(sta-(diag(hypoCenter)*...
     ones(size(sta)))).^2,1))';  % distance between hypocenter and stations
tObs = dist/Vel+t0_actual;%Observations

%% Initialize G Matrix, residuals and itertion step
G =  zeros(10,4);
r0 = zeros(1,10)+1;
iteration = 0;

%% Visualization switch
Visualize = 1;

if Visualize
        scatter(staT(:,1),staT(:,2),'x','LineWidth',2,'MarkerEdgeColor','black')
        hold on
        scatter(hypoCenter(1),hypoCenter(2),100,'*','MarkerEdgeColor','red')
        scatter(hypoCenter0(1),hypoCenter0(2),100,'d','MarkerEdgeColor','green','MarkerFaceColor','green')
        axis equal
        box on
        set(gca,'LineWidth',2,'FontSize',20)
        xlim([-5 15])
end

%% Iteration begins..
while norm(r0) > error_tolerance && iteration < 100
    iteration = iteration +1;
    hypoSta0 =- sta+(diag(hypoCenter0)*...
         ones(size(sta)));
    dist0 = sqrt(sum(abs(sta-(diag(hypoCenter0)*...
         ones(size(sta)))).^2,1))';     
    tCal = dist0/Vel+t0;
    r0 = tObs-tCal;
    G =[1/Vel * hypoSta0'./[dist0 dist0 dist0] ones(10,1)];%Updating Gmatrix
    [U,S,V] = svd(G,'econ');
    X = V*inv(S)*U'*r0;
    [hypoCenter0 ] = [hypoCenter0 ]+X(1:3);
    evolution_path(:,end+1) = hypoCenter0;
    t0 = t0 +X(4);
%     if Visualize
%        
%         scatter(hypoCenter0(1),hypoCenter0(2),100,'d','MarkerEdgeColor','green')
%         plot(evolution_path(1,:),evolution_path(2,:),'red','LineWidth',2);
%         axis equal
%         box on
%         set(gca,'LineWidth',2,'FontSize',20)
%         xlim([-5 15])
%         xlabel('X (km)');ylabel('Y (km)')
%         drawnow
%     end

end


  if Visualize
       
        scatter(evolution_path(1,:),evolution_path(2,:),100,'d','MarkerEdgeColor','green')
        plot(evolution_path(1,:),evolution_path(2,:),'red','LineWidth',2);
        axis equal
        box on
        set(gca,'LineWidth',2,'FontSize',20)
        xlim([-5 15])
        xlabel('X (km)');ylabel('Y (km)')
    end

ndf = 10-4;
sigma_2 = sum(r0.^2)/ndf;
covarianceMatrix = sigma_2*inv(G'*G);
c = 1.96;% 95% normal distribution 
invervalHypoUp = hypoCenter0+c*diag(covarianceMatrix(1:3,1:3));%Upper bound
invervalHypoLp = hypoCenter0-c*diag(covarianceMatrix(1:3,1:3));%Lower bound


