function [] = EKFSLAM_DWA()
clear;
close all;
clc;

Init % 数値初期化

Init_Parameter;

% create LM and Waypoints
[LM, wp] = create_LM_waypoints;

iwp = 1; % index to first waypoint 
Omega = 0; % initial steer angle

LM = LM';

tic;
% Main loop
for i=1 : nSteps
    % Waypointsの従った角速度を計算
    [Omega, iwp] = compute_steering(xTrue, wp, iwp, AT_WAYPOINT, Omega, RATEG, MAXG, dt);

    % perform loops: if final waypoint reached, go back to first
    if iwp==0 & NUMBER_LOOPS > 1, iwp=1; NUMBER_LOOPS= NUMBER_LOOPS-1; end
    
     u = [V Omega]';
% for i=1 : nSteps
    time = time + dt;
    % Input
%     u= DynamicWindowApproach(x,xTrue,Kinematic,GOAL,evalParam,LM,obstacleR);
    
    % Noise
    [Noise11, Noise22] = Noise(nSteps, i);
    % Observation
    [z,xTrue,u] = Observation(xTrue, u, LM, MAX_RANGE, Noise11, Noise22);
    % EKFSLAM
    [xEst_ekf, PEst_ekf, zl_ekf] = EKFSLAM(xEst_ekf, PEst_ekf, z, u);
    % HFSLAM
%     [xEst_hf, PEst_hf, zl_hf] = HFSLAM(xEst_hf, PEst_hf, z, u);

    %Simulation Result
    result.time=[result.time; time];
    result.xTrue=[result.xTrue; xTrue'];
    result.xEst_ekf=[result.xEst_ekf; xEst_ekf(1:3)'];
    result.xEst_hf=[result.xEst_hf; xEst_hf(1:3)'];
    result.u = [result.u; u'];
    
    % 10StepごとにAnimetion
    if rem(time,10)==0 
        Animation(result,xTrue,LM,z,xEst_ekf,zl_ekf);
    end
    
%     if sqrt((xTrue(1)-GOAL(1))^2+(xTrue(2)-GOAL(2))^2) < 0.5 %収束半径
%       disp('Arrive Goal!');
%       break;
%     end
    
end

toc

DrawGraph(result, xEst_ekf);

function Animation(result,xTrue,LM,z,xEst,zl)
%アニメーションを描画する関数
hold off;
plot(result.xTrue(:,1),result.xTrue(:,2),'.b');hold on;
plot(LM(:,1),LM(:,2),'pk','MarkerSize',10);hold on;
% %観測線の表示
% if~isempty(z)
%     for iz=1:length(z(:,1))
%         ray=[xTrue(1:2)';z(iz,3:4)];
%         plot(ray(:,1),ray(:,2),'-g');hold on;
%     end
% end
%SLAMの地図の表示
for il=1:GetnLM(xEst)
    plot(xEst(4+2*(il-1)),xEst(5+2*(il-1)),'r^');hold on;
end
plot(zl(1,:),zl(2,:),'.b');hold on;
plot(result.xEst_ekf(:,1),result.xEst_ekf(:,2),'.r');hold on;
plot(result.xEst_hf(:,1),result.xEst_hf(:,2),'.g');hold on;
arrow=0.5;
x=result.xTrue(end,:);
quiver(x(1),x(2),arrow*cos(x(3)),arrow*sin(x(3)),'ok');hold on;
axis equal;
grid on;
drawnow;

function [lm, wp] = create_LM_waypoints()
    % 事前に設定した経路計画に沿って走行するルートを指定
    %load('example_webmap.mat');
    load('example5.mat');
    
function DrawGraph(result, xEst)
figure(1);
hold off;
x=[result.xTrue(:,1:3) result.xEst_ekf(:,1:3) result.xEst_hf(:,1:3)];
global START GOAL LM
set(gca, 'fontsize', 16, 'fontname', 'times');
plot(x(:,1), x(:,2),'-b','linewidth', 0.5); hold on;
plot(x(:,4), x(:,5),'-r','linewidth', 0.5); hold on;
plot(LM(:,1),LM(:,2),'pk','MarkerSize',10);hold on;%真のランドマークの位置
for il=1:GetnLM(xEst)
   plot(xEst(4+2*(il-1)),xEst(5+2*(il-1)),'r^');hold on;
end
plot(START(1),START(2),'rs');hold on;
plot(GOAL(1),GOAL(2),'bs');hold on;
title('EKF SLAM With Potencial Function', 'fontsize', 16, 'fontname', 'times');
xlabel('X (m)', 'fontsize', 16, 'fontname', 'times');
ylabel('Y (m)', 'fontsize', 16, 'fontname', 'times');
legend('Ground Truth','EKF SLAM','True LM','Estimated LM','Location','Best');
grid on;
axis equal;

figure(2);
subplot(2,1,1)
plot(result.time, result.u(:,1), 'k');
xlabel('Time [s]');ylabel('Velocity');
grid on;

subplot(2,1,2)
plot(result.time, result.u(:,2), 'k');
xlabel('Time [s]');ylabel('Angular Velocity');
grid on;

RMSE = @(x) sqrt(mean(x.^2)); %平均二乗誤差
fprintf('%10s %10s\n','variable','RMSE(ekf)');
for p=1:3
    vname = sprintf('x%d',p);
    fprintf('%10s %10.5f \n',vname,RMSE(x(:,3+p)-x(:,p)));
    fprintf('%10s %10.5f \n',vname,RMSE(x(:,6+p)-x(:,p)));
end

OMEGA = @(x) sum((x-mean(x)).^2)/(length(x)-1);%分散
fprintf('%10s %10s\n','variable','DISPERSION(ekf)');
for p=1:3
    vname = sprintf('x%d',p);
    fprintf('%10s %10.5f \n',vname,OMEGA(x(:,3+p)-x(:,p)));
    fprintf('%10s %10.5f \n',vname,OMEGA(x(:,6+p)-x(:,p)));
end

function n=GetnLM(xEst)
%ランドマークの数を計算する関数
n=(length(xEst)-3)/2;