clear;
close all;
clc;

addpath(genpath('./map'));

% Setup parameter
Init_Parameter;

% create Landmark and Waypoints
[LM, wp] = create_LM_waypoints;

%Landmarkの位置 [x, y]
LM = LM';

[plines, pcount, dtsum, ftag, da_table, iwp, Omega, data] = Init_other_parameter(LM, xEst, xTrue, PEst);
 
%------------------------------------------------------------------------------------
% MAIN LOOP
%------------------------------------------------------------------------------------
tic;
while iwp ~= 0
    %--------------------------------------------------------
    % Control Input
    %--------------------------------------------------------
    % Waypointsの従った角速度を計算
    [Omega, iwp] = compute_steering(xTrue, wp, iwp, AT_WAYPOINT, Omega, RATEG, MAXG, dt);

    % perform loops: if final waypoint reached, go back to first
    if iwp==0 && NUMBER_LOOPS > 1, iwp=1; NUMBER_LOOPS= NUMBER_LOOPS-1; end
    
    % 速度と角速度
    u = [V Omega]';
    %--------------------------------------------------------
    % EKF-SLAM
    %--------------------------------------------------------
    % Observation
    [z,xTrue,xd,u] = Observation(xTrue, xd, u, LM, MAX_RANGE, Qsigma, Rsigma);

    % EKF-SLAM
%     [xEst,PEst,jj] = EKF_SLAM(u,z,xEst,PEst,initP,alpha,jj); 
    
    dtsum= dtsum + dt;
    if dtsum >= DT_OBSERVE
        dtsum= 0;
    end
    %--------------------------------------------------------
    % SEIF-SLAM
    %--------------------------------------------------------
    % SEIF-SLAM(example1.mapに対してのみ有効, MAXRANGE = 25m)
    [xEst,invPEst,Xi,j] = SEIF_SLAM_MK(u,z,xEst,invPEst,Xi,initinvPEst,initXi,alpha,j); 
    PEst = inv(invPEst);
%     [xEst,invPEst,Xi,j] = SEIF_SLAM_cubic(u,z,xEst,invPEst,Xi,initinvPEst,initXi,alpha,j); 
%     PEst = inv(invPEst);

    %--------------------------------------------------------
    % SAVE DATA
    %--------------------------------------------------------
    %Simulation Result
%     result.time=[result.time; time];
    result.xTrue=[result.xTrue; xTrue'];
    result.xd=[result.xd; xd'];
    result.xEst=[result.xEst;xEst(1:3)'];
    result.u=[result.u; u'];
    
    % offline data store
%     data = store_data(data, xEst, PEst, xTrue);
    
%     %pcov = make_covariance_ellipses(xEst, PEst);
%     ptmp= make_covariance_ellipses(xEst(1:3),PEst(1:3,1:3));
%     pcov(:,1:size(ptmp,2))= ptmp;
%     if dtsum==0
%         pcount= pcount+1;
%         if pcount == 100
%             gcf = figure(1);
%              plot(pcov(1,:), pcov(2,:), 'm'); hold on;
%              plot(data.xTrue(1,1:data.i), data.xTrue(2,1:data.i), 'r--', 'LineWidth', 1.5);
%              plot(data.xEst(1,1:data.i), data.xEst(2,1:data.i), 'g--', 'LineWidth', 1.5);
%             % Auto save graph
% %             saveas(gcf, 'EKF_SLAM.png');
% %             figure(2);
% %             imagesc(PEst)
% %             colorbar
%             pcount=0;
%         end
%         
% %         legend('LM','Trajectory', 'Waypoints', 'covariance','True', 'Estimate', 'Location', 'best');
%     end
%     drawnow
    
    %--------------------------------------------------------
    % Animation
    %--------------------------------------------------------
    pcount= pcount+1;
    if pcount == 200
        Animation(result,xTrue,LM,z,xEst, PEst);
        pcount=0;
    end
    if dtsum == 0
        % 状態推定誤差共分散行列の描画
        covstate = make_covariance_ellipses(xEst,PEst);
        plot(covstate(1,:),covstate(2,:), 'm'); 
        drawnow;
    end
end
toc
DrawGraph(result,xEst,LM);

function [plines, pcount, dtsum, ftag, da_table, iwp, Omega, data] = Init_other_parameter(lm, xEst, xTrue, PEst)
    % % initialise other variables and constants
    plines   = []; % for laser line animation
    pcount   = 0;
    dtsum    = 0; % change in time since last observation
    ftag     = 1:size(lm, 2); % identifier for each landmark
    da_table = zeros(1, size(lm, 2)); % data association table 
    data     = initialise_store(xEst, PEst, xTrue); % stored data for off-line
    iwp      = 1; % index to first waypoint 
    Omega    = 0; % initial steer angle
end
    % 初期データの格納
function data= initialise_store(xEst, PEst, xTrue)
    % offline storage initialisation
    data.i             = 1;
    data.xEst          = xEst;
    data.xTrue         = xTrue;
    data.state(1).xEst = xEst;
    data.state(1).PEst = diag(PEst);
end
    % シミュレーションのデータ格納
function data= store_data(data, xEst, PEst, xTrue)
    % add current data to offline storage
    CHUNK= 5000;
    if data.i == size(data.xEst,2) % grow array in chunks to amortise reallocation
        data.xEst  = [data.xEst zeros(3,CHUNK)];
        data.xTrue = [data.xTrue zeros(3,CHUNK)];
    end
    i                  = data.i + 1;
    data.i             = i;
    data.xEst(:,i)     = xEst(1:3);
    data.xTrue(:,i)    = xTrue;
    data.state(i).xEst = xEst;
    data.state(i).PEst = diag(PEst);
end

function p = make_covariance_ellipses(x,P)
    % compute ellipses for plotting state covariances
    N= 10;
    inc= 2*pi/N;
    phi= 0:inc:2*pi;

    lenx= length(x);
    lenf= (lenx-3)/2;
    p= zeros (2,(lenf+1)*(N+2));

    ii=1:N+2;
    p(:,ii)= make_ellipse(x(1:2), P(1:2,1:2), 2, phi);

    ctr= N+3;
    for i=1:lenf
        ii= ctr:(ctr+N+1);
        jj= 2+2*i; jj= jj:jj+1;

        p(:,ii)= make_ellipse(x(jj), P(jj,jj), 2, phi);
        ctr= ctr+N+2;
    end
end

function p= make_ellipse(x,P,s, phi)
    % make a single 2-D ellipse of s-sigmas over phi angle intervals 
    r= sqrtm(P);
    a= s*r*[cos(phi); sin(phi)];
    p(2,:)= [a(2,:)+x(2) NaN];
    p(1,:)= [a(1,:)+x(1) NaN];
end

function n=GetnLM(xEst)
    %ランドマークの数を計算する関数
    n=(length(xEst)-3)/2;
end

function Animation(result,xTrue,LM,z,xEst, PEst)
    %アニメーションを描画する関数
    hold off;
    plot(result.xTrue(:,1),result.xTrue(:,2),'.b');hold on;
    plot(LM(:,1),LM(:,2),'pk','MarkerSize',10);hold on;
    %観測線の表示
    if~isempty(z)
        for iz=1:length(z(:,1))
            ray=[xTrue(1:2)';z(iz,3:4)];
            plot(ray(:,1),ray(:,2),'-g');hold on;
        end
    end

    %SLAMの地図の表示
    for il=1:GetnLM(xEst)
        plot(xEst(4+2*(il-1)),xEst(5+2*(il-1)),'.c');hold on;
    end
    % plot(zl(1,:),zl(2,:),'.b');hold on;
    plot(result.xd(:,1),result.xd(:,2),'.k');hold on;
    plot(result.xEst(:,1),result.xEst(:,2),'.r');hold on;
    arrow=0.5;
    x=result.xEst(end,:);
    quiver(x(1),x(2),arrow*cos(x(3)),arrow*sin(x(3)),'ok');hold on;
    axis equal;
    grid on;

    drawnow;
end

function x = f(x, u)
    % Motion Model
    global dt;
    global PoseSize;
    global LMSize;

    F = horzcat(eye(PoseSize),zeros(PoseSize,LMSize*GetnLM(x)));

    B = [dt*cos(x(3)) 0
         dt*sin(x(3)) 0
         0 dt];

    x = x+F'*B*u;
    x(3) = PI2PI(x(3));%角度補正
end

function [z, x, xd, u] = Observation(x, xd, u, LM ,MAX_RANGE, Qsigma, Rsigma)
    %Calc Observation from noise prameter
    x  = f(x, u);% Ground Truth
    u  = u + Qsigma * randn(2,1);%add Process Noise
    xd = f(xd, u);% Dead Reckoning
    %Simulate Observation
    z = [];
    for iz = 1:length(LM(:, 1))
        % 2021 1/31
        %LMの位置をロボット座標系に変換
        yaw = zeros(3,1);
        yaw(3) = -x(3);
        localLM = HomogeneousTransformation2D(LM(iz,:)-x(1:2)',yaw');
        dx  = localLM(1) - x(1);
        dy  = localLM(2) - x(2);
        d   = norm(localLM);%距離
        angle = PI2PI(atan2(localLM(2),localLM(1)));
    %     if abs(pi / 2 - angle) < pi && d < MAX_RANGE
        if d < MAX_RANGE
            noise = Rsigma * randn(2,1);
            z = [z;[d+noise(1) PI2PI(atan2(localLM(2),localLM(1))) + noise(2) LM(iz,:)]];
        end
    end
end

function [lm, wp] = create_LM_waypoints()
    % 事前に設定した経路計画に沿って走行するルートを指定
    load('example1.mat');
    fig = figure(1);
    plot(lm(1,:),lm(2,:),'b*')
    hold on, axis equal
    plot(wp(1,:),wp(2,:), 'k', wp(1,:),wp(2,:),'k.')
    xlabel('X[m]'), ylabel('Y[m]')
    set(fig, 'name', 'EKF-SLAM Simulator')
end

function DrawGraph(result,xEst,LM)
    %Plot Result
    figure(1);
    hold off;
    x=[ result.xTrue(:,1:2) result.xEst(:,1:2)];
    set(gca, 'fontsize', 16, 'fontname', 'times');
    plot(x(:,1), x(:,2),'-b','linewidth', 4); hold on;
    plot(result.xd(:,1), result.xd(:,2),'-k','linewidth', 4); hold on;
    plot(x(:,3), x(:,4),'-r','linewidth', 4); hold on;
    plot(LM(:,1),LM(:,2),'pk','MarkerSize',10);hold on;%真のランドマークの位置
    %LMの地図の表示
    for il=1:GetnLM(xEst)
        plot(xEst(4+2*(il-1)),xEst(5+2*(il-1)),'.g','MarkerSize',10);hold on;
    end

    title('EKF SLAM Result', 'fontsize', 16, 'fontname', 'times');
    xlabel('X (m)', 'fontsize', 16, 'fontname', 'times');
    ylabel('Y (m)', 'fontsize', 16, 'fontname', 'times');
    legend('Ground Truth','Dead Reckoning','EKF SLAM','True LM','Estimated LM');
    grid on;
    axis equal;
end

function angle=PI2PI(angle)
    %ロボットの角度を-pi~piの範囲に補正する関数
    angle = mod(angle, 2*pi);

    i = find(angle>pi);
    angle(i) = angle(i) - 2*pi;

    i = find(angle<-pi);
    angle(i) = angle(i) + 2*pi;
end

function out = HomogeneousTransformation2D(in, base, mode)
%function out = HomogeneousTransformation2D(in, base,mode)
%HOMOGENEOUSTRANSFORMATION2D 二次元同次変換関数
%   一個の点の変換から，複数の点の一括変換まで可能．
%   レーザやレーダの点群を座標変換するのに便利です．
%   引数をinとbaseだけを入れた場合，各点を回転した後並進する．
%
%   Input1:変換前ベクトル [x_in_1 y_in_1;
%                        x_in_2  y_in_2;
%                               ....]
%           *変換前ベクトルは，３つめ以上の要素が含まれていても，
%            最初の２つを取り出し，残りは無視します．
%   Input2:基準ベクトル(トラックの位置とか) [x_base y_base theta_base]
%   Input3:同次変換モード　この変数は引数として入れなくても動きます．
%           引数を入れない場合，デフォルトで各点を回転した後並進する．
%           mode=0の場合も各点を回転した後並進する．
%           mode=1の場合，並進した後に回転します．
%
%   Output1:変換後ベクトル [x_out y_out;
%                        x_out_2  y_out_2;
%                               ....]

    %回転行列
    Rot=[cos(base(3)) sin(base(3)); -sin(base(3)) cos(base(3))];

    %点数分だけbase座標を配列に格納
    Nin=size(in);
    baseMat=repmat(base(1:2),Nin(1),1);

    % x-y以外のデータが入っていた場合．一回他の変数に置いておく．
    if Nin(2)>=3
        inxy=in(:,1:2);
        inOther=in(:,3:end);
        in=inxy;
    end

    %同次変換
    if nargin==2 || mode==0 %回転→並進
        out=baseMat+in*Rot;
    else %並進→回転
        out=(baseMat+in)*Rot;
    end

    %取り除いた値をくっつける．
    if Nin(2)>=3
        out=[out inOther];
    end
end

function radian = toRadian(degree)
    % degree to radian
    radian = degree/180*pi;
end

function degree = toDegree(radian)
    % radian to degree
    degree = radian/pi*180;
end