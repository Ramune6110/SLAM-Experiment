% シミュレーション時間設定
time = 0;
global endtime;
endtime = 500; % シミュレーション終了時間[sec]
global dt;
dt = 0.5; % シミュレーション刻み時間[sec]
nSteps = ceil((endtime - time)/dt);%シミュレーションのステップ数
 
%計算結果格納用変数
result.time=[];
result.xTrue=[];
result.xEst_ekf=[];
result.xEst_hf=[];
result.u = [];

% State Vector [x y yaw]'
global START
START = [0;0;0];
xEst_ekf=[0;0;0];
xEst_hf=[0;0;0];
global PoseSize;
PoseSize=3;%ロボットの姿勢の状態数[x,y,yaw]
global LMSize;
LMSize=2;%ランドマークの状態量[x,y]
% True State
xTrue=[0;0;0];
 
% % Covariance Matrix for predict
% global R
% R=diag([0.2 0.2 toRadian(1)]).^2;
%  
% % Covariance Matrix for observation
% global Q;
% Q=diag([10 toRadian(30)]).^2;%range[m], Angle[rad]

% Simulation parameter
global Qsigma
Qsigma=diag([0.1 toRadian(20)]).^2;
global Rsigma
Rsigma=diag([0.1 toRadian(1)]).^2;

%ゴールの位置
global GOAL
%GOAL=[40 35];
GOAL=[25 15];

%Landmarkの位置 [x, y]
global LM;
% LM=[0 18;
%     0 10;
%     10 5;
%     20 10;
%     10 20;
%     15 30;
%     20 35;
%     25 20;
%     30 30;
%     30 40;
%     30 35;
%     35 25;
%     35 35;
%     40 20];

% LM=[10 20;
%     10 18;
%     10 16;
%     10 14;
%     10 12;
%     10 10;
%     12 10;
%     14 10;
%     16 10;
%     18 10;
%     20 10];

% LM=[10 20;
%     10 19;
%     10 18;
%     10 17;
%     10 16;
%     10 15;
%     10 14;
%     10 13;
%     10 12;
%     10 11;
%     10 10;
%     11 10;
%     12 10;
%     13 10;
%     14 10;
%     15 10;
%     16 10;
%     17 10;
%     18 10;
%     19 10;
%     20 10];

obstacleR=0.5;%衝突判定用の障害物の半径
MAX_RANGE=10;%最大観測距離
global alpha
alpha=0.01;%ランドマーク識別用マハラノビス距離閾値
global initP
PEst_ekf = eye(3);
PEst_hf = eye(3);
initP=eye(2)*100000;

%ロボットの力学モデル
%[最高速度[m/s],最高回頭速度[rad/s],最高加減速度[m/ss],最高加減回頭速度[rad/ss],
% 速度解像度[m/s],回頭速度解像度[rad/s]]
Kinematic=[3.0,toRadian(40.0),1.0,toRadian(50.0),0.01,toRadian(1)];

%評価関数のパラメータ [heading,dist,velocity,predictDT]
evalParam=[0.1,0.2,0.1,3.0];
%evalParam=[1.0,2.0,1.0,5.0];

%初期化
x = [START(1);START(2);START(3);0;0];

global gamma delta
gamma = ((1*sqrt(80))^(2));
%gamma = 10;
% gamma = ((1*sqrt(85))^(2));
% gamma = ((1*sqrt(90))^(2));
% gamma = ((1*sqrt(100))^(2));
% gamma = ((1*sqrt(150))^(2));
% gamma = (18^(-2));
delta = 10^(-1);

function radian = toRadian(degree)
% degree to radian
radian = degree/180*pi;
end