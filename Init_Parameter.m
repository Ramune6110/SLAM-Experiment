%--------------------------------------------------------
% control parameters
%--------------------------------------------------------
V           = 3; % m/s
MAXG        = 30*pi/180; % radians, maximum steering angle (-MAXG < g < MAXG)
RATEG       = 20*pi/180; % rad/s, maximum rate of change in steer angle
WHEELBASE   = 4; % metres, vehicle wheel-base
DT_CONTROLS = 0.025; % seconds, time interval between control signals

%--------------------------------------------------------
% observation parameters
%--------------------------------------------------------
% MAX_RANGE=20;%�ő�ϑ�����
% MAX_RANGE  = 30;%�ő�ϑ�����(EKF-SLAM example1.map)
MAX_RANGE  = 20;%�ő�ϑ�����(SEIF-SLAM example1.map)
DT_OBSERVE = 8 * DT_CONTROLS; % seconds, time interval between observations

%--------------------------------------------------------
% waypoint proximity
%--------------------------------------------------------
AT_WAYPOINT  = 1.0; % metres, distance from current waypoint at which to switch to next waypoint
NUMBER_LOOPS = 2; % number of loops through the waypoint list

%--------------------------------------------------------
% save data box
%--------------------------------------------------------
result.time  = [];
result.xTrue = [];
result.xd    = [];
result.xEst  = [];
result.z     = [];
result.PEst  = [];
result.u     = [];

%--------------------------------------------------------
% state value
%--------------------------------------------------------
global dt;
dt = DT_CONTROLS; % �V�~�����[�V�������ݎ���[sec]
 
% State Vector [x y yaw]'
xEst = [0 0 0]';

% True State
xTrue = xEst;
 
% Dead Reckoning State
xd = xTrue;

global PoseSize;
PoseSize = length(xEst);%���{�b�g�̎p���̏�Ԑ�[x,y,yaw]

global LMSize;
LMSize = 2;%�����h�}�[�N�̏�ԗ�[x,y]

%--------------------------------------------------------
% Covariance Matrix
%--------------------------------------------------------
% Covariance Matrix for predict
global R;
R = diag([0.0001 0.0001 0.0001]).^2;
 
% Covariance Matrix for observation
global Q;
% Q = diag([0.3^2 (3.0 * pi / 180)^2]);%range[m], Angle[rad]
Q = diag([0.5 0.5]).^2; 

% Simulation parameter
Qsigma = diag([0.3 toRadian(3)]).^2;
Rsigma = diag([0.1 toRadian(1)]).^2;

%--------------------------------------------------------
% EKF-SLAM Parameter
%--------------------------------------------------------
% alpha=0.2;%�����h�}�[�N���ʗp�}�n���m�r�X����臒l
% alpha=0.0001;%�����h�}�[�N���ʗp�}�n���m�r�X����臒l
% alpha = 10^-9;%�����h�}�[�N���ʗp�}�n���m�r�X����臒l(example1.map)
% alpha = 3; %(EKF-SLAM example1.map)
alpha = 0.5; %(SEIF-SLAM example1.map)

PEst  = diag(ones(1,3)*10^-4);
initP = eye(2)*10^9;

jj = [];

%--------------------------------------------------------
% SEIF-SLAM Parameter
%--------------------------------------------------------
global Q1;
Q1 = inv(R);

global R1;
R1 = inv(Q);

j           = [];
invPEst     = inv(PEst);
initinvPEst = inv(initP);    
Xi          = PEst\xEst; % ���x�N�g��
initXi      = initP\[0;0]; % ���x�N�g��

function radian = toRadian(degree)
    % degree to radian
    radian = degree/180*pi;
end