function [xhat,P,NL,y,jx] = EKF_SLAM(U,LMo,xhat,P,initP,alpha,jx)
% ------ EKF-SLAM --------
global LMSize;    
global R;
global dt;
% global Rsigma;
% global LM;
        NL = NumLM(xhat);
%         LMo
if isempty(LMo)
    F = horzcat(eye(3),zeros(3,LMSize*NumLM(xhat))); 
    B = [dt*cos(xhat(3))  0
         dt*sin(xhat(3))  0
              0       dt]; 
    xhat = xhat + F'*B*U;
%     xhat(3) = PI2PI(xhat(3));                                              %角度補正
y = [0; 0];
% disp('観測なし')
else
% disp('観測あり')
    % Predict
    xhat = f(xhat, U);
%     xhat(3) = PI2PI(xhat(3));                                              %角度補正
    [G, Fx] = jacobF(xhat, U);
    P = G'*P*G + Fx'*R*Fx;        
    for ilm = 1: length(LMo(:,1))
        zl=CalcLMPosiFromZ(xhat,LMo(ilm,:));                               %観測値そのものからLMの位置を計算
        xAug=[xhat;zl];
        PAug=[P zeros(length(xhat),LMSize);
              zeros(LMSize,length(xhat)) initP];
        Jx = [jx; LMo(ilm,3)]; 
        jM = [];                                                           % 特徴量格納
        for jrow = 1 : length(Jx)                                          % 既知の特徴量による判別
          if jrow == length(Jx)
             jM = [jM 10^(-1)];
          else
             D = norm(Jx(jrow,:)-LMo(ilm,3));
             jM = [jM D];
          end
        end                                                                % 特徴量の判別終了                                                             
        [mi,num] = min(jM);                                                % ノルムの最小値と最小値の行番号の導出
        if num == length(Jx)                                               % 初検知したランドマークの場合
             jx = Jx;                                                      % 特徴量を反映
             xhat=xAug;                                                    % 状態量を反映
             P=PAug;                                                       % 共分散行列を反映
        end  
    lm=xhat(4+2*(num-1):5+2*(num-1));                                          %対応付けられたランドマークデータの取得
   
    %innovation and update
   [y,S,H]=CalcInnovation(lm,xhat,P,LMo(ilm,1:2),num);
    K = P*H' / S;
    xhat = xhat + K*y;
%     xhat(3) = PI2PI(xhat(3));
    P = (eye(size(xhat,1)) - K*H)*P;
    end
%     xhat(3) = PI2PI(xhat(3));                                              %角度補正
end
    
function [y,S,H]=CalcInnovation(lm,xhat,P,z,LMId)
%対応付け結果からイノベーションを計算する関数
global Q;
delta=lm-xhat(1:2);                                                        %dx,dyの導出(絶対座標系)
q=delta'*delta;                                                            %dx^2,dy^2の導出
zangle=atan2(delta(2),delta(1))-xhat(3);
zp=[sqrt(q) PI2PI(zangle)];                                                %観測値の予測
y=(z-zp)';
H=jacobH(q,delta,xhat,LMId);
S=H*P*H'+Q;

function n = NumLM(xhat)
n = (length(xhat)-3)/2;

function zl=CalcLMPosiFromZ(x,z)
zl=x(1:2)+[z(1)*cos(x(3)+z(2));z(1)*sin(x(3)+z(2))];                       %観測値からLMの位置を計算する関数

function x = f(x, u)
global dt;
global LMSize;

F = horzcat(eye(3),zeros(3,LMSize*NumLM(x))); 
B = [dt*cos(x(3))  0
     dt*sin(x(3))  0
          0       dt]; 
x = x + F'*B*u;
% x(3)=PI2PI(x(3));%角度補正

function H=jacobH(q,delta,x,i)
%観測モデルのヤコビ行列を計算する関数
sq=sqrt(q);
G=[-sq*delta(1) -sq*delta(2) 0 sq*delta(1) sq*delta(2);
    delta(2)    -delta(1)   -1 -delta(2)    delta(1)];
G=G/q;
F=[eye(3) zeros(3,2*NumLM(x));
   zeros(2,3) zeros(2,2*(i-1)) eye(2) zeros(2,2*NumLM(x)-2*i)];
H=G*F;

function [G, Fx] = jacobF(x, u)
global dt;
global LMSize;

Fx = horzcat(eye(3),zeros(3,LMSize*NumLM(x))); 
JF = [0 0 -dt*u(1)*sin(x(3));
      0 0  dt*u(1)*cos(x(3));
      0 0         0         ];
G = eye(length(x)) + Fx'*JF*Fx;

function angle=PI2PI(angle)
%ロボットの角度を-pi~piの範囲に補正する関数
angle = mod(angle, 2*pi);

i = find(angle>pi);
angle(i) = angle(i) - 2*pi;

i = find(angle<=-pi);
angle(i) = angle(i) + 2*pi;