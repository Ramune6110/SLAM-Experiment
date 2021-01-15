function [xEst,PEst,zl] = HFSLAM(xEst,PEst,z,u)
global LMSize dt R Q alpha initP gamma delta

if isempty(z)
    F = horzcat(eye(3),zeros(3,LMSize*GetnLM(xEst))); 
    B = [dt*cos(xEst(3))  0
         dt*sin(xEst(3))  0
              0       dt]; 
    xEst = xEst + F'*B*u;
    xEst(3) = PI2PI(xEst(3));                                              
    disp('観測なし')
else
    disp('観測あり')
    % ------ EKF SLAM --------
    % Predict
    xEst = f(xEst, u);
    [G,Fx]=jacobF(xEst, u);
    PEst= G'*PEst*G + Fx'*R*Fx;
    % Update
    for iz=1:length(z(:,1))%それぞれの観測値に対して観測値をランドマークとして追加
        zl=CalcLMPosiFromZ(xEst,z(iz,:));%観測値そのものからLMの位置を計算
        xAug=[xEst;zl];%状態ベクトルと共分散行列の追加
        PAug=[PEst zeros(length(xEst),LMSize);
              zeros(LMSize,length(xEst)) initP];
        
        mdist=[];%マハラノビス距離のリスト
        for il=1:GetnLM(xAug) %それぞれのランドマークについて
            if il==GetnLM(xAug)
                mdist=[mdist alpha];%新しく追加した点の距離はパラメータ値を使う
            else
                lm=xAug(4+2*(il-1):5+2*(il-1));
                [y,S,~]=CalcInnovation(lm,xAug,PAug,z(iz,1:2),il);
                mdist=[mdist y'/S*y];%マハラノビス距離の計算
            end
        end
        
        %マハラノビス距離が最も近いものに対応付け
        [~,I]=min(mdist);
      
        %一番距離が小さいものが追加したものならば、その観測値をランドマークとして採用
        if I==GetnLM(xAug)
            %disp('New LM')
            xEst=xAug;
            PEst=PAug;
        end
        
        lm=xEst(4+2*(I-1):5+2*(I-1));%対応付けられたランドマークデータの取得
        %イノベーションの計算
        [y,S,H]=CalcInnovation(lm,xEst,PEst,z(iz,1:2),I);
        PS = (1 + delta)*(eye(length(xEst)) + H'/Q*H*PEst - gamma^(-2)*eye(length(xEst))*PEst);
        K = PEst/PS*H'/Q;
        xEst = xEst + K*y;
        PEst = (eye(size(xEst,1)) - K*H)*PEst;   
    end
end
    xEst(3)=PI2PI(xEst(3));%角度補正
end

function [y,S,H]=CalcInnovation(lm,xEst,PEst,z,LMId)
%対応付け結果からイノベーションを計算する関数
global Q;
delta=lm-xEst(1:2);
q=delta'*delta;
zangle=atan2(delta(2),delta(1))-xEst(3);
zp=[sqrt(q) PI2PI(zangle)];%観測値の予測
y=(z-zp)';
H=jacobH(q,delta,xEst,LMId);
S=H*PEst*H'+Q;
end

function n=GetnLM(xEst)
%ランドマークの数を計算する関数
n=(length(xEst)-3)/2;
end

function zl=CalcLMPosiFromZ(x,z)
%観測値からLMの位置を計算する関数
zl=x(1:2)+[z(1)*cos(x(3)+z(2));z(1)*sin(x(3)+z(2))];
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

x= x+F'*B*u;
x(3)=PI2PI(x(3));%角度補正
end

function [G,Fx]=jacobF(x, u)
% 運動モデルのヤコビ行列の計算関数
global dt;
global PoseSize;
global LMSize;

Fx = horzcat(eye(PoseSize),zeros(PoseSize,LMSize*GetnLM(x)));
 
jF=[0 0 -dt*u(1)*sin(x(3))
    0 0 dt*u(1)*cos(x(3))
    0 0 0];

G=eye(length(x))+Fx'*jF*Fx;
end

function H=jacobH(q,delta,x,i)
%観測モデルのヤコビ行列を計算する関数
sq=sqrt(q);
G=[-sq*delta(1) -sq*delta(2) 0 sq*delta(1) sq*delta(2);
    delta(2)    -delta(1)   -1 -delta(2)    delta(1)];
G=G/q;
F=[eye(3) zeros(3,2*GetnLM(x));
   zeros(2,3) zeros(2,2*(i-1)) eye(2) zeros(2,2*GetnLM(x)-2*i)];
H=G*F;
end

function angle=PI2PI(angle)
%ロボットの角度を-pi~piの範囲に補正する関数
angle = mod(angle, 2*pi);

i = find(angle>pi);
angle(i) = angle(i) - 2*pi;

i = find(angle<-pi);
angle(i) = angle(i) + 2*pi;
end