function [z,x,u] = Observation(x, u, LM, MAX_RANGE,Noise11,Noise22)
global Qsigma Rsigma;
x=f(x, u);% Ground Truth
u=u+sqrt(Qsigma)*(Noise11);%add Process Noise
%Simulate Observation
z=[];
for iz=1:length(LM(:,1))
    %LMの位置をロボット座標系に変換
    yaw=zeros(3,1);
    yaw(3)=-x(3);
    localLM=HomogeneousTransformation2D(LM(iz,:)-x(1:2)',yaw');
    d=norm(localLM);%距離
    if d<MAX_RANGE %観測範囲内
        noise=sqrt(Rsigma)*(Noise22);
        z=[z;[d+noise(1) PI2PI(atan2(localLM(2),localLM(1))+noise(2)) LM(iz,:)]];
    end
end
end

function n=GetnLM(xEst)
%ランドマークの数を計算する関数
n=(length(xEst)-3)/2;
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

function out = HomogeneousTransformation2D(in, base, mode)
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

function angle=PI2PI(angle)
%ロボットの角度を-pi~piの範囲に補正する関数
angle = mod(angle, 2*pi);

i = find(angle>pi);
angle(i) = angle(i) - 2*pi;

i = find(angle<-pi);
angle(i) = angle(i) + 2*pi;
end