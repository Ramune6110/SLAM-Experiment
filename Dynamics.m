function u = Dynamics(u,V,omega)
global dt
b = 0.073;       %機体半幅
d = 0.01;        %CORとCOGの距離
l = 0.087;       %機体半長
M = 0.480;       %機体質量
I = 0.00594;     %機体慣性モーメント
g =9.81;         %重力加速度

Ka = 1.2;
Kb = 1;
Kc = 2;

%土壌係数
myuu = 0.2;
Rx = myuu*M*g;

%走行抵抗
Rfy = (1/2)*myuu*M*g;
fy = Rfy/l;
Mr = 2*fy*(l^2-d^2);

%Vd,faiはポテンシャル関数の指令値
%psi=現在の角度
%V=現在の並進速度
ad = Ka*(u(1) - V)/dt;
omegad = Kb*(fai-psi)/dt;
omegaA = Kc*(u(2) - omega)/dt;

%各クローラの速度の計算
v1 = V - b*omega;
v2 = V + b*omega;
 
%出すべき推進力 
F1 =  (1/2)*M*ad+Rx*sign(v1)-(I*omegaA+Mr*sign(omega))/(2*b);
F2 =  (1/2)*M*ad+Rx*sign(v2)+(I*omegaA+Mr*sign(omega))/(2*b);

%G1の計算
if v1 ~= 0
    G1 = F1-Rx*sign(v1);
elseif (v1 == 0) && (abs(F1)<Rx)
    G1 = 0;
elseif (v1 == 0) && (abs(F1)>Rx)
    G1 = F1-Rx*sign(F1);
end   

%G2の計算
if v2 ~= 0
    G2 = F2-Rx*sign(v2);
elseif (v2 == 0) && (abs(F2)<Rx)
    G2 = 0;
elseif (v2 == 0) && (abs(F2)>Rx)
    G2 = F2-Rx*sign(F2);
end 

%G3の計算
Mf = b*(G2-G1);

if omega ~= 0
  G3 = Mf-Mr*sign(omega);
elseif (omega == 0)&& (abs(Mf)<=Mr)
  G3 = 0;
elseif (omega == 0)&& (abs(Mf)>Mr)  
  G3 = Mf-Mr*sign(Mf);
end

%運動方程式
a = (1/M)*(G1+G2);
psid = (1/I)*G3;

%並進,回転速度の計算
V = V + a*dt;
omega = omega + psid*dt;
u = [V; omega];