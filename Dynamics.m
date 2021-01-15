function u = Dynamics(u,V,omega)
global dt
b = 0.073;       %�@�̔���
d = 0.01;        %COR��COG�̋���
l = 0.087;       %�@�̔���
M = 0.480;       %�@�̎���
I = 0.00594;     %�@�̊������[�����g
g =9.81;         %�d�͉����x

Ka = 1.2;
Kb = 1;
Kc = 2;

%�y��W��
myuu = 0.2;
Rx = myuu*M*g;

%���s��R
Rfy = (1/2)*myuu*M*g;
fy = Rfy/l;
Mr = 2*fy*(l^2-d^2);

%Vd,fai�̓|�e���V�����֐��̎w�ߒl
%psi=���݂̊p�x
%V=���݂̕��i���x
ad = Ka*(u(1) - V)/dt;
omegad = Kb*(fai-psi)/dt;
omegaA = Kc*(u(2) - omega)/dt;

%�e�N���[���̑��x�̌v�Z
v1 = V - b*omega;
v2 = V + b*omega;
 
%�o���ׂ����i�� 
F1 =  (1/2)*M*ad+Rx*sign(v1)-(I*omegaA+Mr*sign(omega))/(2*b);
F2 =  (1/2)*M*ad+Rx*sign(v2)+(I*omegaA+Mr*sign(omega))/(2*b);

%G1�̌v�Z
if v1 ~= 0
    G1 = F1-Rx*sign(v1);
elseif (v1 == 0) && (abs(F1)<Rx)
    G1 = 0;
elseif (v1 == 0) && (abs(F1)>Rx)
    G1 = F1-Rx*sign(F1);
end   

%G2�̌v�Z
if v2 ~= 0
    G2 = F2-Rx*sign(v2);
elseif (v2 == 0) && (abs(F2)<Rx)
    G2 = 0;
elseif (v2 == 0) && (abs(F2)>Rx)
    G2 = F2-Rx*sign(F2);
end 

%G3�̌v�Z
Mf = b*(G2-G1);

if omega ~= 0
  G3 = Mf-Mr*sign(omega);
elseif (omega == 0)&& (abs(Mf)<=Mr)
  G3 = 0;
elseif (omega == 0)&& (abs(Mf)>Mr)  
  G3 = Mf-Mr*sign(Mf);
end

%�^��������
a = (1/M)*(G1+G2);
psid = (1/I)*G3;

%���i,��]���x�̌v�Z
V = V + a*dt;
omega = omega + psid*dt;
u = [V; omega];