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
%     xhat(3) = PI2PI(xhat(3));                                              %�p�x�␳
y = [0; 0];
% disp('�ϑ��Ȃ�')
else
% disp('�ϑ�����')
    % Predict
    xhat = f(xhat, U);
%     xhat(3) = PI2PI(xhat(3));                                              %�p�x�␳
    [G, Fx] = jacobF(xhat, U);
    P = G'*P*G + Fx'*R*Fx;        
    for ilm = 1: length(LMo(:,1))
        zl=CalcLMPosiFromZ(xhat,LMo(ilm,:));                               %�ϑ��l���̂��̂���LM�̈ʒu���v�Z
        xAug=[xhat;zl];
        PAug=[P zeros(length(xhat),LMSize);
              zeros(LMSize,length(xhat)) initP];
        Jx = [jx; LMo(ilm,3)]; 
        jM = [];                                                           % �����ʊi�[
        for jrow = 1 : length(Jx)                                          % ���m�̓����ʂɂ�锻��
          if jrow == length(Jx)
             jM = [jM 10^(-1)];
          else
             D = norm(Jx(jrow,:)-LMo(ilm,3));
             jM = [jM D];
          end
        end                                                                % �����ʂ̔��ʏI��                                                             
        [mi,num] = min(jM);                                                % �m�����̍ŏ��l�ƍŏ��l�̍s�ԍ��̓��o
        if num == length(Jx)                                               % �����m���������h�}�[�N�̏ꍇ
             jx = Jx;                                                      % �����ʂ𔽉f
             xhat=xAug;                                                    % ��ԗʂ𔽉f
             P=PAug;                                                       % �����U�s��𔽉f
        end  
    lm=xhat(4+2*(num-1):5+2*(num-1));                                          %�Ή��t����ꂽ�����h�}�[�N�f�[�^�̎擾
   
    %innovation and update
   [y,S,H]=CalcInnovation(lm,xhat,P,LMo(ilm,1:2),num);
    K = P*H' / S;
    xhat = xhat + K*y;
%     xhat(3) = PI2PI(xhat(3));
    P = (eye(size(xhat,1)) - K*H)*P;
    end
    xhat(3) = PI2PI(xhat(3));                                              %�p�x�␳
end
    
function [y,S,H]=CalcInnovation(lm,xhat,P,z,LMId)
%�Ή��t�����ʂ���C�m�x�[�V�������v�Z����֐�
global Q;
delta=lm-xhat(1:2);                                                        %dx,dy�̓��o(��΍��W�n)
q=delta'*delta;                                                            %dx^2,dy^2�̓��o
zangle=atan2(delta(2),delta(1))-xhat(3);
zp=[sqrt(q) PI2PI(zangle)];                                                %�ϑ��l�̗\��
y=(z-zp)';
H=jacobH(q,delta,xhat,LMId);
S=H*P*H'+Q;

function n = NumLM(xhat)
n = (length(xhat)-3)/2;

function zl=CalcLMPosiFromZ(x,z)
zl = x(1:2)+[z(1)*cos(x(3)+z(2));z(1)*sin(x(3)+z(2))];                       %�ϑ��l����LM�̈ʒu���v�Z����֐�

function x = f(x, u)
global dt;
global LMSize;

F = horzcat(eye(3),zeros(3,LMSize*NumLM(x))); 
B = [dt*cos(x(3))  0
     dt*sin(x(3))  0
          0       dt]; 
x = x + F'*B*u;
x(3) = PI2PI(x(3));%�p�x�␳

function H=jacobH(q,delta,x,i)
%�ϑ����f���̃��R�r�s����v�Z����֐�
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
%���{�b�g�̊p�x��-pi~pi�͈̔͂ɕ␳����֐�
angle = mod(angle, 2*pi);

i = find(angle>pi);
angle(i) = angle(i) - 2*pi;

i = find(angle<=-pi);
angle(i) = angle(i) + 2*pi;