function [xhat,Omega,Xi,jx] = SEIF_SLAM_cubic(U,z,xhat,Omega_chi,Xi_chi,initOmega,initXi,jx)
global R1;
global Q1;
global LMSize

    % Motion Update Step
    [G,JF,Fx,Delta] = jacobF(xhat, U);                                     % ���R�r�s��̌v�Z(��ԕ�����)
    Ps = Fx'*(inv(JF)-eye(3))*Fx;                                          % ���R�r�s��̋t�s����������邽�߂̕ϐ� G^-1 = I + Ps
    lambda = Ps'*Omega_chi + Omega_chi*Ps + Ps'*Omega_chi*Ps;              % P.358(�m�����{�e�B�N�X)
    Ph = Omega_chi + lambda;                                               % ���s��,�y��,EKF�̋����U�̗\���X�V�̑�1���ڂ̕ϐ�
    kappa = Ph*Fx'*inv(Q1 + Fx*Ph*Fx')*Fx*Ph;                              % ���s��,�y��,EKF�̋����U�̗\���X�V�̑�2���ڂ̕ϐ�
    Omega_bar = Ph - kappa;                                                % �\���X�V�ɂ�������s��
    Xi_bar = Xi_chi + (lambda - kappa)*xhat + Omega_bar*Fx'*Delta;         % �\���X�V�ɂ�������x�N�g��
    xhat_bar = f(xhat,U);                                                  % ��ԗ\��
        
    Xi = Xi_bar;
    Omega = Omega_bar;
    xhat = xhat_bar;                                                       % ��ԃx�N�g���̓��o
%         xhat = inv(Omega)*Xi;                                            % ��ԃx�N�g���̓��o
%     xhat(3) = PI2PI(xhat(3));                                              % �p�x�␳
    if isempty(z)~=1
    % Measurement Update Step
     numM = [];
     zM = [];
    for ILm = 1: length(z(:,1))                                            % ���m���������h�}�[�N���ׂČv�Z
        zl=CalcLMPosiFromZ(xhat,z(ILm,:));                                 % �ϑ��l���̂��̂���LM�̈ʒu���v�Z
        Jx = [jx; z(ILm,3)];                                               % �����ʒǉ�
        Xhatbar = [xhat; zl];                                              % �����h�}�[�N�ǉ�
        Xibar = [Xi; initXi];                                              % ���x�N�g���ǉ�(�����h�}�[�N�ǉ��ɂ��킹��)
        Omegabar=[Omega zeros(length(Omega),LMSize);
                  zeros(LMSize,length(Omega)) initOmega];                 % ���s��ǉ�(�����h�}�[�N�ǉ��ɂ��킹��)              
        jM = [];                                                           % �����ʊi�[

        for jrow = 1 : length(Jx)                                          % ���m�̓����ʂɂ�锻��
          if jrow == length(Jx)
             jM = [jM 10^(-1)];
          else
             D = norm(Jx(jrow,:)-z(ILm,3));
             jM = [jM D];
          end
        end                                                                % �����ʂ̔��ʏI��                                                             
        [mi,num] = min(jM);                                                % �m�����̍ŏ��l�ƍŏ��l�̍s�ԍ��̓��o
                                                                           % �ǉ����������h�}�[�N�ɂ�����V���Ȓ�`
       numM=[numM;num];
       zM = [zM;z(ILm,3)];
          if num == length(Jx)                                             % �����m���������h�}�[�N�̏ꍇ
             jx = Jx;                                                      % �����ʂ𔽉f
             xhat = Xhatbar;                                               % ��ԗʂ𔽉f
             Xi = Xibar;                                                   % ���x�N�g���𔽉f
             Omega = Omegabar;                                             % ���s��𔽉f
          end                                                              % ��`�̏I��
    lm = xhat(4+2*(num-1):5+2*(num-1));                                    % �Ή��t����ꂽ�����h�}�[�N�f�[�^�̎擾
    delta = lm - xhat(1:2);                                                % dx,dy�̓��o(��΍��W�n)
    q = delta'*delta;                                                      % dx^2,dy^2�̓��o
    zhat = [sqrt(q) PI2PI(atan2(delta(2),delta(1))-xhat(3))];              % �ϑ��l�̗\��
    y = (z(ILm,1:2) - zhat)';                                              % �C�m�x�[�V����
    H = jacobH(q,delta,xhat,num);                                          % ���R�r�s��v�Z(�ϑ�������)
    Xi = Xi + H'*R1*(y + H*xhat);                                          % �v���X�V(���x�N�g��)
    Omega = Omega + H'*R1*H;                                               % �v���X�V(���s��)
    end
    
    if NumLM(xhat) < 4                                                     % �a���X�V���s��Ȃ��ꍇ
    Omegaxx = Omega(1:3,1:3);
    Omegaxm = Omega(1:3,4:length(Omega));
    Xhatlm = xhat(4:length(xhat));
    Xix = Xi(1:3);
    Xhats = inv(Omegaxx)*(Xix-(Omegaxm*Xhatlm));
    xhat = [Xhats;Xhatlm];                                                 % ��ԃx�N�g���̓��o
%     xhat(3) = PI2PI(xhat(3));                                              % �p�x�␳
    else                                                                   % �a���X�V���s���ꍇ
    % Sparsification Step
    FF = eye(length(xhat(:,1)));
    Fx = FF(1:3,:);
    FmplusM = [];
    Fm0M = [];
    E = 1 : length(jx);
    for ii = 1 : length(jx)
        CC=0;
        for iii = 1 : length(numM)
            if jx(ii,1)==zM(iii,1)
               il=numM(iii,1);
            else
               CC = CC +1;
            end
            if CC == length(numM)
               ill=E(1,ii);
               il = [];
            else
                ill = [];
            end
        end
               FmplusM = [FmplusM;FF(4+2*(il-1):5+2*(il-1),:)];
            if isempty(ill)~=1
               Fm0M=[Fm0M;FF(4+2*(ill-1):5+2*(ill-1),:)];
            end
    end
    Fm0=Fm0M;
    Fmplus=FmplusM;

    Fx_m0 = [Fx; Fm0];                                                     % projection matrix active to passive(m0)
    F_x_mplus_m0 = [Fx; Fm0; Fmplus];                                       % projection matrix active(m+)
    
    Omega0 = (F_x_mplus_m0')*F_x_mplus_m0*Omega*(F_x_mplus_m0')*F_x_mplus_m0;
%     Omega1 = Omega0 - Omega0*(Fm0')*inv(Fm0*Omega0*Fm0')*Fm0*Omega0;      
%     Omega2 = Omega0 - Omega0*(Fx_m0')*inv(Fx_m0*Omega0*Fx_m0')*Fx_m0*Omega0;
%     Omega3 = Omega - Omega*(Fx')*inv(Fx*Omega*Fx')*Fx*Omega;
%     Omega_chi = Omega1 - Omega2 + Omega3;                                  % �a���X�V�ɂ�������s��
    Omega_chi = Omega - Omega0*(Fm0')*inv(Fm0*Omega0*Fm0')*Fm0*Omega0 +...
                Omega0*(Fx_m0')*inv(Fx_m0*Omega0*Fx_m0')*Fx_m0*Omega0 -...
                Omega*(Fx')*inv(Fx*Omega*Fx')*Fx*Omega;                    % �a���X�V�ɂ�������s��
    Xi_chi = Xi + (Omega_chi - Omega)*xhat;                                % �a���X�V�ɂ�������x�N�g��
%     disp('�a������')
    Xi = Xi_chi;                                                           % �a���X�V�ɂ�������x�N�g���̖߂�l�ϐ���`
    Omega = Omega_chi;                                                     % �a���X�V�ɂ�������s��̖߂�l�ϐ���`
    Omegaxx = Omega(1:3,1:3);
    Omegaxm = Omega(1:3,4:length(Omega));
    Xhatlm = xhat(4:length(xhat));
    Xix = Xi(1:3);
    Xhats = inv(Omegaxx)*(Xix-(Omegaxm*Xhatlm));
    xhat = [Xhats;Xhatlm];                                                 % ��ԃx�N�g���̓��o
%     xhat(3) = PI2PI(xhat(3));                                              % �p�x�␳
    end
%     disp('�ϑ�����')
    end

    
function x = f(x, u)
global dt;
global LMSize;

Fx = horzcat(eye(3),zeros(3,LMSize*NumLM(x))); 
B = [dt*cos(x(3))  0
     dt*sin(x(3))  0
          0       dt];
Delta = B*u;
% Delta = [-(u(1)/u(2))*sin(x(3)) + (u(1)/u(2))*sin(x(3)+(u(2)*dt));
%           (u(1)/u(2))*cos(x(3)) - (u(1)/u(2))*cos(x(3)+(u(2)*dt));
%                                    x(3)*dt                       ];
x = x + Fx'*Delta;

function [G, Jf, Fx, Delta] = jacobF(x, u)
global dt;
global LMSize;
B = [dt*cos(x(3))  0
     dt*sin(x(3))  0
          0       dt];
Delta = B*u;
Fx = horzcat(eye(3),zeros(3,LMSize*NumLM(x))); 
JF = [0 0 -dt*u(1)*sin(x(3));
      0 0  dt*u(1)*cos(x(3));
      0 0         0         ];
% JF = [0 0 (u(1)/u(2))*cos(x(3)) - (u(1)/u(2))*cos(x(3)+(u(2)*dt));
%       0 0 (u(1)/u(2))*sin(x(3)) - (u(1)/u(2))*sin(x(3)+(u(2)*dt));
%       0 0                        0                                ];
Jf = eye(3) - JF;
G = eye(length(x)) + Fx'*JF*Fx;

function zl=CalcLMPosiFromZ(x,z)
zl=x(1:2)+[z(1)*cos(PI2PI(x(3)+z(2)));z(1)*sin(PI2PI(x(3)+z(2)))];                       %�ϑ��l����LM�̈ʒu���v�Z����֐�

function H=jacobH(q,delta,x,I)
%�ϑ����f���̃��R�r�s����v�Z����֐�
sq=sqrt(q);
G=[-sq*delta(1) -sq*delta(2) 0 sq*delta(1) sq*delta(2);
    delta(2)    -delta(1)   -1 -delta(2)    delta(1)];
G=G/q;
F=[eye(3) zeros(3,2*NumLM(x));
   zeros(2,3) zeros(2,2*(I-1)) eye(2) zeros(2,2*NumLM(x)-2*I)];
H=G*F;

function n = NumLM(x)
n = (length(x)-3)/2;

function angle=PI2PI(angle)
%���{�b�g�̊p�x��-pi~pi�͈̔͂ɕ␳����֐�
angle = mod(angle, 2*pi);

i = find(angle>pi);
angle(i) = angle(i) - 2*pi;

i = find(angle<=-pi);
angle(i) = angle(i) + 2*pi;