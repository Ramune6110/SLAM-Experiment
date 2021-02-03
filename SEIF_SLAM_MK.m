function [xEst,Omega,Xi,jx] = SEIF_SLAM_MK(u,z,xEst,Omega_chi,Xi_chi,initOmega,initXi,alpha,jx)
    % ------ SEIF-SLAM --------
    if isempty(z)                                                          % �ϑ��l�������Ȃ��ꍇ�͏�Ԃ̂ݍX�V����
         [xEst, Omega, Xi] = SEIF_motion_update(xEst, Omega_chi, Xi_chi, u);
    else
        %--------------------------------------------------------
        % Motion Update
        %--------------------------------------------------------
        [xEst_bar, Omega_bar, Xi_bar] = SEIF_motion_update(xEst, Omega_chi, Xi_chi, u);

        %--------------------------------------------------------
        % Measurement Update
        %--------------------------------------------------------
        numM = [];
        zM   = [];
        
        [Xi, Omega, xEst, numM, zM, jx] = SEIF_measurement_update(Xi_bar, Omega_bar, xEst_bar, z, numM, zM, initOmega, initXi, alpha, jx);
    
        if NumLM(xEst) < 4                                                 % �a���X�V���s��Ȃ��ꍇ
            %--------------------------------------------------------
            % Update state estimate
            %--------------------------------------------------------
            [Xi, Omega, xEst] = SEIF_update_state_estimate(Xi, Omega, xEst);
            
        else                                                               % �a���X�V���s���ꍇ
            %--------------------------------------------------------
            % Sparsification
            %--------------------------------------------------------
            [Xi_chi, Omega_chi, jx] = SEIF_sparsification(Xi, Omega, xEst, numM, zM, jx);
            
            %--------------------------------------------------------
            % Update state estimate
            %--------------------------------------------------------
            [Xi, Omega, xEst] = SEIF_update_state_estimate(Xi_chi, Omega_chi, xEst);
       
        end
    end
end

%--------------------------------------------------------
% Motion Update
%--------------------------------------------------------
function [xEst_bar, Omega_bar, Xi_bar] = SEIF_motion_update(xEst, Omega_chi, Xi_chi, u)
    global Q1;
    
    [JF,Fx,Delta] = jacobF(xEst, u);                                       % ���R�r�s��̌v�Z(��ԕ�����)
    Ps        = Fx' * (inv(eye(3) + JF) - eye(3)) * Fx;                    % ���R�r�s��̋t�s����������邽�߂̕ϐ� G^-1 = I + Ps
    lambda    = Ps' * Omega_chi + Omega_chi * Ps + Ps' * Omega_chi * Ps;   % P.358(�m�����{�e�B�N�X)
    Ph        = Omega_chi + lambda;                                        % ���s��,�y��,EKF�̋����U�̗\���X�V�̑�1���ڂ̕ϐ�
    kappa     = Ph * Fx' * inv(Q1 + Fx * Ph * Fx') * Fx * Ph;              % ���s��,�y��,EKF�̋����U�̗\���X�V�̑�2���ڂ̕ϐ�
    Omega_bar = Ph - kappa;                                                % �\���X�V�ɂ�������s��
    Xi_bar    = Xi_chi + (lambda - kappa) * xEst + Omega_bar * Fx' * Delta;% �\���X�V�ɂ�������x�N�g��
    xEst_bar  = xEst + Fx' * Delta;                                        % ���O����l
end

%--------------------------------------------------------
% Measurement Update
%--------------------------------------------------------
function [Xi, Omega, xEst, numM, zM, jx] = SEIF_measurement_update(Xi, Omega, xEst, z, numM, zM, initOmega, initXi, alpha, jx)
    global R1;
    global LMSize
    
    for ilm = 1: length(z(:,1))                                            % ���m���������h�}�[�N���ׂČv�Z
        zl       = CalcLMPosiFromZ(xEst,z(ilm,:));                         % �ϑ��l���̂��̂���LM�̈ʒu���v�Z
        Jx       = [jx; z(ilm,3)];                                         % �����ʒǉ�
        xEstbar  = [xEst; zl];                                             % �����h�}�[�N�ǉ�
        Xibar    = [Xi; initXi];                                           % ���x�N�g���ǉ�(�����h�}�[�N�ǉ��ɂ��킹��)
        Omegabar = [Omega zeros(length(Omega),LMSize);
                    zeros(LMSize,length(Omega)) initOmega];                % ���s��ǉ�(�����h�}�[�N�ǉ��ɂ��킹��)              
        jM = [];                                                           % �����ʊi�[

        for il = 1 : length(Jx)                                            % ���m�̓����ʂɂ�锻��
          if il == length(Jx)
             jM = [jM alpha];
          else
             lm = xEstbar(4+2*(il-1):5+2*(il-1));
             [y,S,~] = CalcInnovation(lm,xEstbar,inv(Omegabar),z(ilm,1:2),il);
             MD = y' / S * y;                                              %�}�n���m�r�X����
             jM = [jM MD];
          end
        end                                                                % �����ʂ̔��ʏI��                                                             
        [~,num] = min(jM);                                                 % �m�����̍ŏ��l�ƍŏ��l�̍s�ԍ��̓��o
                                                                           % �ǉ����������h�}�[�N�ɂ�����V���Ȓ�`
        numM = [numM; num];
        zM   = [zM; z(ilm,3)];
        
        if num == length(Jx)                                               % �����m���������h�}�[�N�̏ꍇ
            jx    = Jx;                                                    % �����ʂ𔽉f
            xEst  = xEstbar;                                               % ��ԗʂ𔽉f
            Xi    = Xibar;                                                 % ���x�N�g���𔽉f
            Omega = Omegabar;                                              % ���s��𔽉f
        end                                                                % ��`�̏I��
        
        lm    = xEst(4+2*(num-1):5+2*(num-1));                             % �Ή��t����ꂽ�����h�}�[�N�f�[�^�̎擾
        delta = lm - xEst(1:2);                                            % dx,dy�̓��o(��΍��W�n)
        q     = delta' * delta;                                            % dx^2,dy^2�̓��o
        zhat  = [sqrt(q) PI2PI(atan2(delta(2),delta(1))-xEst(3))];         % �ϑ��l�̗\��
        y     = (z(ilm,1:2) - zhat)';                                      % �C�m�x�[�V����
        H     = jacobH(q,delta,xEst,num);                                  % ���R�r�s��v�Z(�ϑ�������)
        Xi    = Xi + H' * R1 * (y + H * xEst);                             % �v���X�V(���x�N�g��)
        Omega = Omega + H' * R1 * H;                                       % �v���X�V(���s��)
    end
end

%--------------------------------------------------------
% Sparsification
%--------------------------------------------------------
function [Xi_chi, Omega_chi, jx] = SEIF_sparsification(Xi, Omega, xEst, numM, zM, jx)
    FF = eye(length(xEst(:,1)));
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

    Omega_chi = Omega - Omega0*(Fm0')*inv(Fm0*Omega0*Fm0')*Fm0*Omega0 +...
                Omega0*(Fx_m0')*inv(Fx_m0*Omega0*Fx_m0')*Fx_m0*Omega0 -...
                Omega*(Fx')*inv(Fx*Omega*Fx')*Fx*Omega;                    % �a���X�V�ɂ�������s��
    Xi_chi = Xi + (Omega_chi - Omega)*xEst;                                % �a���X�V�ɂ�������x�N�g��
end

%--------------------------------------------------------
% Update state estimate
%--------------------------------------------------------
function [Xi, Omega, xEst] = SEIF_update_state_estimate(Xi, Omega, xEst)
    Omegaxx = Omega(1:3,1:3);
    Omegaxm = Omega(1:3,4:length(Omega));
    xEstlm  = xEst(4:length(xEst));
    Xix     = Xi(1:3);
    xEsts   = inv(Omegaxx)*(Xix-(Omegaxm*xEstlm));
    xEst    = [xEsts;xEstlm];                                                 % ��ԃx�N�g���̓��o
%     xEst(3) = PI2PI(xEst(3));                                              % �p�x�␳
end

function [y,S,H]=CalcInnovation(lm,xEst,P,z,LMId)
    %�Ή��t�����ʂ���C�m�x�[�V�������v�Z����֐�
    global Q;
    delta=lm-xEst(1:2);                                                        %dx,dy�̓��o(��΍��W�n)
    q=delta'*delta;                                                            %dx^2,dy^2�̓��o
    zangle=atan2(delta(2),delta(1))-xEst(3);
    zp=[sqrt(q) PI2PI(zangle)];                                                %�ϑ��l�̗\��
    y=(z-zp)';
    H=jacobH(q,delta,xEst,LMId);
    S=H*P*H'+Q;
    S = (S+S')*0.5; % make symmetric
end

function [JF, Fx, Delta] = jacobF(x, u)
    global dt;
    global LMSize;
    
    B = [dt*cos(x(3))  0
         dt*sin(x(3))  0
              0       dt];
    
    Fx = horzcat(eye(3),zeros(3,LMSize*NumLM(x))); 
    Delta = B*u;
    JF = [0 0 -dt*u(1)*sin(x(3));
          0 0  dt*u(1)*cos(x(3));
          0 0         0         ];
end

function zl=CalcLMPosiFromZ(x,z)
    zl=x(1:2)+[z(1)*cos(PI2PI(x(3)+z(2)));z(1)*sin(PI2PI(x(3)+z(2)))];                       %�ϑ��l����LM�̈ʒu���v�Z����֐�
end

function H=jacobH(q,delta,x,I)
    %�ϑ����f���̃��R�r�s����v�Z����֐�
    sq=sqrt(q);
    G=[-sq*delta(1) -sq*delta(2) 0 sq*delta(1) sq*delta(2);
        delta(2)    -delta(1)   -1 -delta(2)    delta(1)];
    G=G/q;
    F=[eye(3) zeros(3,2*NumLM(x));
       zeros(2,3) zeros(2,2*(I-1)) eye(2) zeros(2,2*NumLM(x)-2*I)];
    H=G*F;
end

function n = NumLM(x)
    n = (length(x)-3)/2;
end

function angle=PI2PI(angle)
    %���{�b�g�̊p�x��-pi~pi�͈̔͂ɕ␳����֐�
    angle = mod(angle, 2*pi);

    i = find(angle>pi);
    angle(i) = angle(i) - 2*pi;

    i = find(angle<=-pi);
    angle(i) = angle(i) + 2*pi;
end