function [xEst,PEst,jx] = EKF_SLAM(u,z,xEst,PEst,initP,alpha,jx)
    % ------ EKF-SLAM --------
    global R;
    
    if isempty(z)                                                          % �ϑ��l�������Ȃ��ꍇ�͏�Ԃ̂ݍX�V����
        xEst = f(xEst, u);
    else
        %--------------------------------------------------------
        % Predict
        %--------------------------------------------------------
        [xEst, PEst] = predict(xEst, PEst, u, R);
        %--------------------------------------------------------
        % Update
        %--------------------------------------------------------
        [xEst, PEst, jx] = update(xEst, PEst, initP, z, alpha, jx);
    end
end

%--------------------------------------------------------
% Predict
%--------------------------------------------------------
function [xEst, PEst] = predict(xEst, PEst, u, R)
   xEst = f(xEst, u);
   [G, Fx] = jacobF(xEst, u);
   PEst = G'*PEst*G + Fx'*R*Fx;    
end

%--------------------------------------------------------
% Update
%--------------------------------------------------------
function [xEst, PEst, jx] = update(xEst, PEst, initP, z, alpha, jx)
    global LMSize;    
  
    for ilm = 1: length(z(:,1))
        zl = CalcLMPosiFromZ(xEst,z(ilm,:));                           %�ϑ��l���̂��̂���LM�̈ʒu���v�Z
        xAug = [xEst;zl];
        PAug = [PEst zeros(length(xEst),LMSize);
                zeros(LMSize,length(xEst)) initP];
        Jx = [jx; z(ilm,3)]; 
        mdist = [];                                                    %�}�n���m�r�X�����̃��X�g
        for il = 1 : length(Jx)                                        % ���m�̓����ʂɂ�锻��
          if il == length(Jx)
             mdist = [mdist alpha];
          else
             lm = xAug(4+2*(il-1):5+2*(il-1));
             [y,S,~] = CalcInnovation(lm,xAug,PAug,z(ilm,1:2),il);
             MD = y' / S * y;                                          %�}�n���m�r�X����
             mdist = [mdist MD];
          end
        end                                                            % �����ʂ̔��ʏI��                                                             
        [~,num] = min(mdist);                                          % �m�����̍ŏ��l�ƍŏ��l�̍s�ԍ��̓��o
        if num == length(Jx)                                           % �����m���������h�}�[�N�̏ꍇ
             jx   = Jx;                                                % �����ʂ𔽉f
             xEst = xAug;                                              % ��ԗʂ𔽉f
             PEst = PAug;                                              % �����U�s��𔽉f
        end  
        lm = xEst(4+2*(num-1):5+2*(num-1));                            %�Ή��t����ꂽ�����h�}�[�N�f�[�^�̎擾

        %innovation and update
        [y,S,H] = CalcInnovation(lm,xEst,PEst,z(ilm,1:2),num);
        K = PEst*H' / S;
        xEst = xEst + K*y;
        PEst = PEst - K * S * K';
    end
    %     xEst(3) = PI2PI(xEst(3));                                              %�p�x�␳
end

function [y,S,H]=CalcInnovation(lm,xEst,PEst,z,LMId)
    %�Ή��t�����ʂ���C�m�x�[�V�������v�Z����֐�
    global Q;
    delta = lm - xEst(1:2);                                                %dx,dy�̓��o(��΍��W�n)
    q = delta' * delta;                                                    %dx^2,dy^2�̓��o
    zangle = atan2(delta(2),delta(1)) - xEst(3);
    zp = [sqrt(q) PI2PI(zangle)];                                          %�ϑ��l�̗\��
    y = (z - zp)';
    H = jacobH(q,delta,xEst,LMId);
    S = H * PEst * H' + Q;
    S = (S + S') * 0.5; % make symmetric
end

function n = NumLM(xEst)
    n = (length(xEst) - 3) / 2;
end

function zl=CalcLMPosiFromZ(x,z)
    zl = x(1:2)+[z(1)*cos(x(3)+z(2));z(1)*sin(x(3)+z(2))];                 %�ϑ��l����LM�̈ʒu���v�Z����֐�
end

function x = f(x, u)
    global dt;
    global LMSize;

    F = horzcat(eye(3),zeros(3,LMSize*NumLM(x))); 
    B = [dt*cos(x(3))  0
         dt*sin(x(3))  0
              0       dt]; 
    x = x + F'*B*u;
    x(3) = PI2PI(x(3));%�p�x�␳
end

function H=jacobH(q,delta,x,i)
    %�ϑ����f���̃��R�r�s����v�Z����֐�
    sq = sqrt(q);
    G = [-sq*delta(1) -sq*delta(2) 0 sq*delta(1) sq*delta(2);
         delta(2)    -delta(1)   -1 -delta(2)    delta(1)];
    G = G / q;
    F = [eye(3) zeros(3,2*NumLM(x));
         zeros(2,3) zeros(2,2*(i-1)) eye(2) zeros(2,2*NumLM(x)-2*i)];
    H = G * F;
end

function [G, Fx] = jacobF(x, u)
    global dt;
    global LMSize;

    Fx = horzcat(eye(3),zeros(3,LMSize*NumLM(x))); 
    JF = [0 0 -dt*u(1)*sin(x(3));
          0 0  dt*u(1)*cos(x(3));
          0 0         0         ];
    G = eye(length(x)) + Fx'*JF*Fx;
end

function angle=PI2PI(angle)
    %���{�b�g�̊p�x��-pi~pi�͈̔͂ɕ␳����֐�
    angle = mod(angle, 2*pi);

    i = find(angle>pi);
    angle(i) = angle(i) - 2*pi;

    i = find(angle<=-pi);
    angle(i) = angle(i) + 2*pi;
end