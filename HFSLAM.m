function [xEst,PEst,zl] = HFSLAM(xEst,PEst,z,u)
global LMSize dt R Q alpha initP gamma delta

if isempty(z)
    F = horzcat(eye(3),zeros(3,LMSize*GetnLM(xEst))); 
    B = [dt*cos(xEst(3))  0
         dt*sin(xEst(3))  0
              0       dt]; 
    xEst = xEst + F'*B*u;
    xEst(3) = PI2PI(xEst(3));                                              
    disp('�ϑ��Ȃ�')
else
    disp('�ϑ�����')
    % ------ EKF SLAM --------
    % Predict
    xEst = f(xEst, u);
    [G,Fx]=jacobF(xEst, u);
    PEst= G'*PEst*G + Fx'*R*Fx;
    % Update
    for iz=1:length(z(:,1))%���ꂼ��̊ϑ��l�ɑ΂��Ċϑ��l�������h�}�[�N�Ƃ��Ēǉ�
        zl=CalcLMPosiFromZ(xEst,z(iz,:));%�ϑ��l���̂��̂���LM�̈ʒu���v�Z
        xAug=[xEst;zl];%��ԃx�N�g���Ƌ����U�s��̒ǉ�
        PAug=[PEst zeros(length(xEst),LMSize);
              zeros(LMSize,length(xEst)) initP];
        
        mdist=[];%�}�n���m�r�X�����̃��X�g
        for il=1:GetnLM(xAug) %���ꂼ��̃����h�}�[�N�ɂ���
            if il==GetnLM(xAug)
                mdist=[mdist alpha];%�V�����ǉ������_�̋����̓p�����[�^�l���g��
            else
                lm=xAug(4+2*(il-1):5+2*(il-1));
                [y,S,~]=CalcInnovation(lm,xAug,PAug,z(iz,1:2),il);
                mdist=[mdist y'/S*y];%�}�n���m�r�X�����̌v�Z
            end
        end
        
        %�}�n���m�r�X�������ł��߂����̂ɑΉ��t��
        [~,I]=min(mdist);
      
        %��ԋ��������������̂��ǉ��������̂Ȃ�΁A���̊ϑ��l�������h�}�[�N�Ƃ��č̗p
        if I==GetnLM(xAug)
            %disp('New LM')
            xEst=xAug;
            PEst=PAug;
        end
        
        lm=xEst(4+2*(I-1):5+2*(I-1));%�Ή��t����ꂽ�����h�}�[�N�f�[�^�̎擾
        %�C�m�x�[�V�����̌v�Z
        [y,S,H]=CalcInnovation(lm,xEst,PEst,z(iz,1:2),I);
        PS = (1 + delta)*(eye(length(xEst)) + H'/Q*H*PEst - gamma^(-2)*eye(length(xEst))*PEst);
        K = PEst/PS*H'/Q;
        xEst = xEst + K*y;
        PEst = (eye(size(xEst,1)) - K*H)*PEst;   
    end
end
    xEst(3)=PI2PI(xEst(3));%�p�x�␳
end

function [y,S,H]=CalcInnovation(lm,xEst,PEst,z,LMId)
%�Ή��t�����ʂ���C�m�x�[�V�������v�Z����֐�
global Q;
delta=lm-xEst(1:2);
q=delta'*delta;
zangle=atan2(delta(2),delta(1))-xEst(3);
zp=[sqrt(q) PI2PI(zangle)];%�ϑ��l�̗\��
y=(z-zp)';
H=jacobH(q,delta,xEst,LMId);
S=H*PEst*H'+Q;
end

function n=GetnLM(xEst)
%�����h�}�[�N�̐����v�Z����֐�
n=(length(xEst)-3)/2;
end

function zl=CalcLMPosiFromZ(x,z)
%�ϑ��l����LM�̈ʒu���v�Z����֐�
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
x(3)=PI2PI(x(3));%�p�x�␳
end

function [G,Fx]=jacobF(x, u)
% �^�����f���̃��R�r�s��̌v�Z�֐�
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
%�ϑ����f���̃��R�r�s����v�Z����֐�
sq=sqrt(q);
G=[-sq*delta(1) -sq*delta(2) 0 sq*delta(1) sq*delta(2);
    delta(2)    -delta(1)   -1 -delta(2)    delta(1)];
G=G/q;
F=[eye(3) zeros(3,2*GetnLM(x));
   zeros(2,3) zeros(2,2*(i-1)) eye(2) zeros(2,2*GetnLM(x)-2*i)];
H=G*F;
end

function angle=PI2PI(angle)
%���{�b�g�̊p�x��-pi~pi�͈̔͂ɕ␳����֐�
angle = mod(angle, 2*pi);

i = find(angle>pi);
angle(i) = angle(i) - 2*pi;

i = find(angle<-pi);
angle(i) = angle(i) + 2*pi;
end