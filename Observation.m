function [z,x,u] = Observation(x, u, LM, MAX_RANGE,Noise11,Noise22)
global Qsigma Rsigma;
x=f(x, u);% Ground Truth
u=u+sqrt(Qsigma)*(Noise11);%add Process Noise
%Simulate Observation
z=[];
for iz=1:length(LM(:,1))
    %LM�̈ʒu�����{�b�g���W�n�ɕϊ�
    yaw=zeros(3,1);
    yaw(3)=-x(3);
    localLM=HomogeneousTransformation2D(LM(iz,:)-x(1:2)',yaw');
    d=norm(localLM);%����
    if d<MAX_RANGE %�ϑ��͈͓�
        noise=sqrt(Rsigma)*(Noise22);
        z=[z;[d+noise(1) PI2PI(atan2(localLM(2),localLM(1))+noise(2)) LM(iz,:)]];
    end
end
end

function n=GetnLM(xEst)
%�����h�}�[�N�̐����v�Z����֐�
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
x(3)=PI2PI(x(3));%�p�x�␳
end

function out = HomogeneousTransformation2D(in, base, mode)
%��]�s��
Rot=[cos(base(3)) sin(base(3)); -sin(base(3)) cos(base(3))];

%�_��������base���W��z��Ɋi�[
Nin=size(in);
baseMat=repmat(base(1:2),Nin(1),1);

% x-y�ȊO�̃f�[�^�������Ă����ꍇ�D��񑼂̕ϐ��ɒu���Ă����D
if Nin(2)>=3
    inxy=in(:,1:2);
    inOther=in(:,3:end);
    in=inxy;
end

%�����ϊ�
if nargin==2 || mode==0 %��]�����i
    out=baseMat+in*Rot;
else %���i����]
    out=(baseMat+in)*Rot;
end
    
%��菜�����l����������D
if Nin(2)>=3
    out=[out inOther];
end
end

function angle=PI2PI(angle)
%���{�b�g�̊p�x��-pi~pi�͈̔͂ɕ␳����֐�
angle = mod(angle, 2*pi);

i = find(angle>pi);
angle(i) = angle(i) - 2*pi;

i = find(angle<-pi);
angle(i) = angle(i) + 2*pi;
end