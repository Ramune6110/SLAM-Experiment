function u = DynamicWindowApproach(x,xTrue,model,goal,evalParam,ob,R)
x(1) = xTrue(1);
x(2) = xTrue(2);
x(3) = xTrue(3);
%Dynamic Window[vmin,vmax,��min,��max]�̍쐬
Vr=CalcDynamicWindow(x,model);
%�]���֐��̌v�Z
evalDB=Evaluation(x,Vr,goal,ob,R,model,evalParam);

if isempty(evalDB)
    disp('no path to goal!!');
    u=[0;0];return;
end

%�e�]���֐��̐��K��
evalDB=NormalizeEval(evalDB);

%�ŏI�]���l�̌v�Z
feval=[];
for id=1:length(evalDB(:,1))
    feval=[feval;evalParam(1:3)*evalDB(id,3:5)'];
end
evalDB=[evalDB feval];

[~,ind]=max(feval);%�ł��]���l���傫�����͒l�̃C���f�b�N�X���v�Z
u=evalDB(ind,1:2)';%�]���l���������͒l��Ԃ�

function evalDB=Evaluation(x,Vr,goal,ob,R,model,evalParam)
%�e�p�X�ɑ΂��ĕ]���l���v�Z����֐�
evalDB=[];

for vt=Vr(1):model(5):Vr(2)
    for ot=Vr(3):model(6):Vr(4)
        %�O�Ղ̐���
        xt=GenerateTrajectory(x,vt,ot,evalParam(4),model);
        %�e�]���֐��̌v�Z
        heading=CalcHeadingEval(xt,goal);
        dist=CalcDistEval(xt,ob,R);
        vel=abs(vt);
        evalDB=[evalDB;[vt ot heading dist vel]]; 
    end
end

function EvalDB=NormalizeEval(EvalDB)
%�]���l�𐳋K������֐�
if sum(EvalDB(:,3))~=0
    EvalDB(:,3)=EvalDB(:,3)/sum(EvalDB(:,3));
end
if sum(EvalDB(:,4))~=0
    EvalDB(:,4)=EvalDB(:,4)/sum(EvalDB(:,4));
end
if sum(EvalDB(:,5))~=0
    EvalDB(:,5)=EvalDB(:,5)/sum(EvalDB(:,5));
end

function x=GenerateTrajectory(x,vt,ot,evaldt,~)
%�O�Ճf�[�^���쐬����֐�
global dt;
time=0;
u=[vt;ot];%���͒l
while time<=evaldt
    time=time+dt;%�V�~�����[�V�������Ԃ̍X�V
    x=f(x,u);%�^�����f���ɂ�鐄��
end

function dist=CalcDistEval(x,ob,R)
%��Q���Ƃ̋����]���l���v�Z����֐�

dist=2;
for io=1:length(ob(:,1))
    disttmp=norm(ob(io,:)-x(1:2)')-R;%�p�X�̈ʒu�Ə�Q���Ƃ̃m�����덷���v�Z
    if dist>disttmp%�ŏ��l��������
        dist=disttmp;
    end
end

function heading=CalcHeadingEval(x,goal)
%heading�̕]���֐����v�Z����֐�

theta=toDegree(x(3));%���{�b�g�̕���
goalTheta=toDegree(atan2(goal(2)-x(2),goal(1)-x(1)));%�S�[���̕���

if goalTheta>theta
    targetTheta=goalTheta-theta;%�S�[���܂ł̕��ʍ���[deg]
else
    targetTheta=theta-goalTheta;%�S�[���܂ł̕��ʍ���[deg]
end

heading=180-targetTheta;

function Vr=CalcDynamicWindow(x,model)
%���f���ƌ��݂̏�Ԃ���DyamicWindow���v�Z
global dt;
%�ԗ����f���ɂ��Window
Vs=[0 model(1) -model(2) model(2)];

%�^�����f���ɂ��Window
Vd=[x(4)-model(3)*dt x(4)+model(3)*dt x(5)-model(4)*dt x(5)+model(4)*dt];

%�ŏI�I��Dynamic Window�̌v�Z
Vtmp=[Vs;Vd];
%[vmin,vmax,��min,��max]
Vr=[max(Vtmp(:,1)) min(Vtmp(:,2)) max(Vtmp(:,3)) min(Vtmp(:,4))];

function x = f(x, u)
% Motion Model
global dt;
 
F = [1 0 0 0 0
     0 1 0 0 0
     0 0 1 0 0
     0 0 0 0 0
     0 0 0 0 0];
 
B = [dt*cos(x(3)) 0
    dt*sin(x(3)) 0
    0 dt
    1 0
    0 1];

x= F*x+B*u;

function degree = toDegree(radian)
% radian to degree
degree = radian/pi*180;
