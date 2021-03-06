function [xEst,PEst,jx] = EKF_SLAM(u,z,xEst,PEst,initP,alpha,jx)
    % ------ EKF-SLAM --------
    global R;
    
    if isempty(z)                                                          % 観測値が何もない場合は状態のみ更新する
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
        zl = CalcLMPosiFromZ(xEst,z(ilm,:));                           %観測値そのものからLMの位置を計算
        xAug = [xEst;zl];
        PAug = [PEst zeros(length(xEst),LMSize);
                zeros(LMSize,length(xEst)) initP];
        Jx = [jx; z(ilm,3)]; 
        mdist = [];                                                    %マハラノビス距離のリスト
        for il = 1 : length(Jx)                                        % 既知の特徴量による判別
          if il == length(Jx)
             mdist = [mdist alpha];
          else
             lm = xAug(4+2*(il-1):5+2*(il-1));
             [y,S,~] = CalcInnovation(lm,xAug,PAug,z(ilm,1:2),il);
             MD = y' / S * y;                                          %マハラノビス距離
             mdist = [mdist MD];
          end
        end                                                            % 特徴量の判別終了                                                             
        [~,num] = min(mdist);                                          % ノルムの最小値と最小値の行番号の導出
        if num == length(Jx)                                           % 初検知したランドマークの場合
             jx   = Jx;                                                % 特徴量を反映
             xEst = xAug;                                              % 状態量を反映
             PEst = PAug;                                              % 共分散行列を反映
        end  
        lm = xEst(4+2*(num-1):5+2*(num-1));                            %対応付けられたランドマークデータの取得

        %innovation and update
        [y,S,H] = CalcInnovation(lm,xEst,PEst,z(ilm,1:2),num);
        K = PEst*H' / S;
        xEst = xEst + K*y;
        PEst = PEst - K * S * K';
    end
    %     xEst(3) = PI2PI(xEst(3));                                              %角度補正
end

function [y,S,H]=CalcInnovation(lm,xEst,PEst,z,LMId)
    %対応付け結果からイノベーションを計算する関数
    global Q;
    delta = lm - xEst(1:2);                                                %dx,dyの導出(絶対座標系)
    q = delta' * delta;                                                    %dx^2,dy^2の導出
    zangle = atan2(delta(2),delta(1)) - xEst(3);
    zp = [sqrt(q) PI2PI(zangle)];                                          %観測値の予測
    y = (z - zp)';
    H = jacobH(q,delta,xEst,LMId);
    S = H * PEst * H' + Q;
    S = (S + S') * 0.5; % make symmetric
end

function n = NumLM(xEst)
    n = (length(xEst) - 3) / 2;
end

function zl=CalcLMPosiFromZ(x,z)
    zl = x(1:2)+[z(1)*cos(x(3)+z(2));z(1)*sin(x(3)+z(2))];                 %観測値からLMの位置を計算する関数
end

function x = f(x, u)
    global dt;
    global LMSize;

    F = horzcat(eye(3),zeros(3,LMSize*NumLM(x))); 
    B = [dt*cos(x(3))  0
         dt*sin(x(3))  0
              0       dt]; 
    x = x + F'*B*u;
    x(3) = PI2PI(x(3));%角度補正
end

function H=jacobH(q,delta,x,i)
    %観測モデルのヤコビ行列を計算する関数
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
    %ロボットの角度を-pi~piの範囲に補正する関数
    angle = mod(angle, 2*pi);

    i = find(angle>pi);
    angle(i) = angle(i) - 2*pi;

    i = find(angle<=-pi);
    angle(i) = angle(i) + 2*pi;
end