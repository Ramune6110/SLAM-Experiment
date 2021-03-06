function [xEst,Omega,Xi,jx] = SEIF_SLAM_MK(u,z,xEst,Omega_chi,Xi_chi,initOmega,initXi,alpha,jx)
    % ------ SEIF-SLAM --------
    if isempty(z)                                                          % 観測値が何もない場合は状態のみ更新する
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
    
        if NumLM(xEst) < 4                                                 % 疎化更新を行わない場合
            %--------------------------------------------------------
            % Update state estimate
            %--------------------------------------------------------
            [Xi, Omega, xEst] = SEIF_update_state_estimate(Xi, Omega, xEst);
            
        else                                                               % 疎化更新を行う場合
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
    
    [JF,Fx,Delta] = jacobF(xEst, u);                                       % ヤコビ行列の計算(状態方程式)
    Ps        = Fx' * (inv(eye(3) + JF) - eye(3)) * Fx;                    % ヤコビ行列の逆行列を消去するための変数 G^-1 = I + Ps
    lambda    = Ps' * Omega_chi + Omega_chi * Ps + Ps' * Omega_chi * Ps;   % P.358(確率ロボティクス)
    Ph        = Omega_chi + lambda;                                        % 情報行列,及び,EKFの共分散の予測更新の第1項目の変数
    kappa     = Ph * Fx' * inv(Q1 + Fx * Ph * Fx') * Fx * Ph;              % 情報行列,及び,EKFの共分散の予測更新の第2項目の変数
    Omega_bar = Ph - kappa;                                                % 予測更新における情報行列
    Xi_bar    = Xi_chi + (lambda - kappa) * xEst + Omega_bar * Fx' * Delta;% 予測更新における情報ベクトル
    xEst_bar  = xEst + Fx' * Delta;                                        % 事前推定値
end

%--------------------------------------------------------
% Measurement Update
%--------------------------------------------------------
function [Xi, Omega, xEst, numM, zM, jx] = SEIF_measurement_update(Xi, Omega, xEst, z, numM, zM, initOmega, initXi, alpha, jx)
    global R1;
    global LMSize
    
    for ilm = 1: length(z(:,1))                                            % 検知したランドマークすべて計算
        zl       = CalcLMPosiFromZ(xEst,z(ilm,:));                         % 観測値そのものからLMの位置を計算
        Jx       = [jx; z(ilm,3)];                                         % 特徴量追加
        xEstbar  = [xEst; zl];                                             % ランドマーク追加
        Xibar    = [Xi; initXi];                                           % 情報ベクトル追加(ランドマーク追加にあわせて)
        Omegabar = [Omega zeros(length(Omega),LMSize);
                    zeros(LMSize,length(Omega)) initOmega];                % 情報行列追加(ランドマーク追加にあわせて)              
        jM = [];                                                           % 特徴量格納

        for il = 1 : length(Jx)                                            % 既知の特徴量による判別
          if il == length(Jx)
             jM = [jM alpha];
          else
             lm = xEstbar(4+2*(il-1):5+2*(il-1));
             [y,S,~] = CalcInnovation(lm,xEstbar,inv(Omegabar),z(ilm,1:2),il);
             MD = y' / S * y;                                              %マハラノビス距離
             jM = [jM MD];
          end
        end                                                                % 特徴量の判別終了                                                             
        [~,num] = min(jM);                                                 % ノルムの最小値と最小値の行番号の導出
                                                                           % 追加したランドマークにおける新たな定義
        numM = [numM; num];
        zM   = [zM; z(ilm,3)];
        
        if num == length(Jx)                                               % 初検知したランドマークの場合
            jx    = Jx;                                                    % 特徴量を反映
            xEst  = xEstbar;                                               % 状態量を反映
            Xi    = Xibar;                                                 % 情報ベクトルを反映
            Omega = Omegabar;                                              % 情報行列を反映
        end                                                                % 定義の終了
        
        lm    = xEst(4+2*(num-1):5+2*(num-1));                             % 対応付けられたランドマークデータの取得
        delta = lm - xEst(1:2);                                            % dx,dyの導出(絶対座標系)
        q     = delta' * delta;                                            % dx^2,dy^2の導出
        zhat  = [sqrt(q) PI2PI(atan2(delta(2),delta(1))-xEst(3))];         % 観測値の予測
        y     = (z(ilm,1:2) - zhat)';                                      % イノベーション
        H     = jacobH(q,delta,xEst,num);                                  % ヤコビ行列計算(観測方程式)
        Xi    = Xi + H' * R1 * (y + H * xEst);                             % 計測更新(情報ベクトル)
        Omega = Omega + H' * R1 * H;                                       % 計測更新(情報行列)
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
                Omega*(Fx')*inv(Fx*Omega*Fx')*Fx*Omega;                    % 疎化更新における情報行列
    Xi_chi = Xi + (Omega_chi - Omega)*xEst;                                % 疎化更新における情報ベクトル
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
    xEst    = [xEsts;xEstlm];                                                 % 状態ベクトルの導出
%     xEst(3) = PI2PI(xEst(3));                                              % 角度補正
end

function [y,S,H]=CalcInnovation(lm,xEst,P,z,LMId)
    %対応付け結果からイノベーションを計算する関数
    global Q;
    delta=lm-xEst(1:2);                                                        %dx,dyの導出(絶対座標系)
    q=delta'*delta;                                                            %dx^2,dy^2の導出
    zangle=atan2(delta(2),delta(1))-xEst(3);
    zp=[sqrt(q) PI2PI(zangle)];                                                %観測値の予測
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
    zl=x(1:2)+[z(1)*cos(PI2PI(x(3)+z(2)));z(1)*sin(PI2PI(x(3)+z(2)))];                       %観測値からLMの位置を計算する関数
end

function H=jacobH(q,delta,x,I)
    %観測モデルのヤコビ行列を計算する関数
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
    %ロボットの角度を-pi~piの範囲に補正する関数
    angle = mod(angle, 2*pi);

    i = find(angle>pi);
    angle(i) = angle(i) - 2*pi;

    i = find(angle<=-pi);
    angle(i) = angle(i) + 2*pi;
end