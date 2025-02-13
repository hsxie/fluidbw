function RC = RASY(Y, AL)
    % 参数类型定义
    d2p = 'double'; % 双精度类型
    RC = zeros(2, 2, 'like', complex(1.0, 0.0)); % 初始化复数矩阵
    
    % 变量计算
    Y2 = Y^2;
    COT = cos(pi*Y) / sin(pi*Y);
    C = 7.2e35; % 替代原代码中的 1.E99
    
    % 初始化变量
    PN = -Y / AL;
    PYN = PN;
    A = 1 / (AL * sqrt(2 * pi * AL));
    QN = pi * Y2 * COT * A;
    QYN = QN * (2 - Y * pi * COT) - Y * pi^2 * Y2 * A;
    
    P = PN + QN;
    PY = PYN + QYN;
    PP = -PN - 1.5 * QN;
    PPY = -PYN - 1.5 * QYN;
    AY = abs(Y) + 2;
    
    T = 0; % 初始化T
    
    % 迭代循环
    for N = 1:100
        M = N - 1;
        PYN = (PYN * (M^2 - Y2) - 2 * Y2 * PN) / ((2*M + 1) * AL);
        PN = PN * (M^2 - Y2) / ((2*M + 1) * AL);
        QYN = (QYN * ((M + 0.5)^2 - Y2) - 2 * Y2 * QN) / (2 * N * AL);
        QN = QN * ((M + 0.5)^2 - Y2) / (2 * N * AL);
        
        if M >= AY
            C = N * (abs(PN) + abs(QN));
            if C <= 1e-7 * abs(PP) || C >= T
                break;
            end
        end
        
        P = P + PN + QN;
        PY = PY + PYN + QYN;
        PP = PP - (N + 1) * PN - (N + 1.5) * QN;
        PPY = PPY - (N + 1) * PYN - (N + 1.5) * QYN;
        T = C;
    end
    
    % 输出结果
    RC(1,1) = P + PN + QN;
    RC(2,1) = PY + PYN + QYN;
    RC(1,2) = PP + P;
    RC(2,2) = PPY + PY;
end