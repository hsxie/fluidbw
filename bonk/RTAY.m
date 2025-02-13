function RC = RTAY(Y, AL)
    % 初始化复数矩阵
    RC = zeros(2, 2, 'like', complex(1.0, 0.0));
    Y2 = Y^2;
    
    % 初始值计算
    PN = Y / (Y2 - 1);
    PYN = -Y*(Y2 + 1) / (Y2 - 1)^2;
    RC(1,1:2) = PN;
    RC(2,1:2) = PYN;
    
    % 泰勒级数迭代
    for I = 2:100
        COT = (2*I - 1) / (Y2 - I^2) * AL;
        PYN = COT * (PYN - 2*Y2/(Y2 - I^2) * PN);
        PN = COT * PN;
        RC(1,1) = RC(1,1) + PN;
        RC(2,1) = RC(2,1) + PYN;
        RC(1,2) = RC(1,2) + I*PN;
        RC(2,2) = RC(2,2) + I*PYN;
        
        T = abs(PN) * 1e8;
        if T < abs(RC(1,1))
            break;
        end
    end
end