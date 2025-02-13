function RC = RYLA(Y, AL)
    % 根据条件选择计算方法
    RC = zeros(2, 2, 'like', complex(1.0, 0.0));
    AY = abs(Y);
    
    if AL < 4
        RC = RTAY(Y, AL);
        return;
    end
    
    if AY^2 > 75*AL
        RC = RTAY(Y, AL);
        return;
    end
    
    if AY > 40 + AL/3
        RC = RTAY(Y, AL);
        return;
    end
    
    if 3*(AL - 10) > AY && AY^2 < 15*AL
        RC = RASY(Y, AL);
    else
        RC = RINT(Y, AL);
    end
end