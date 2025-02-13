% 21-08-31 17:34
function y=Gamman(n,x)
%     y=besseli(n,x).*exp(-x);
    % y=besseli(n,x,1);
    y=besseli(n,x,1).* exp(-x + abs(real(x))); % 25-02-07 16:20
    % y=GamnLpole(n,x);
end