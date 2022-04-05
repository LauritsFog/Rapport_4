function [dF] = FDGModeldF(t,C,K,CA)
    dF = [K(1)*CA-(K(2)+K(3))*C(1)+K(4)*C(2);
          K(3)*C(1)-K(4)*C(2)];
end