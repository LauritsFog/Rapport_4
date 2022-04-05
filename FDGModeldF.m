function [dF] = FDGModeldF(t,C,K,CA)
    dF = [K(1)*CA-(K(2)+K(3))*C(1)+K(4)*C(2);
          K(3)*C(1)-(K(4)+K(5))*C(2)+K(6)*C(3);
          K(5)*C(2)-(K(6)+K(7))*C(3)+K(8)*C(4);
          K(7)*C(3)-(K(8)+K(9))*C(4)+K(10)*C(5);
          K(9)*C(4)-K(10)*C(5)];
end