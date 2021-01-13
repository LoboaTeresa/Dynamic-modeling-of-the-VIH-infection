function [dy]=VIH_eqs(t,y)

%Paramteros

global lambda k d p delta0 c t1 t2 inT1 inT2 beta;
%T es y(1); I=y(2); V=y(3)
fun=@(t)( ( beta ./ ( 1 + k.*exp(-(t-t1)./inT1) ) ) - ( beta ./ ( 1 + k.*exp(-(t-t2)./inT2) ) ) );

dy=[ lambda - k.*y(3).*y(1) - d.*y(1) ;
    ( k.*y(3).*y(1) - delta0.*y(3).*y(2) );
     ( p.*y(2) - c.*y(3) )];
end
