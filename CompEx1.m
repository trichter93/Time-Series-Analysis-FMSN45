%% Computer Exercise 1 - Time Series analysis

%% Task 1

A1 = [ 1 -1.79 0.84 ] ;
C1 = [ 1 -0.18 -0.11 ] ;

A2 = [ 1 -1.79 ] ;
C2 = [ 1 -0.18 -0.11 ] ;
%%
ARMA_poly1 = idpoly(A1, [], C1);
ARMA_poly2 = idpoly(A2, [], C2);
%%

