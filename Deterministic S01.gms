$OFFLISTING
$EOLCOM //
$INLINECOM { }
$TITLE Water  Management
* Options:
OPTION LIMROW = 0; OPTION LIMCOL = 0;
OPTION SOLPRINT = OFF; OPTION SYSOUT = OFF;


set i vertices nodes / 1, 2, 3, 4, 5, 6, 7, 8, 9, 10  /,
    j arcs edges / 12, 13, 14, 25, 26, 27, 35, 36, 37, 45, 46, 47,
                    58, 59, 510, 68, 69, 610, 78, 79, 710/,
    S indices of realizations / S01*S02 /
    ;
    
parameter c(j) costs / 12 15.0, 13 13.0, 14 7.0, 25 10.0, 26 2.0, 27 4.0, 35 8.0,
                        36 7.0, 37 4.0, 45 10.0, 46 5.0, 47 4.0, 58 5.0, 59 3.0,
                        510 6.0, 68 7.0, 69 2.0, 610 5.0, 78 4.0, 79 10.0, 710 11.0/;
                                                
table bxi(i,S) scenario values for b(i)
    S01 S02
1   40  40
2   18  18
3   0   0
4   0   0
5   0   0
6   0   0
7   0   0
8  -20 -15
9  -25 -11
10 -10 -17;


table a(i,j) incidence node edge matrix
    12  13  14  25  26  27  35  36  37  45  46  47  58  59  510  68  69  610  78  79  710
1   1   1   1
2  -1           1   1   1
3      -1                   1   1   1
4          -1                           1   1   1
5              -1          -1          -1           1    1    1
6                  -1          -1          -1                     1   1    1
7                      -1          -1          -1                              1   1   1
8                                                  -1            -1           -1
9                                                       -1           -1           -1
10                                                           -1           -1          -1;




variables z;
positive variables x(j);
equations objfunc, constr;
objfunc..     z =E= sum(j, c(j)*x(j));
constr(i)..   sum(j, a(i,j) * x(j) ) =L= bxi(i,"S01");
model water / objfunc, constr /;
solve water minimizing z using LP;
display i, j, c, bxi, a, z.L, x.L;



$ontext
min z
           z =E= 3 * x("12") + 6 * x("13") + 2 * x("23");
constr("1")..   +1 * x("12") + 1 * x("13") + 0 * x("23")  =L= b("1");
constr("2")..   -1 * x("12") + 0 * x("13") + 1 * x("23")  =L= b("2");
                 0 * x("12") - 1 * x("13") - 1 * x("23")  =L= b("3");
$offtext