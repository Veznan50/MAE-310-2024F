%% Matlab mesh
%% plate_mesh, Created by Gmsh
%% ASCII
clear msh;
msh.nbNod = 46;
msh.POS = [
1 -1 0;
1 1 0;
-1 1 0;
-1 -1 0;
-0.7 -1 0;
-1 -0.7 0;
-0.7878679656440357 -0.7878679656440357 0;
-0.7057644158966322 -0.9414729033066735 0;
-0.7228361403694678 -0.8851949699938777 0;
-0.7505591163718729 -0.8333289300003772 0;
-0.8333289301203975 -0.7505591162916778 0;
-0.8851949705866085 -0.7228361401239506 0;
-0.9414729038091568 -0.705764415796682 0;
-1 -0.2750000000008892 0;
-1 0.1499999999978956 0;
-1 0.5749999999991542 0;
-0.5000000000013867 1 0;
-2.750244476601438e-12 1 0;
0.499999999998614 1 0;
1 0.5000000000013867 0;
1 2.750244476601438e-12 0;
1 -0.499999999998614 0;
0.5750000000008677 -1 0;
0.1500000000021044 -1 0;
-0.2749999999991541 -1 0;
0.5530330085900215 0.5530330085900215 0;
0.1060660171802601 0.1060660171802601 0;
-0.3409009742318874 -0.3409009742318874 0;
-0.6103682259530723 0.5735588960504502 0;
-0.7207364519047037 0.1471177921006458 0;
-0.8311046778569304 -0.2793233118478963 0;
-0.2212987426482009 0.5692909649691015 0;
-0.4425974852935426 0.138581929938106 0;
-0.6638962279400757 -0.2921271050928416 0;
0.166667767469633 0.5623602209276399 0;
-0.1666644650591853 0.1247204418553394 0;
-0.4999966975897913 -0.3129193372181286 0;
0.5623602209080131 0.1666677675017143 0;
0.5692909649085705 -0.2212987424958964 0;
0.5735588960267415 -0.6103682258253744 0;
0.1247204418162788 -0.1666644649977964 0;
0.1385819298174308 -0.4425974849944357 0;
0.1471177920538114 -0.7207364516520813 0;
-0.3129193372778477 -0.4999966974990862 0;
-0.2921271052761204 -0.663896227494156 0;
-0.2793233119215645 -0.831104677479377 0;
];
msh.MAX = max(msh.POS);
msh.MIN = min(msh.POS);
msh.LINES =[
 5 8 0
 8 9 0
 9 10 0
 10 7 0
 7 11 0
 11 12 0
 12 13 0
 13 6 0
 6 14 0
 14 15 0
 15 16 0
 16 3 0
 3 17 0
 17 18 0
 18 19 0
 19 2 0
 2 20 0
 20 21 0
 21 22 0
 22 1 0
 1 23 0
 23 24 0
 24 25 0
 25 5 0
 2 26 0
 26 27 0
 27 28 0
 28 7 0
];
msh.QUADS =[
 3 17 29 16 0
 16 29 30 15 0
 15 30 31 14 0
 14 31 13 6 0
 17 18 32 29 0
 29 32 33 30 0
 30 33 34 31 0
 31 34 12 13 0
 18 19 35 32 0
 32 35 36 33 0
 33 36 37 34 0
 34 37 11 12 0
 19 2 26 35 0
 35 26 27 36 0
 36 27 28 37 0
 37 28 7 11 0
 2 26 38 20 0
 20 38 39 21 0
 21 39 40 22 0
 22 40 23 1 0
 26 27 41 38 0
 38 41 42 39 0
 39 42 43 40 0
 40 43 24 23 0
 27 28 44 41 0
 41 44 45 42 0
 42 45 46 43 0
 43 46 25 24 0
 28 7 10 44 0
 44 10 9 45 0
 45 9 8 46 0
 46 8 5 25 0
];
msh.PNT =[
 1 0
 2 0
 3 0
 4 0
 5 0
 6 0
 7 0
];
