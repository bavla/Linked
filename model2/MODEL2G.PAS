{STRAN-PC/MODEL 2 - generalized equivalences by local optimization }

PROGRAM Model2g;
{M 25000, 0, 0 }

USES Crt, Dos, Stenv;

{
   PROGRAM PACKAGE      S T R A N  / PC
   PROGRAM UNIT         GENERALIZED EQUIVALENCES
   PROGRAMMING LANGUAGE     pascal
   PROGRAMMED BY        Vladimir BATAGELJ

   PROGRAM UNIT DESCRIPTION

  (c) 1990 Vladimir Batagelj

  version 0:  november 1990
  version 1:  26. nov  1990  zacelo delati / lokalna optimizacija
                             ave_mean/line_ave
  version 2:  29. nov  1990  predelava BLOCK v REGULA
  version 3:   5. dec  1990  top ten
  version 4:   9. dec  1990  premene
  version 5:  14. may  1993  ANA file
  version 6:  29. may  1993  zacetne razvrstitve, ask, glava strani
                             MODEL
  version 7:   5. jun  1993  average, median
  version 8:   7. jun  1993  active -> off, dia, outdia, all
  version 9:  10. oct  1993  dia/out-density
  version 10: 16. jul  1994  kanonski zapis razvrstitve
  version 11: 27. jan  1996  modelska matrika
  version 12: 28. jan  1996  omejitve
  version 13: 17. mar  1996  drobni popravki, don't care
  version 14: 16. may  1998  razlicne kazni za posamezne omejitve
  version 15: 24. may  1998  velike matrike
  version 16:  7. aug  1998  non-null, symmetric;
                             odstrani kanonske (tezave pri omejitvah)
  version 17: 18. aug  1998  numax=130, pem=50
                             omejitev  i notin Cj

  Ideje:
  16. jul 1994
     - omejitev na stevilo skupin z eno enoto
     - tabu z utezmi
     - kazenski parametri za neizpolnjevanje omejitev
     - nastavitveni zaslon
     - kanonski zapis razvrstitve
  18. maj 1998
     - tip  not null
     - omejitvi min / max na stevilo enot v posamezni skupini
     - velike matrike
     - seznamska predstavitev skupin
}

LABEL 14;

CONST
   numax  = 130;  { max # of units }
   clmax  = 30;   { max # of clusters }
   bignum = 1E10; { infinity }
   tmax   = 14;   { max # of types }
   pem    = 50;   { max # of constraints }

TYPE
   vertex   = 0 .. numax;
   partition = ARRAY [1..clmax] OF SET OF 1..numax;
   element   = RECORD ib: LongInt; pb: real; sb: partition END;
   heading   = String[40];
   uname  = STRING[10];
   tipi   = SET OF 0 .. tmax;
   index  = ARRAY [ 1 .. numax ] OF byte;
   names  = ARRAY [ 0 .. numax ] OF uname;
{   matrix = ARRAY [ 1 .. numax, 1 .. numax ] OF integer; }
   matrix = ARRAY [ 0 .. 0 ] OF integer;
   modmtx = ARRAY [ 1 .. clmax, 1 .. clmax ] OF integer;
   tipmtx = ARRAY [ 1 .. clmax, 1 .. clmax ] OF tipi;
{   table  = ARRAY [ 1 .. numax*numax ] OF integer; }
   table  = ARRAY [ 0 .. 0 ] OF integer;
{   relation = ARRAY [ 1 .. numax, 1 .. numax ] OF byte; }
   relation = ARRAY [ 0 .. 0 ] OF byte;
   vector = ARRAY [ 0 .. tmax ] OF real;
   counters = ARRAY [ 0 .. numax ] OF integer;

VAR
   w : vector;
   tab : ^table;
   active : ARRAY[ 0 .. tmax ] OF boolean;
   prior  : ARRAY[ 0 .. tmax ] OF byte;
   seq    : ARRAY[ 0 .. tmax+1 ] OF integer;
   tina   : ARRAY[ 0 .. tmax ] OF string[15];
   til    : ARRAY[ 0 .. tmax ] OF string[ 3];
   head : heading;
   num,  cnt, ibest, cnt_max : LongInt ;
   s_p : partition;
   x: element;
   top_ten: ARRAY [0..10] OF element;
   t, tit: string;
   c, s, z, cn, sn, zn, ind, p_p : index;
   pet, pel, pec : index; {constraints}
   penality : counters;
   pen, penal, leng, m : integer;
   penerr : real;
   name : names;
   a : ^matrix;
   block_mod, block_val, block_err, block_pen : modmtx;
   block_ave : ARRAY [ 1 .. clmax, 1 .. clmax ] OF real;
   block_tip : tipmtx;
   tipij : tipi;
   rel : ^relation;
   c_read, test, all, fd, ok, continue, found, RandomCls,
   symm, ask, save, compute, transpose, finish, ReadANA, c_list  : boolean;
   mf, md, c_max, c_new, s_max, r_max, z_max, NumUnits, debug, max_t,
   t_print, sim, step, u, i, j, k, ij, p, q, tmp, nnn, nw, v, lc, ict, er,
   limit, nclust, indana, savelev, select, selval, form, clm : integer;
   c_total, s_total, pbest, density, diadensity, outdensity : real;
   ch : char;
   kontra, keyclu : string[5];
   prgnam : string[10];
   title : ARRAY[0..1] OF String[70];
   netit : string[70];
   filnet, fillst, filana, filmdl: String[40];
   ana, net, mdl: text;
{
   c[i] -- i-th unit belong to cluster c[i]
   b[i] -- i-th unit belong to cluster b[i]
   bn[i] = # of units in the i-th cluster of  b
}

FUNCTION max(a,b:integer):integer;
BEGIN
   IF a < b THEN max := b ELSE max := a;
END { max };

FUNCTION min(a,b:integer):integer;
BEGIN
   IF a < b THEN min := a ELSE min := b;
END { min };

FUNCTION amax(a,b:real):real;
BEGIN
   IF a < b THEN amax := b ELSE amax := a;
END { amax };

PROCEDURE sort( VAR tab: table; len: integer );
VAR i, j, k, tmp : integer;
BEGIN
   j := len; k := 1;
   WHILE j >1 DO BEGIN
      FOR i := 1 TO j-1 DO IF tab[i] > tab[j] THEN
         BEGIN tmp := tab[i]; tab[i] := tab[j]; tab[j] := tmp; k := j END;
      j := k; k := 1;
   END;
END;

PROCEDURE read_names;
VAR k, i: integer; ch: char; un: uname;
BEGIN
   k := 0;
   WHILE k < NumUnits DO BEGIN
      IF eof(net) THEN BEGIN {k := NumUnits + 1}
         ClrScr; GotoXY(20,15); write('Error in NAMES');
         readln; Halt;
      END ELSE BEGIN
         k := k + 1; un := '';
         FOR i := 1 TO 10 DO BEGIN
            IF eoln(net) THEN ch := ' ' ELSE read(net,ch);
            un := un + ch;
         END;
         name[k] := un;
         IF debug > 3 THEN writeln( k:5, ' - ', un );
         IF k MOD 7 = 0 THEN readln(net);
      END;
   END;
END;

PROCEDURE Read_Network;
LABEL 14;
VAR num, k, i, j, r, nvtxp : integer; ch:char;
    u, v : vertex;
BEGIN
   IF (form < 0) OR (form > 3) OR (form = 2) THEN Error_Cl(6);
   IF form < 2 THEN BEGIN
      ij := 0;
      FOR i := 1 TO NumUnits DO FOR j := 1 TO NumUnits DO BEGIN
         ij := ij + 1; read(net, a^[ij]);
      END;
      symm := FALSE;
   END ELSE IF form = 3 THEN BEGIN
      nvtxp := NumUnits + 1; symm := sim <> 0;
      mf := 2; md := 0;
      writeln(lst,'Network with ',NumUnits:1, ' vertices');
      ij := 0;
      FOR u := 1 TO NumUnits DO FOR v := 1 TO NumUnits DO
         BEGIN ij := ij + 1; a^[ij] := 0 END;
      v := 0;
      WHILE NOT eof(net) DO BEGIN
         read(net,num);
         IF abs(num) > NumUnits THEN BEGIN
            writeln(lst);
            writeln(lst,'number ',num:1,' in network out of range');
            writeln('number ',num:1,' in network out of range');
            close(lst); readln; halt;
         END;
         IF num = 0 THEN GOTO 14;
         IF num < 0 THEN BEGIN
            v := -num; k := 0;
            IF t_print > 0 THEN
               BEGIN writeln(lst); write(lst,'vertex',v:4,' : ') END;
         END ELSE BEGIN
            k := k + 1;
            IF (k > 1) AND (k MOD 10 = 1) AND (t_print > 0) THEN
               BEGIN writeln(lst); write(lst,' ':13) END;
            IF t_print > 0 THEN write(lst,num:4);
            a^[(v-1)*NumUnits+num] := 1; IF symm THEN a^[(num-1)*NumUnits+v] := 1;
         END;
      END;
   END;
14:
   writeln(lst); ij := 0;
   FOR i := 1 TO NumUnits DO FOR j := 1 TO NumUnits DO BEGIN
      ij := ij + 1;
      IF a^[ij] > 0 THEN rel^[ij] := 1 ELSE rel^[ij] := 0;
   END;
   FOR i := 1 TO NumUnits DO ind[i] := i;
END {read_network};

PROCEDURE Read_Network_Old;
VAR i,j,r:integer; ch:char;
BEGIN
   REPEAT read(net,ch); write(ch) UNTIL ch = ' ';
   REPEAT read(net,ch); write(ch) UNTIL ch = ' ';
   writeln; readln(net, NumUnits, r_max);
   readln(net, title[0]);
   Read_Names;
   r_max := 1;                                                { <--- **** }
   FOR r := 1 TO r_max DO BEGIN
      REPEAT readln(net,title[r]) UNTIL title[r][1] = '*';
      readln(net,title[r]); writeln(title[r]); ij := 0;
      FOR i := 1 TO NumUnits DO FOR j := 1 TO NumUnits DO BEGIN
         ij := ij+1; read(net, a^[ij]);
         IF a^[ij] > 0 THEN rel^[ij] := 1 ELSE rel^[ij] := 0;
      END;
      readln(net);
   END;
   FOR i := 1 TO NumUnits DO ind[i] := i;
END;

PROCEDURE Print_Matrix(VAR a:matrix; VAR ind:index);
VAR i,j: integer;
BEGIN
   newlin(-(NumUnits+3),0); newlin(3,0);
   writeln(lst); writeln(lst,'MATRIX : ', title[1]); writeln(lst);
   write(lst,' ':14);
   FOR i := 1 TO NumUnits DO write(lst,ind[i]:3); writeln(lst);
   FOR i := 1 TO NumUnits DO BEGIN
      write(lst, name[ind[i]], ind[i]:3, ' ');
      FOR j := 1 TO NumUnits DO write(lst,a[(ind[i]-1)*NumUnits+ind[j]]:3);
      writeln(lst); newlin(1,0);
   END;
   writeln(lst);
END;

PROCEDURE Print_Network(ind:index);
VAR r: integer;
BEGIN
   FOR r := 1 TO r_max DO Print_Matrix(a^,ind);
END;

PROCEDURE Print_ModMatrix;
VAR i,j,r: integer;
BEGIN
   r := 1;
   writeln(lst); writeln(lst,'MODEL MATRIX : ', title[r]); writeln(lst);
   newlin(4,0); write(lst,' ':6);
   FOR i := 1 TO c_max DO write(lst,i:4); writeln(lst);
   FOR i := 1 TO c_max DO BEGIN
      newlin(1,0); write(lst, i:6, ' ');
      FOR j := 1 TO c_max DO write(lst,til[block_mod[i,j]]:4);
      writeln(lst);
   END;
   writeln(lst);
END;

PROCEDURE Print_ErrMatrix;
VAR i,j,r: integer;
BEGIN
   r := 1;
   writeln(lst); writeln(lst,'ERROR MATRIX : ', title[r]); writeln(lst);
   newlin(4,0); write(lst,' ':7);
   FOR i := 1 TO c_max DO write(lst,i:4); writeln(lst);
   FOR i := 1 TO c_max DO BEGIN
      newlin(1,0); write(lst, i:6, ' ');
      FOR j := 1 TO c_max DO write(lst,block_err[i,j]:4);
      writeln(lst);
   END;
   writeln(lst);
END;

PROCEDURE Print_ValMatrix;
VAR i,j,r: integer;
BEGIN
   r := 1;
   writeln(lst); writeln(lst,'VALUE MATRIX : ', title[r]); writeln(lst);
   newlin(4,0); write(lst,' ':7);
   FOR i := 1 TO c_max DO write(lst,i:4); writeln(lst);
   FOR i := 1 TO c_max DO BEGIN
      newlin(1,0); write(lst, i:6, ' ');
      FOR j := 1 TO c_max DO write(lst,' ',block_val[i,j]:3);
      writeln(lst);
   END;
   writeln(lst);
END;

PROCEDURE Print_AveMatrix;
VAR i,j,r: integer;
BEGIN
   r := 1;
   writeln(lst); writeln(lst,'VALUE MATRIX : ', title[r]); writeln(lst);
   newlin(4,0); write(lst,' ':7);
   FOR i := 1 TO c_max DO write(lst,i:9); writeln(lst);
   FOR i := 1 TO c_max DO BEGIN
      newlin(1,0); write(lst, i:6, ' ');
      FOR j := 1 TO c_max DO write(lst,' ',block_ave[i,j]:8:3);
      writeln(lst);
   END;
   writeln(lst);
END;

PROCEDURE Print_Values;
VAR r: integer;
BEGIN
   FOR r := 1 TO r_max DO Print_ValMatrix;
END;

{ ======================================================= }

PROCEDURE block_error(ic,jc:integer; VAR d: counters; VAR p_err: vector);
VAR i, j, ib, jb, r, k, block : integer;
   sum, suma, prod, delta, diff, syme : integer;
   delt, average : real;
   dsum, tsum, diag, cmax, cmin, rmax, rmin, vsum,
   rpoz, cpoz, rnul, cnul, nr, nc, num,
   median, ind, jnd : integer;
   q : vector;
   rsum, csum : counters;

BEGIN
   block_mod[ic,jc] := 0;
   tipij := block_tip[ic,jc];
   FOR i := 0 TO tmax DO active[i] := i in tipij;

{ compute block (ic,jc) characteristics:
  tsum  - sum of block elements
  dsum  - sum of diagonal elements
  nr    - # of rows
  rpoz  - # of non-null rows
  rmin  - min row sum
  rmax  - max row sum
  rsum  - row sums
  nc    - # of cols
  cpoz  - # of non-null cols
  cmin  - min col sum
  cmax  - max col sum
  csum  - col sums
  syme  - symmetric error
}
   nr := 0; nc := 0; dsum := 0; tsum := 0; rpoz := 0; cpoz := 0;
   rmax := 0; cmax := 0; rmin := 999; cmin := 999; syme := 0;
   FOR i := 1 TO NumUnits DO IF c[i] = ic THEN BEGIN
      sum := 0; nr := nr + 1;
      FOR j := 1 TO NumUnits DO IF c[j] = jc THEN BEGIN
         sum := sum + rel^[(i-1)*NumUnits+j];
         syme := syme + abs(rel^[(i-1)*NumUnits+j] - rel^[(j-1)*NumUnits+i]);
         IF i = j THEN dsum := dsum + rel^[(i-1)*NumUnits+j];
      END;
      IF sum > 0 THEN rpoz := rpoz + 1;
      IF rmax < sum THEN BEGIN rmax := sum; ind := i END;
      IF rmin > sum THEN rmin := sum;
      rsum[nr] := sum; tsum := tsum + sum;
   END;
   IF ic=jc THEN syme := syme DIV 2;
   FOR j := 1 TO NumUnits DO IF c[j] = jc THEN BEGIN
      sum := 0; nc := nc + 1;
      FOR i := 1 TO NumUnits DO IF c[i] = ic THEN
         sum := sum + rel^[(i-1)*NumUnits+j];
      IF sum > 0 THEN cpoz := cpoz + 1;
      IF cmax < sum THEN BEGIN cmax := sum; jnd := j END;
      IF cmin > sum THEN cmin := sum;
      csum[nc] := sum;
   END;
   IF test THEN BEGIN
      writeln(lst, ic:3, nr:3, rmin:3, rmax:3, rpoz:3, dsum:5 );
      writeln(lst, jc:3, nc:3, cmin:3, cmax:3, cpoz:3, tsum:5 );
   END;
   rnul := nr - rpoz; cnul := nc - cpoz;
   suma := cn[ic] + cn[jc]; prod := cn[ic] * cn[jc];
   IF ic = jc THEN BEGIN
      IF 2*dsum <= nr THEN diag := dsum ELSE diag := nr-dsum;
      density := diadensity
   END ELSE density := outdensity;

{ 0 null }
   IF active[0] THEN BEGIN
      delta := tsum;
      IF ic = jc THEN delta := delta + diag - dsum;
      d[0] := delta;
      IF select = 2 THEN p_err[0] := (1+delta)*w[0]*block_pen[ic,jc]/prod
         ELSE p_err[0] := delta*w[0]*block_pen[ic,jc];
   END ELSE p_err[0] := 999999.;

{ 1 complete }
   IF active[1] THEN BEGIN
      IF ic=jc THEN BEGIN
         delta := prod - tsum + dsum + diag - nr; d[1] := delta;
         IF select = 2 THEN p_err[1] := (1+delta)*w[1]*block_pen[ic,jc]/prod
            ELSE p_err[1] := delta*w[1]*block_pen[ic,jc];
      END ELSE BEGIN
         delta := prod - tsum;
         IF ic = jc THEN delta := delta + dsum + diag - nr;
         d[1] := delta;
         IF select = 2 THEN p_err[1] := (1+delta)*w[1]*block_pen[ic,jc]/prod
            ELSE p_err[1] := delta*w[1]*block_pen[ic,jc];
      END;
   END ELSE p_err[1] := 999999.;

{ 2 in-dominant }
   IF active[2] AND (limit <= nc) THEN BEGIN
      IF ic=jc THEN BEGIN
         delta := (nc - rmax)*nr;
         IF dsum = 0 THEN delta := delta - nr;
         d[2] := delta;
         IF select = 2 THEN p_err[2] := (1+delta)*w[2]*block_pen[ic,jc]/prod
            ELSE p_err[2] := delta*w[2]*block_pen[ic,jc];
      END ELSE BEGIN
         delta := (nc - rmax)*nr; d[2] := delta;
         IF select = 2 THEN p_err[2] := (1+delta)*w[2]*block_pen[ic,jc]/prod
            ELSE p_err[2] := delta*w[2]*block_pen[ic,jc];
      END;
   END ELSE p_err[2] := 999999.;

{ 3 out-dominant }
   IF active[3] AND (limit <= nr) THEN BEGIN
      IF ic=jc THEN BEGIN
         delta := (nr - cmax)*nc;
         IF dsum = 0 THEN delta := delta - nc;
         d[3] := delta;
         IF select = 2 THEN p_err[3] := (1+delta)*w[3]*block_pen[ic,jc]/prod
            ELSE p_err[3] := delta*w[3]*block_pen[ic,jc];
      END ELSE BEGIN
         delta := (nr - cmax)*nc; d[3] := delta;
         IF select = 2 THEN p_err[3] := (1+delta)*w[3]*block_pen[ic,jc]/prod
            ELSE p_err[3] := delta*w[3]*block_pen[ic,jc];
      END;
   END ELSE p_err[3] := 999999.;

{ 4 regular }
   IF active[4] THEN BEGIN
      delta := cnul*nr + rnul*cpoz; d[4] := delta;
      IF select = 2 THEN p_err[4] := (1+delta)*w[4]*block_pen[ic,jc]/prod
         ELSE p_err[4] := delta*w[4]*block_pen[ic,jc];
   END ELSE p_err[4] := 999999.;

{ 5 in-regular }
   IF active[5] THEN BEGIN
      delta := rnul*nc;
      IF nr < limit THEN p_err[5] := 999999.
      ELSE BEGIN
         d[5] := delta;
         IF select = 2 THEN p_err[5] := (1+delta)*w[5]*block_pen[ic,jc]/prod
            ELSE p_err[5] := delta*w[5]*block_pen[ic,jc];
      END;
   END ELSE p_err[5] := 999999.;

{ 6 out-regular }
   IF active[6] THEN BEGIN
      delta := cnul*nr;
      IF nc < limit THEN p_err[6] := 999999.
      ELSE BEGIN
         d[6] := delta;
         IF select = 2 THEN p_err[6] := (1+delta)*w[6]*block_pen[ic,jc]/prod
            ELSE p_err[6] := delta*w[6]*block_pen[ic,jc];
      END;
   END ELSE p_err[6] := 999999.;

{ 7 in-function }
   IF active[7] AND (limit <= nr) THEN BEGIN
      delta := tsum - rpoz + rnul*nc; d[7] := delta;
      IF select = 2 THEN p_err[7] := (1+delta)*w[7]*block_pen[ic,jc]/prod
         ELSE p_err[7] := delta*w[7]*block_pen[ic,jc];
   END ELSE p_err[7] := 999999.;

{ 8 out-function }
   IF active[8] AND (limit <= nc) THEN BEGIN
      delta := tsum - cpoz + cnul*nr; d[8] := delta;
      IF select = 2 THEN p_err[8] := (1+delta)*w[8]*block_pen[ic,jc]/prod
         ELSE p_err[8] := delta*w[8]*block_pen[ic,jc];
   END ELSE p_err[8] := 999999.;

{ 9 density }
   IF active[9] THEN BEGIN
      IF ic = jc THEN BEGIN
         IF dsum > 0 THEN delta := Round(0.5+density*prod)
            ELSE delta := Round(0.5+density*(prod-nr));
         delta := delta - tsum;
      END ELSE delta := Round(0.5+density*prod) - tsum;
      IF delta < 0 THEN delta := 0;
      d[9] := delta;
      IF select = 2 THEN p_err[9] := (1+delta)*w[9]*block_pen[ic,jc]/prod
         ELSE p_err[9] := delta*w[9]*block_pen[ic,jc];
   END ELSE p_err[9] := 999999.;

{ 10 do not care }
   IF active[10] THEN BEGIN
      d[10] := 0; p_err[10] := 0;
   END ELSE p_err[10] := 999999.;

{ 11 non-null }
   IF active[11] THEN BEGIN
      delta := 1;
      IF tsum > 0 THEN delta := 0;
      d[11] := delta;
      IF select = 2 THEN p_err[11] := (1+delta)*w[11]*block_pen[ic,jc]/prod
         ELSE p_err[11] := delta*w[11]*block_pen[ic,jc];
   END ELSE p_err[11] := 999999.;

{ 12 symmetric }
   IF active[12] THEN BEGIN
      delta := syme;
      d[12] := delta;
      IF select = 2 THEN p_err[12] := (1+delta)*w[12]*block_pen[ic,jc]/prod
         ELSE p_err[12] := delta*w[12]*block_pen[ic,jc];
   END ELSE p_err[12] := 999999.;

END { block_error };


PROCEDURE Compute_value( ic, jc, block : integer );
VAR i, j, ib, jb, r, k, nr, nc, prod, poz : integer;
   sum, suma, delta, diff : integer;
   delt, average : real;
   dsum, tsum, diag, cmax, cmin, rmax, rmin, vsum,
   rpoz, cpoz, rnul, cnul, num,
   median, ind, jnd : integer;
   q, p_err : vector;
   rsum, csum, d : counters;
BEGIN
   IF selval <> 1 THEN BEGIN
         block_val[ic,jc] := 14;
   END ELSE BEGIN
      IF selval = 1 THEN BEGIN {average}
      { compute block (ic,jc) characteristics }
         vsum := 0; poz := 0; nc := 0; nr := 0;
         FOR i := 1 TO NumUnits DO IF c[i] = ic THEN BEGIN
            nr := nr + 1;
            FOR j := 1 TO NumUnits DO IF c[j] = jc THEN BEGIN
               nc := nc + 1;
               vsum := vsum + a^[(i-1)*NumUnits+j];
               IF a^[(i-1)*NumUnits+j] > 0 THEN poz := poz + 1;
            END;
         END;
         prod := nc*nr;
         CASE block OF
         0 : average := vsum/prod;
         1 : IF ic <> jc THEN average := vsum/prod
             ELSE IF dsum > 0 THEN average := vsum/prod
             ELSE IF prod > nr THEN average := vsum/(prod-nr)
             ELSE average := 0.;
         2 : average := vsum/nc;
         3 : average := vsum/nr;
         4 : average := vsum/max(nr,nc);
         5 : average := vsum/nr;
         6 : average := vsum/nc;
         7 : average := vsum/nr;
         8 : average := vsum/nc;
         9 : IF ic <> jc THEN average := vsum/prod/outdensity
             ELSE IF dsum > 0 THEN average := vsum/prod/diadensity
             ELSE average := vsum/(prod-nr)/diadensity;
        10 : average := vsum/prod;
        11 : average := vsum;
        12 : IF poz > 0 THEN average := vsum/poz ELSE average := 0;
         END;
         block_ave[ic,jc] := average;
      END ELSE IF selval = 2 THEN BEGIN {median}
      { compute block (ic,jc) characteristics }
         CASE block OF
         0 : median := 0;
         1, 10 : BEGIN { complete }
                k := 0;
                FOR i := 1 TO NumUnits DO IF c[i] = ic THEN
                   FOR j := 1 TO NumUnits DO IF c[j] = jc THEN
                      BEGIN k := k+1; tab^[k] := a^[(i-1)*NumUnits+j] END;
                sort(tab^,k);
                IF ic <> jc THEN median := tab^[1+k DIV 2]
                ELSE IF dsum > 0 THEN median := tab^[1+k DIV 2]
                ELSE median := tab^[(1+k+nr) DIV 2];
             END;
         2 : BEGIN { in-dominant }
                k := 0;
                FOR j := 1 TO NumUnits DO IF c[j] = jc THEN
                   BEGIN k := k+1; tab^[k] := a^[(ind-1)*NumUnits+j] END;
                sort(tab^,k);
                IF ic <> jc THEN median := tab^[1+k DIV 2]
                ELSE IF dsum > 0 THEN median := tab^[1+k DIV 2]
                ELSE median := tab^[(1+k+1) DIV 2];
             END;
         3 : BEGIN { out-dominant }
                k := 0;
                FOR i := 1 TO NumUnits DO IF c[i] = ic THEN
                   BEGIN k := k+1; tab^[k] := a^[(i-1)*NumUnits+jnd] END;
                sort(tab^,k);
                IF ic <> jc THEN median := tab^[1+k DIV 2]
                ELSE IF dsum > 0 THEN median := tab^[1+k DIV 2]
                ELSE median := tab^[(1+k+1) DIV 2];
             END;
         4 : BEGIN { regular }
                k := 0;
                FOR i := 1 TO NumUnits DO IF c[i] = ic THEN BEGIN
                   rmax := 0;
                   FOR j := 1 TO NumUnits DO IF c[j] = jc THEN
                      IF rmax < a^[(i-1)*NumUnits+j] THEN rmax := a^[(i-1)*NumUnits+j];
                   k := k+1; tab^[k] := rmax
                END;
                FOR j := 1 TO NumUnits DO IF c[j] = jc THEN BEGIN
                   cmax := 0;
                   FOR i := 1 TO NumUnits DO IF c[i] = ic THEN
                      IF cmax < a^[(i-1)*NumUnits+j] THEN cmax := a^[(i-1)*NumUnits+j];
                   k := k+1; tab^[k] := cmax
                END;
                sort(tab^,k); median := tab^[1+k DIV 2];
             END;
         5 : BEGIN { in-regular }
                k := 0;
                FOR i := 1 TO NumUnits DO IF c[i] = ic THEN BEGIN
                   rmax := 0;
                   FOR j := 1 TO NumUnits DO IF c[j] = jc THEN
                      IF rmax < a^[(i-1)*NumUnits+j] THEN rmax := a^[(i-1)*NumUnits+j];
                   k := k+1; tab^[k] := rmax
                END;
                sort(tab^,k); median := tab^[1+k DIV 2];
             END;
         6 : BEGIN { out-regular }
                k := 0;
                FOR j := 1 TO NumUnits DO IF c[j] = jc THEN BEGIN
                   cmax := 0;
                   FOR i := 1 TO NumUnits DO IF c[i] = ic THEN
                      IF cmax < a^[(i-1)*NumUnits+j] THEN cmax := a^[(i-1)*NumUnits+j];
                   k := k+1; tab^[k] := cmax
                END;
                sort(tab^,k); median := tab^[1+k DIV 2];
             END;
         7 : BEGIN { in-function }
                k := 0;
                FOR i := 1 TO NumUnits DO IF c[i] = ic THEN BEGIN
                   rmax := 0;
                   FOR j := 1 TO NumUnits DO IF c[j] = jc THEN
                      IF rmax < a^[(i-1)*NumUnits+j] THEN rmax := a^[(i-1)*NumUnits+j];
                   k := k+1; tab^[k] := rmax
                END;
                sort(tab^,k); median := tab^[1+k DIV 2];
             END;
         8 : BEGIN { out-function }
                k := 0;
                FOR j := 1 TO NumUnits DO IF c[j] = jc THEN BEGIN
                   cmax := 0;
                   FOR i := 1 TO NumUnits DO IF c[i] = ic THEN
                      IF cmax < a^[(i-1)*NumUnits+j] THEN cmax := a^[(i-1)*NumUnits+j];
                   k := k+1; tab^[k] := cmax
                END;
                sort(tab^,k); median := tab^[1+k DIV 2];
             END;
         9 : BEGIN { density }
                k := 0;
                FOR i := 1 TO NumUnits DO IF c[i] = ic THEN
                   FOR j := 1 TO NumUnits DO IF c[j] = jc THEN
                      BEGIN k := k+1; tab^[k] := a^[(i-1)*NumUnits+j] END;
                sort(tab^,k);
                IF dsum > 0 THEN median := tab^[k - Trunc(density*k/2)]
                ELSE median := tab^[(1+k+nr+Round((1-density))*(k-nr)) DIV 2]
             END;
         END;
         block_val[ic,jc] := median;
      END;
   END;
END { Compute_Value};


FUNCTION Pair_Block_Error(ic,jc: integer): real;
VAR i, j, ib, jb, r, k, block, blocx : integer;
   sum, suma, prod, delta, diff, syme : integer;
   delt, delx, average : real;
   dsum, tsum, diag, cmax, cmin, rmax, rmin, vsum,
   rpoz, cpoz, rnul, cnul, nr, nc, num,
   median, ind, jnd : integer;
   q, era, erb : vector;
   rsum, csum, da, db : counters;

BEGIN

   block_error(ic,jc,  da, era );
   delt := 999998.; block := -1;
   FOR i := 1 TO max_t+1 DO BEGIN
      k := seq[i];
      IF k >= 0 THEN IF era[k] < delt THEN
         BEGIN delt := era[k]; block := k END;
   END;
   IF test THEN
      writeln(lst,' error[ ',ic:2,', ',jc:2,' ] = ',delt:10:3, block:5);
   IF ic = jc THEN BEGIN { diagonal block }
      block_mod[ic,jc] := block; block_err[ic,jc] := da[block];
      Pair_Block_Error := delt;
   END ELSE BEGIN { non-diagonal blocks }
      block_error(jc,ic,  db, erb );

      delx := 999998.; blocx := -1;
      FOR i := 1 TO max_t+1 DO BEGIN
         k := seq[i];
         IF k >= 0 THEN IF erb[k] < delx THEN
            BEGIN delx := erb[k]; blocx := k END;
      END;
      IF 2*era[12] < era[block]+erb[blocx] THEN
         BEGIN block := 12; blocx := 12 END;

      block_mod[ic,jc] := block; block_err[ic,jc] := da[block];
      block_mod[jc,ic] := blocx; block_err[jc,ic] := db[blocx];
      Pair_Block_Error := era[block]+erb[blocx];
      IF test THEN
         writeln(lst,' error[ ',jc:2,', ',ic:2,' ] = ',delx:10:3, blocx:5);
   END;
   IF compute THEN BEGIN
      Compute_value( ic, jc, block );
      IF ic <> jc THEN Compute_value( jc, ic, blocx );
   END;

END { Pair_Block_Error };

FUNCTION total_error: real;
VAR ic, jc, i : integer;
    total: real;
BEGIN
 { errors in pairs of blocks }
   total := 0;
   FOR ic := 1 TO c_max DO FOR jc := ic TO c_max DO
      total := total + pair_block_error(ic,jc);
   penerr := 0;
 { add penalties }
   FOR i := 1 TO pen DO
      CASE pet[i] OF
       1 : IF c[pel[i]] <> pec[i]    THEN penerr := penerr + penality[i];
       2 : IF c[pel[i]] =  pec[i]    THEN penerr := penerr + penality[i];
       3 : IF c[pel[i]] <> c[pec[i]] THEN penerr := penerr + penality[i];
       4 : IF c[pel[i]] =  c[pec[i]] THEN penerr := penerr + penality[i];
       5 : IF cn[pec[i]] < pel[i]    THEN penerr := penerr +
                                     (pel[i] - cn[pec[i]])*penality[i];
       6 : IF cn[pec[i]] > pel[i]    THEN penerr := penerr +
                                     (cn[pec[i]] - pel[i])*penality[i];
      END;
   total_error := total + penerr;
END { total_error };

PROCEDURE list(s: partition);
VAR i,j: integer;
BEGIN
   FOR i:=1 TO c_max DO BEGIN
      write(lst,' { ');
      FOR j:=1 TO NumUnits DO BEGIN
         IF j IN s[i] THEN write(lst,j:1,' ')
      END;
      write(lst,'}')
   END;
   write(lst,'  ');
END; { list }

PROCEDURE printresult;
VAR i : integer;
BEGIN
   newlin(-15,0); newlin(4,0);
   writeln(lst); writeln(lst,' T O P   T E N'); writeln(lst);
   FOR i := 1 TO 10 DO WITH top_ten[i] DO
      IF pb < bignum THEN BEGIN
         newlin(1,0); write(lst,ib:6,':');
         list(sb); writeln(lst,'P = ',pb:8:3)
      END;
   writeln(lst);
   s_p := top_ten[1].sb;
   k := 0;
   FOR j := 1 TO s_max DO BEGIN
      NewLin(1,0); writeln(lst,'CLUSTER ',j:1); sn[j] := 0;
      m := 0;
      FOR i := 1 TO NumUnits DO IF i IN s_p[j] THEN BEGIN
         k := k+1; ind[k] := i; s[i] := j;
         m := m+1; write(lst,i:5,' ',name[i]);
         IF m MOD 5 = 0 THEN
            BEGIN NewLin(1,0); writeln(lst) END;
      END;
      IF m MOD 5 <> 0 THEN
         BEGIN NewLin(1,0); writeln(lst) END;
   END;
   FOR i := 1 TO NumUnits DO sn[s[i]] := sn[s[i]] + 1;
   c := s; cn := sn;
   Print_Network(ind);
   compute := TRUE;
   s_total := total_error;
   Print_ModMatrix; Print_ErrMatrix;
   IF selval = 1 THEN Print_AveMatrix;
   IF selval = 2 THEN Print_ValMatrix;
   writeln(lst,'final error =', s_total:12:3);
   compute := FALSE;
END; { printresult }

PROCEDURE Read_Clustering;
VAR prgnam, memode: String;
BEGIN
   REPEAT
      indana := 0; kontra := '';
      REPEAT
     IF eof(ana) THEN BEGIN
            ReadANA := FALSE; close(ana); append(ana);
     END ELSE BEGIN
        read(ana,ch);
        IF ch = '*' THEN BEGIN
           kontra := '';
           FOR i := 1 TO 5 DO BEGIN
              IF eoln(ana) THEN ch := ' ' ELSE read(ana,ch);
          kontra := Concat(kontra,ch);
           END;
        END;
        readln(ana);
     END;
      UNTIL (NOT ReadANA) OR (kontra=keyclu) ;
      IF kontra=keyclu THEN BEGIN
         readln(ana,c_tit); prgnam := '';
     FOR i := 1 TO 10 DO BEGIN
        IF eoln(ana) THEN ch := ' ' ELSE read(ana,ch);
        prgnam := Concat(prgnam,ch);
     END;
     memode := '';
     FOR i := 1 TO 10 DO BEGIN
        IF eoln(ana) THEN ch := ' ' ELSE read(ana,ch);
        memode := Concat(memode,ch);
     END;
     readln(ana, nnn, nclust );
         ClrScr; HighVideo;
         gotoxy(15,7);  write('***  S T R A N  -  MODEL 2');
         gotoxy(15,8);  write('blockmodeling by local optimization');
         NormVideo;
         gotoxy(15,10); write( 'Reading new clustering from ANA');
     gotoxy(15,11); write( c_tit );
     gotoxy(15,12); write( nclust:5, ' clusters,  program: ', prgnam,
            ',  method: ', memode );
     gotoxy(15,13); write( 'choose ANA-option' );
     gotoxy(15,14); write( '   1 - optimize the clustering' );
     gotoxy(15,15); write( '   2 - skip the clustering' );
     gotoxy(15,16); write( '   3 - stop *.ANA file search' );
         indana := ReadInt(15,17,'enter your choice = ',1,3);
     IF indana <= 1 THEN BEGIN indana := 1; c_new := -1 END
     ELSE IF indana > 3 THEN indana := 3;
     ReadANA := indana < 3;
     IF indana = 3 THEN
            BEGIN close(ana); append(ana); indana := 0 END;
         IF indana = 1 THEN BEGIN
            write( ' Read ANA,' );
        FOR i := 1 TO NumUnits DO read(ana, z[i]);
            readln(ana); cnt_max := cnt + 1
         END;
      END;
   UNTIL indana <= 1;
   Newlin(100,0);
END;

PROCEDURE Save_Clustering;
BEGIN
   IF ask THEN BEGIN
      GoToXY(1,21); ClrEol; HighVideo; Write('Saving clustering:');
      NormVideo; Write(c_tit); HighVideo;
      GoToXY(1,22); Write('New title or <ENTER> : '); NormVideo;
      GoToXY(1,23); Readln(t);
      IF Length(t) = 0 THEN t := c_tit;
      IF Length(t) > 70 THEN t := Copy( t, 1,70 );
      FOR i := 21 TO 23 DO BEGIN GoToXY(1,i); ClrEol END;
   END ELSE t := tit;
   writeln( ana, '*CLUSTERING', ' ':9, c_dat );
   writeln( ana, t );
   writeln( ana, prgnam:10, memode:10, NumUnits:5, s_max:5,
      0:5, 0:5, s_total:15:7 );
   FOR i := 1 TO NumUnits DO BEGIN
      write( ana, s[i]:4 );
      IF i MOD 20 = 0 THEN writeln( ana )
   END;
   writeln( ana );
END;

PROCEDURE Random_Initial;
BEGIN
   cnt_max := ReadInt(15,16,
      'Number of random repetitions              = ',0,maxint-1);
   savelev := ReadInt(15,18,
      'Save: 0 - no, 1 - ask, 2 - opt, 3 - all   = ',0,3);
   save := savelev = 3;
   ClrScr; HighVideo; RandomCls := TRUE; ask := savelev = 1;
   memode := 'RANDOM     '; c_tit := tit;
   gotoxy(20,7);  write('***  S T R A N  -  MODEL 2');
   gotoxy(20,8);  write('blockmodeling by local optimization'); NormVideo;
   gotoxy(20,10); write('NETWORK FILE : ', f_dat );
   gotoxy(20,11); write(NumUnits:1, ' UNITS READ' );
   gotoxy(20,12); write('NUMBER OF CLUSTERS = ', c_new:7);
   gotoxy(20,14); write('CLUSTERING         = ');
   gotoxy(20,15); write('CLUSTERING *       = ');
   gotoxy(20,16); write('P*                 = ');
   IF c_new = c_max THEN BEGIN
      cnt_max := cnt_max + cnt - 1;
      FOR u := 1 TO 10 DO IF top_ten[u].pb < bignum THEN BEGIN
         GotoXY(55,6+u); write(top_ten[u].ib:6,top_ten[u].pb:12:3);
      END;
   END ELSE BEGIN
      IF cnt > 1 THEN BEGIN
         newlin(-5,0); newlin(5,0); writeln(lst);
         writeln(lst,' NUMBER OF UNITS       = ', NumUnits:6 );
         writeln(lst,' NUMBER OF CLUSTERS    = ', c_max:6 );
         writeln(lst,' NUMBER OF CLUSTERINGS = ', cnt-1:6 );
         printresult; Flush(lst);
      END;
      c_max := c_new; { ask := TRUE; }
      cnt := 1; ibest := 0; pbest := 2*bignum; found := FALSE;
      FOR i:=1 TO 10 DO top_ten[i].pb := pbest;
   END;
   gotoxy(41,14); write(ibest:7);
   IF pbest < bignum THEN BEGIN gotoxy(41,16); write(pbest:11:3) END;
END;

PROCEDURE Read_Initial;
BEGIN
   c_read := TRUE; ask := c_new = 0; memode := 'READ      ';
   z_max := 0; ClrScr;
   FOR i := 1 TO NumUnits DO zn[i] := 0;
   IF c_new = 0 THEN
      BEGIN GotoXY(20,14); write('Enter clustering ') END;
   FOR i := 1 TO NumUnits DO BEGIN
      IF c_new = 0 THEN BEGIN
         GotoXY(25,16); write('c[', i:3);
         z[i] := ReadInt(31,16,' ] = ',0,NumUnits);
      END;
      IF z_max < z[i] THEN z_max := z[i];
      zn[z[i]] := zn[z[i]] + 1;
   END;
   GotoXY(25,19); write('NUMBER OF CLUSTERS    = ', z_max:6);
   Delay(500); cnt_max := 0;
   IF z_max <> c_max THEN BEGIN
      IF cnt > 1 THEN BEGIN
         newlin(-5,0); newlin(5,0); writeln(lst);
         writeln(lst,' NUMBER OF UNITS       = ', NumUnits:6 );
         writeln(lst,' NUMBER OF CLUSTERS    = ', c_max:6 );
         writeln(lst,' NUMBER OF CLUSTERINGS = ', cnt-1:6 );
         printresult; Flush(lst);
      END;
      cnt := 1; ibest := 0; pbest := 2*bignum; found := FALSE;
      FOR i:=1 TO 10 DO top_ten[i].pb := pbest;
   END;
   c := z; cn := zn; c_max := z_max;
END;

PROCEDURE Generate_Initial;
BEGIN
   c_read := FALSE;
   CASE cnt MOD 3 OF
   0 : FOR i := 1 TO NumUnits DO
          BEGIN c[i] := i MOD c_max + 1; cn[c[i]] := cn[c[i]] + 1 END;
{ **** 12. jan 1991 - finta: enakomerno velike skupine so slaba
  za~etna razvrstitev }
   1 : BEGIN
          FOR i := 1 TO c_max-1 DO
             BEGIN c[i] := i; cn[i] := cn[i] + 1 END;
          FOR i := c_max TO NumUnits DO
             BEGIN c[i] := c_max; cn[c_max] := cn[c_max] + 1 END;
       END;
   2 : BEGIN
          FOR i := 1 TO c_max DO
             BEGIN c[i] := i; cn[i] := cn[i] + 1 END;
          FOR i := c_max+1 TO NumUnits DO BEGIN
             c[i] := 1+Trunc(c_max*random); cn[c[i]] := cn[c[i]] + 1
          END;
       END;
   END;
   FOR i := NumUnits DOWNTO 2 DO BEGIN
      j := 1 + trunc(i*random); tmp := c[j]; c[j] := c[i]; c[i] := tmp;
   END;
END;

PROCEDURE Local_Optimization;
LABEL 8;
BEGIN
   compute := TRUE; s_total := total_error; compute := FALSE;
   IF c_list THEN BEGIN
      Print_ModMatrix; Print_ErrMatrix;
      IF ReadAna THEN BEGIN
        IF selval = 1 THEN Print_AveMatrix;
        IF selval = 2 THEN Print_ValMatrix;
      END;
   END;
   {*** NewLin(1,0); write(lst,cnt:5,'. initial error =', s_total:12:3);}
   REPEAT
      step := 0; u := 0;
      REPEAT
         u := u MOD NumUnits + 1;
         FOR p := 1 TO s_max DO BEGIN
            q := s[u];
            IF (q <> p) AND (sn[q] > 1) THEN BEGIN
               c := s; cn := sn;
               c[u] := p; cn[p] := cn[p] + 1; cn[q] := cn[q] - 1;
               c_total := total_error;
               IF test THEN Print_ModMatrix;
               IF test THEN writeln(lst,'try ', u:3, ' : ', q:3,
                  ' --> ', p:2, '  error =', c_total:12:3);
               IF c_total < s_total THEN BEGIN
                  s_total := c_total; s := c; sn := cn; step := -p;
                  IF all THEN writeln(lst,'move unit ', u:3,
                     ' from ', q:3, ' to ', p:2, '  error =', s_total:12:3);
                  GotoXY(20,19); write('move unit ', u:3, ' from ', q:3,
                     ' to ', p:2, '      error =', s_total:12:3);
                  IF test THEN
                     FOR j := 1 TO s_max DO BEGIN
                        writeln(lst,'CLUSTER ',j:1,' / ',sn[j]:1);
                        FOR q := 1 TO NumUnits DO IF s[q] = j THEN
                           writeln(lst,q:5,' ',name[q]);
                     END;
               END;
            END;
         END;
         step := step + c_max;
      UNTIL step >= NumUnits*c_max;
      finish := TRUE;
      IF transpose THEN BEGIN
         cn := sn;
         FOR u := 1 TO NumUnits-1 DO FOR v := u+1 TO NumUnits DO
         IF s[u] <> s[v] THEN BEGIN
            c := s; c[u] := s[v]; c[v] := s[u];
            c_total := total_error;
            IF test THEN Print_ModMatrix;
            IF test THEN writeln(lst,'try ', u:3, '/', s[u]:1,
               ' <--> ', v:2, '/', s[v]:1, '  error =', c_total:12:3);
            IF c_total < s_total THEN BEGIN
               s_total := c_total; s := c; sn := cn;
               IF all THEN writeln(lst,'transpose units ', u:3, '/',
                  s[v]:1, ' and ', v:3, '/', s[u]:1,
                  '  error =', s_total:12:3);
               GotoXY(20,19); write('transpose units ', u:3, '/',
                  s[v]:1, ' and ', v:3, '/', s[u]:1,
                  '   error =', s_total:12:3);
               IF test THEN
                  FOR j := 1 TO s_max DO BEGIN
                     writeln(lst,'CLUSTER ',j:1,' / ',sn[j]:1);
                     FOR q := 1 TO NumUnits DO IF s[q] = j THEN
                        writeln(lst,q:5,' ',name[q]);
                  END;
               finish := FALSE; GOTO 8;
            END;
         END;
      END;
    8:
   UNTIL finish;

   k := 0;
   FOR i := 1 TO s_max DO s_p[i] := [];
   FOR j := 1 TO s_max DO BEGIN
      IF c_list THEN BEGIN
         NewLin(1,0); writeln(lst,'CLUSTER ',j:1,' / ',sn[j]:1);
      END;
      FOR i := 1 TO NumUnits DO IF s[i] = j THEN BEGIN
         k := k+1; ind[k] := i;
         IF c_list THEN BEGIN
            NewLin(1,0); writeln(lst,i:5,' ',name[i]);
         END;
         s_p[j] := s_p[j] + [i];
      END;
   END;
   c := s; cn := sn;
   IF c_list THEN Print_Network(ind);
   compute := TRUE; s_total := total_error; compute := FALSE;
   IF c_list THEN BEGIN
      Print_ModMatrix; Print_ErrMatrix;
      IF selval = 1 THEN Print_AveMatrix;
      IF selval = 2 THEN Print_ValMatrix;
   END;
   IF NOT RandomCls THEN BEGIN
      gotoxy(20,16);
      write('P*                 = ',s_total:11:3); readln;
      NewLin(1,0); writeln(lst,'final error =', s_total:12:3);
   END ELSE IF all THEN BEGIN
      newlin(1,0);
      write(lst,cnt:6,':'); list(s_p); writeln(lst,'P = ',s_total:8:3);
   END;
   {*** writeln(lst,'  final error =', s_total:12:3);}
   IF s_total < top_ten[10].pb THEN BEGIN
      WITH top_ten[0] DO BEGIN ib:=cnt; pb:=s_total; sb:=s_p END;
      x := top_ten[10]; u := 10;
      WHILE top_ten[u-1].pb > s_total DO BEGIN
     top_ten[u] := top_ten[u-1];
         IF (top_ten[u].pb < bignum) AND RandomCls THEN BEGIN
            GotoXY(55,6+u); write(top_ten[u].ib:6,top_ten[u].pb:12:3);
         END;
         u := u - 1
      END;
      top_ten[u] := top_ten[0];
      IF RandomCls THEN BEGIN
         GotoXY(55,6+u); write(top_ten[u].ib:6,top_ten[u].pb:12:3);
      END;
   END;
   IF s_total < 0.000001 THEN BEGIN
      IF NOT found THEN BEGIN
         found := TRUE; newlin(-15,0); newlin(3,0); writeln(lst);
         writeln(lst,' E Q U I V A L E N C E S'); writeln(lst);
      END;
      newlin(1,0); write(lst,cnt:6,':'); list(s_p); writeln(lst);
   END;
   IF savelev = 2 THEN save := s_total <= pbest + 0.000001;
   IF s_total < pbest THEN BEGIN
      pbest := s_total; ibest := cnt;
      IF RandomCls THEN BEGIN
         gotoxy(41,15); write(cnt:7); gotoxy(41,16); write(s_total:11:3)
      END;
   END;
END;

PROCEDURE PrintWeights(i:integer);
BEGIN
   GoToXY(15,8+i); ClrEol; write(i:2, w[i]:9:3, ' ':2);
   write(tina[i]);
   IF i = 9 THEN write(' = ',diadensity:8:5, outdensity:9:5);
END;

PROCEDURE weights;
VAR i: integer;
BEGIN
   ClrScr; GoToXY(15,1);  HighVideo;
   Writeln('***  STRAN / MODEL 2 - blockmodeling by local optimization ');
   GoToXY(15,2); write('W E I G H T S'); NormVideo;
   select := ReadInt(15,4,
           'error type:  1 - constant,  2 - size  = ',1,2 );
   limit  := ReadInt(15,5,
           'minimal dom/fun/reg size  = ',1,NumUnits );
{   selval := ReadInt(15,6,
           'values: 0 - without,  1 - average,  2 - median = ',0,2 ); }
   selval := ReadInt(15,6,
           'values: 0 - without,  1 - average  = ',0,1 );
   IF select = 2 THEN writeln(lst,'error type = size')
      ELSE  writeln(lst,'error type = constant');
{   penality := ReadInt(15,7,
           'penality                  = ',1,999);    }
   writeln(lst,'minimal dom/fun/par size  = ',limit:5 );
   writeln(lst,'averaging rule no/ave/med = ',selval:5 );
{   writeln(lst,'penality                  = ',penality:5 );  }
   writeln(lst); newlin(4,0);
   diadensity := 0.75; outdensity := 0.5; density := 1;
   w[0]  := 1; til[0]  := ' - '; tina[0]  := 'null';
   w[1]  := 1; til[1]  := 'com'; tina[1]  := 'complete';
   w[2]  := 1; til[2]  := 'rdo'; tina[2]  := 'row-dominant';
   w[3]  := 1; til[3]  := 'cdo'; tina[3]  := 'col-dominant';
   w[4]  := 1; til[4]  := 'reg'; tina[4]  := 'regular';
   w[5]  := 1; til[5]  := 'rre'; tina[5]  := 'row-regular';
   w[6]  := 1; til[6]  := 'cre'; tina[6]  := 'col-regular';
   w[7]  := 1; til[7]  := 'rfn'; tina[7]  := 'row-function';
   w[8]  := 1; til[8]  := 'cfn'; tina[8]  := 'col-function';
   w[9]  := 1; til[9]  := 'den'; tina[9]  := 'density';
   w[10] := 1; til[10] := 'dnc'; tina[10] := 'do not care';
   w[11] := 1; til[11] := 'one'; tina[11] := 'non-null';
   w[12] := 1; til[12] := 'sym'; tina[12] := 'symmetric';
   max_t := 12;
   FOR i := 0 TO max_t DO BEGIN prior[i] := i+1; seq[i+1] := i END;
   FOR i := 0 TO max_t DO PrintWeights(i);
   GoToXY(15,21); write('Priorities: ',prior[0]:3);
   FOR i := 1 TO max_t DO write(',',prior[i]:3);
   REPEAT
      i := ReadInt(15,22,
           '-1 - stop, 13 - priority, i - change i-th weight  = ',-1,max_t+1 );
      IF i >= 0 THEN IF i <= max_t THEN BEGIN
            w[i] := ReadReal(15,23, 'new weight = ', 0, 999999. );
            IF i = 9 THEN BEGIN
               diadensity :=
                  ReadReal(15,23, 'dia-density = ', 0.00001, 1.0 );
               outdensity :=
                  ReadReal(15,23, 'out-density = ', 0.00001, 1.0 );
            END;
            GoToXY(15,23); ClrEol;
         PrintWeights(i);
      END ELSE BEGIN
         REPEAT
            GoToXY(27,21); ClrEol;
            FOR i := 0 TO max_t DO BEGIN
               prior[i] := ReadInt(15,23,'priority  = ', 0,max_t+1 );
               GoToXY(27+3*i,21); write(prior[i]:2,',');
            END;
            FOR i := 0 TO max_t+1 DO seq[i] := -1; ok := TRUE;
            FOR i := 0 TO max_t DO IF prior[i] = 0 THEN active[i] := FALSE
               ELSE IF seq[prior[i]] = -1 THEN
                  BEGIN seq[prior[i]] := i END
               ELSE BEGIN seq[prior[i]] := 99; ok := FALSE END;
            FOR i := 0 TO max_t DO PrintWeights(i);
            GoToXY(15,21); write('Priorities: ',prior[0]:3);
            FOR i := 1 TO max_t DO write(',',prior[i]:3);
         UNTIL ok;
      END;
   UNTIL i < 0;
   newlin(3,0); writeln(lst); writeln(lst,'Weights, Priorities, Sequence'); writeln(lst);
   FOR i := 0 TO max_t DO BEGIN
      newlin(1,0); writeln(lst,i:5,' : ',w[i]:9:3, prior[i]:5, seq[i+1]:5);
   END;
   writeln(lst,'density = ',diadensity:8:5, outdensity:9:5);
   writeln(lst);
END;

BEGIN
   r_max := 1; t_print := 1;
   prgtit := 'M O D E L 2g'; ok := TRUE;
   Stran( 'blockmodeling by local optimization','2 Apr 00','2000');
   Get_Env(TRUE,er); asklev := 1;
   head := 'LST - list file'; St_file( lst, app, f_lst, head, ok );
   head := 'NET - network data file'; St_file( net, res, f_dat, head, ok );
   head := 'MDL - model data file'; St_file( mdl, res, f_rel, head, ok );
   ReadAna := TRUE;
   head := 'ANA - analysis file'; St_file( ana, res, f_ana, head, ReadAna );
   ClrScr; debug := 1; all := FALSE; test := FALSE;
   method := 'MODEL 2'; memode := '          ';
   prgnam := 'MODEL 2   '; linmax := 60; pagcnt := 0;
   newlin(100,0);
   newlin(1,0); writeln(lst,'NET - file = ', f_dat);
   newlin(1,0); writeln(lst,'LST - file = ', f_lst);
   newlin(1,0); writeln(lst,'ANA - file = ', f_ana);
   newlin(1,0); writeln(lst,'MDL - file = ', f_rel);
   newlin(1,0); writeln(lst); flush(lst);

{   Read Network Heading }

   IF eof(net) THEN BEGIN
      ClrScr; writeln('Empty NET file'); readln; halt(1);
   END;
   REPEAT read(net,ch); write(ch) UNTIL ch = ' ';
   REPEAT read(net,ch); write(ch) UNTIL ch = ' ';
   writeln; readln(net, NumUnits, r_max);
   IF NumUnits >= numax THEN Error_cl(5);
   readln(net, c_tit); netit := c_tit;
   writeln(c_tit); delay(500);
   Read_Names;
   FOR i := 1 TO NumUnits DO ind[i] := i;
   
   writeln('Memory before = ', Memavail );
   leng := sqr(NumUnits)+1;
   GetMem(tab,2*leng);
   GetMem(a,2*leng);
   GetMem(rel,leng);
   writeln('Memory after = ', Memavail );
   readln;

   found := TRUE; { lin := 3; rf := 3; rd := 0;}
   { diagonal := 0; general := 0; normal := 0; }
   { Read relation heading }
   ch := '.'; WHILE NOT (Eof(net) OR (ch='*')) DO read(net,ch);
   found := ch = '*'; { choice := -2;} writeln;
   IF found THEN BEGIN
      WHILE NOT (Eof(net) OR (ch= ' ')) DO
         BEGIN read(net,ch); write(ch) END;
      WHILE NOT (Eof(net) OR (ch<> ' ')) DO
         BEGIN read(net,ch); write(ch) END;
      WHILE NOT (Eof(net) OR (ch= ' ')) DO
         BEGIN read(net,ch); write(ch) END;
      Delay(1000);
      IF NOT Eof(net) THEN readln(net,form,sim);
      IF NOT Eof(net) THEN readln(net,title[1]);
      IF debug > 0 THEN writeln(title[1]);
      found := NOT eof(net);
   END;
   Read_Network;
   r_max := 1;
   Print_Network(ind);

   pen := 0;
   IF Eof(mdl) THEN BEGIN
      clm := clmax;
      FOR i := 1 TO clmax DO FOR j := 1 TO clmax DO BEGIN
         block_tip[i,j] := [0,1,4]; block_pen[i,j] := 1;
      END;
   END ELSE BEGIN
      readln(mdl); readln(mdl,clm); ch := ' ';
      IF (clm < 2) OR (clm > clmax) THEN clm := clmax;
      FOR i := 1 TO clm DO FOR j := 1 TO clm DO BEGIN
         block_tip[i,j] := [0]; block_pen[i,j] := 1;
      END;
      writeln(lst); writeln(lst,'Model matrix definition');
      WHILE NOT Eof(mdl) AND (ch=' ') DO BEGIN
         read(mdl,ch);
         IF ch=' ' THEN BEGIN
            read(mdl,p,q,penal); tipij := [];
            write(lst,p:5,q:5,penal:7); flush(lst);
            WHILE NOT eoln(mdl) DO BEGIN
               read(mdl,k); write(lst,k:3); tipij := tipij + [k]
            END;
            readln(mdl); writeln(lst);
            IF p=0 THEN
               CASE q OF
                  0 : FOR i := 1 TO clm DO BEGIN 
                         block_tip[i,i] := tipij; block_pen[i,i] := penal;
                      END;
                  1 : FOR i := 1 TO clm-1 DO FOR j := i+1 TO clm DO BEGIN
                         block_tip[i,j] := tipij; block_pen[i,j] := penal;
                      END;
                  2 : FOR i := 2 TO clm DO FOR j := 1 TO i-1 DO BEGIN
                         block_tip[i,j] := tipij; block_pen[i,j] := penal;
                      END;
                  3 : FOR i := 1 TO clm DO FOR j := 1 TO clm DO BEGIN
                         block_tip[i,j] := tipij; block_pen[i,j] := penal;
          END;
               END
            ELSE BEGIN
               block_tip[p,q] := tipij; block_pen[p,q] := penal;
            END;
         END;
      END;
      writeln(lst);
      IF NOT Eof(mdl) THEN BEGIN
         read(mdl,ch);
         IF (ch='c') OR (ch='C') THEN BEGIN
            readln(mdl); ch := ' ';
            writeln(lst); writeln(lst,'Constraints');
            WHILE NOT Eof(mdl) AND (ch=' ') DO BEGIN
               read(mdl,ch);
               IF ch=' ' THEN BEGIN
                  pen := pen+1;
                  IF pen <= pem THEN BEGIN
                     readln(mdl,pet[pen],penality[pen],pel[pen],pec[pen]);
                     writeln(lst,pen:5,pet[pen]:5,penality[pen]:5,pel[pen]:5,pec[pen]:3);
                  END ELSE BEGIN
                     ch := '*';
                     writeln(lst,'To many constraints (>',pem:1,')');
                  END;
                  flush(lst);
               END;
            END;
         END;
      END;
   END;
   IF test THEN BEGIN
      writeln(lst);
      writeln(lst,'Block model matrix');
      writeln(lst);
      FOR i := 1 TO clm DO BEGIN
         FOR j := 1 TO clm DO BEGIN
            tipij := block_tip[i,j]; k := 0;
            FOR p := 0 TO tmax DO
               IF p in tipij THEN k := 2*k + 1 ELSE k := 2*k;
            write(lst,k:6)
         END;
         writeln(lst);
      END;
      writeln(lst); writeln(lst);

      writeln(lst);
      writeln(lst,'Block penalty matrix');
      writeln(lst);
      FOR i := 1 TO clm DO BEGIN
         FOR j := 1 TO clm DO BEGIN
            write(lst,block_pen[i,j]:5)
         END;
         writeln(lst);
      END;
      writeln(lst); writeln(lst);
   END;

   ClrScr; HighVideo;
   Writeln('***  STRAN / MODEL 2 - blockmodeling by local optimization ');
   Write('Title = '); NormVideo; Writeln(c_tit); Writeln; HighVideo;
   Writeln('New title or <ENTER> : '); NormVideo; Readln(t);
   IF Length(t) > 0 THEN c_tit := t;
   IF Length(c_tit) > 70 THEN c_tit := Copy( c_tit, 1,70 );
   title[1] := c_tit;
   tit := c_tit;
   weights;
   cnt := 0; cnt_max := 0; c_max := 0; Randomize;
   keyclu := 'CLUST'; ask := FALSE; compute := FALSE; save := FALSE;

   WHILE TRUE DO BEGIN
      IF ReadANA THEN Read_Clustering;
      cnt := cnt + 1; RandomCls := NOT ReadANA;
      IF cnt > cnt_max THEN BEGIN
         ClrScr; HighVideo;
         gotoxy(15,7);  write('***  S T R A N  -  MODEL 2');
         gotoxy(15,8);  write('blockmodeling by local optimization');
         NormVideo;
         debug := ReadInt(15,12,
            '0 - stop,  1 - short list,  2 - list all  = ',0,4 );
         IF debug = 0 THEN GOTO 14;
         all := debug > 1; test := debug > 3;
{         transpose := 1 = ReadInt(15,14,
            '0 - relocate,  1 - relocate & transpose   = ',0,1); }
         transpose := TRUE;
         IF indana = 0 THEN c_new := ReadInt(15,14,
            'Number of clusters (0 - read clustering)  = ',0,clm)
            ELSE c_new := -1;
         IF c_new > 0 THEN Random_Initial;
      END;
      FOR i := 1 TO NumUnits DO cn[i] := 0;
      IF c_new <= 0 THEN Read_Initial ELSE Generate_Initial;
      c_list := c_read OR all OR test;
      IF RandomCls THEN BEGIN gotoxy(41,14); write(cnt:7); END;
      k := 0;
      FOR j := 1 TO c_max DO BEGIN
         IF c_list THEN BEGIN
            NewLin(1,0); writeln(lst,'CLUSTER ',j:1,' / ',cn[j]:1);
         END;
         m := 0;
         FOR i := 1 TO NumUnits DO IF c[i] = j THEN BEGIN
            k := k+1; ind[k] := i;
            IF c_list THEN BEGIN
               m := m+1; write(lst,i:5,' ',name[i]);
               IF m MOD 5 = 0 THEN
                  BEGIN NewLin(1,0); writeln(lst) END;
            END;
         END;
         IF c_list AND (m MOD 5 <> 0) THEN
            BEGIN NewLin(1,0); writeln(lst) END;
      END;
      s := c; s_max := c_max; sn := cn;
      IF c_list THEN Print_Network(ind);
      Local_Optimization;
      IF ask THEN BEGIN
         GoToXY(20,21); ClrEol;
         write( 'save clustering (Yes/No) ? ' ); readln( ch );
         save := (ch = 'y') OR (ch = 'Y');
      END;
      IF save  THEN Save_Clustering;
   END;
14:IF RandomCls THEN BEGIN
      newlin(-5,0); newlin(5,0); writeln(lst);
      writeln(lst,' NUMBER OF UNITS       = ', NumUnits:6 );
      writeln(lst,' NUMBER OF CLUSTERS    = ', c_max:6 );
      writeln(lst,' NUMBER OF CLUSTERINGS = ', cnt-1:6 );
      printresult;
   END;
   ClrScr; GotoXY(20,15); Write('M O D E L  2 - end'); GoToXY(1,24);
   Put_Env; Close(env); Close(net); Close(lst); Close(ana);
   Delay(1000); ClrScr;
END.


