
{
STRAN-PC - servisne rutine
Copyright (c)1990 Vladimir Batagelj
All rights reserved

CLUSE-PC/CL_UT0 - servisne rutine
(c)1988 Vladimir Batagelj

version *:
version 1:  11. jun 1989  Print_Env --> Change_Env
}

Unit STENV;

(***********************************************)
Interface
(***********************************************)

Uses Crt, Dos;

{
Stran/PC control environment
}

{$I control.pas}

   CONST
   {$IFDEF Demo}
      numax  = 12;   { max # of units }
      clmax  = 12;   { max # of clusters  = numax/2 }
      rmax   =  3;   { max # of relations }
      clver  = '  -  d e m o'; { version }
   {$ENDIF}
   {$IFDEF Study}
      numax  = 30;    { max # of units }
      clmax  = 30;    { max # of clusters  = numax/2 }
      rmax   =  5;    { max # of relations }
      clver  = '  -  student'; { version }
   {$ENDIF}
   {$IFDEF Work}
      numax  = 100;   { max # of units }
      clmax  = 50;    { max # of clusters  = numax/2 }
      rmax   =  1;
      disdim = 5050;  { distance matrix dimension  = numax*(numax+1)/2 }
      clver  = '  -  research'; { version }
   {$ENDIF}
   camax  = 9;                  { max # of cases / matrices }
   list   = 1000;
   bignum = 1E30;               { infinity }

TYPE
   Path_Name = STRING[40];
   file_action = ( res, rew, app );
   heading = String[40];
VAR

{ *Cluse/PC control environment variables }

   c_tit  : STRING[70]; 	{ title }
   c_fil  : Path_Name;		{ default file name }
   c_head : STRING[16]; 	{ CLUSE header }
   a_dat  : STRING;		{ date of last access }
   a_tim  : STRING;     	{ time of last access }
   nunits : LongInt;	        { number of units }
   nuvars ,			{ number of variables }
   ncases ,                     { number of cases / matrices }
   distyp ,			{ type of dissimilarity }
   asstyp ,                     { type of association coefficient }
   statyp ,                     { standardization type }
   source ,			{ 0-RAW/1-DAT }
   linmax ,			{ lines/page }
   width  ,			{ page width }
   skip   ,			{ 0-compact/1-spacing }
   debug  ,			{ debug listing }
   lstlev ,                     { list level }
   asklev : word;		{ ask level }
   f_env  ,     		{ ENV file }
   f_inf  ,      		{ INF file }
   f_raw  ,     		{ RAW file }
   f_dat  ,     		{ DAT file }
   f_nam  ,     		{ NAM file }
   f_var  ,     		{ VAR file }
   f_lst  ,		        { LST file }
   f_dbg  ,     		{ DBG file }
   f_ana  ,     		{ ANA file }
   f_rel  : Path_Name;          { REL file }
   f_dis                	{ DIS files }
                   : ARRAY [1..CaMax] of Path_Name;

{ CLUSE input/output environment }

   cl_ini : STRING[10];
   c_dat  : STRING;		{ current date }
   c_tim  : STRING;     	{ current time }
   dmin   , dmax   : real;	{ scale }
   pagcnt ,			{ page count }
   lincnt ,			{ line count }
   lleng  : integer;		{ line length }
   method ,
   memode ,
   prgtit : STRING;		{ program title }
   err_msg : ARRAY [ 1..15 ] OF STRING[40];
   inf    ,                     { CLUSE info file }
   env    ,                     { environment file }
   lst	  ,			{ output file }
   dbg	  : text;		{ debug file  }

   PROCEDURE WCenter(Y:integer;s:STRING);

   PROCEDURE Stran( des_name, version, c_r_y: STRING );

   FUNCTION ReadInt(X,Y:integer;prompt:STRING;a,b:integer):integer;

   FUNCTION ReadLong(X,Y:integer;prompt:STRING;a,b:LongInt):LongInt;

   FUNCTION ReadReal(X,Y:integer;prompt:STRING;a,b:real):real;

   PROCEDURE date( VAR datum, time : STRING );

   PROCEDURE ClrLine( X, Y: integer );

   PROCEDURE Error_cl( err: integer );

   PROCEDURE Select_File( head: heading; VAR name: Path_Name );

   PROCEDURE St_file( VAR f: text; rerwap: file_action;
      VAR fnam: Path_Name; head: heading; VAR ok: Boolean );

   PROCEDURE Ini_Env;

   PROCEDURE Get_Env( ask: Boolean; VAR error: integer );

   PROCEDURE Put_Env;

   PROCEDURE newlin( num, kind: integer );

   FUNCTION locecl( i, j, n: integer ): LongInt;

(***********************************************)
Implementation
(***********************************************)

   PROCEDURE WCenter(Y:integer;s:STRING);
   BEGIN
      GoToXY(1,Y); ClrEol; GoToXY(41-Length(s) DIV 2,Y); write(s);
   END;

   PROCEDURE Stran( des_name, version, c_r_y: STRING );
      VAR any_key: char;
   BEGIN
      ClrScr; HighVideo; WCenter(12,'S  T  R  A  N   /  PC'+clver);
      NormVideo; WCenter(13,'network structure analysis programs' ); HighVideo;
      WCenter(15,prgtit+'   - '+des_name ); NormVideo;
    {  WCenter(17,'version: '+version+'   (c)'+c_r_y+' Vladimir Batagelj'); }
      WCenter(16,'version: '+version);
      WCenter(18,'Copyright (c)'+c_r_y+' Vladimir Batagelj, All rights reserved');
      any_key := ReadKey;
   END;

   FUNCTION ReadInt(X,Y:integer;prompt:STRING;a,b:integer):integer;
      VAR answer:STRING; er,iv:integer; v:real; test,exit:Boolean;
   BEGIN
      test := a<b;
      IF test THEN BEGIN
         GoToXY(X,Y+1); ClrEol; write('[',a:5,',',b:5,']');
      END;
      REPEAT
         GoToXY(X,Y); ClrEol; write(prompt); readln(answer);
         Val(answer,v,er); exit := (Abs(v)<maxint) AND (er=0);
         IF exit THEN BEGIN
            iv := Round(v);
            IF test THEN exit := (a <= iv) AND (iv <= b);
         END;
      UNTIL exit;
      GoToXY(X,Y+1); ClrEol; ReadInt := iv;
   END {ReadInt};

   FUNCTION ReadLong(X,Y:integer;prompt:STRING;a,b:LongInt):LongInt;
      VAR answer:STRING; er:integer; iv:LongInt; v:real; test,exit:Boolean;
   BEGIN
      test := a<b;
      IF test THEN BEGIN
         GoToXY(X,Y+1); ClrEol; write('[',a:5,',',b:5,']');
      END;
      REPEAT
         GoToXY(X,Y); ClrEol; write(prompt); readln(answer);
         Val(answer,v,er); exit := (er=0);
         IF exit THEN BEGIN
            iv := Round(v);
            IF test THEN exit := (a <= iv) AND (iv <= b);
         END;
      UNTIL exit;
      GoToXY(X,Y+1); ClrEol; ReadLong := iv;
   END {ReadLong};

   FUNCTION ReadReal(X,Y:integer;prompt:STRING;a,b:real):real;
      VAR answer:STRING; er:integer; v:real; test,exit:Boolean;
   BEGIN
      test := a<b;
      IF test THEN BEGIN
         GoToXY(X,Y+1); ClrEol; write('[',a:9:4,',',b:9:4,']');
      END;
      REPEAT
         GoToXY(X,Y); ClrEol; write(prompt); readln(answer);
         Val(answer,v,er); exit := er=0;
         IF exit THEN IF test THEN exit := (a <= v) AND (v <= b);
      UNTIL exit;
      GoToXY(X,Y+1); ClrEol; ReadReal := v;
   END {ReadReal};

   PROCEDURE date( VAR datum, time : STRING );
   TYPE
      char3   = PACKED ARRAY [ 1 .. 3 ] OF char;
   VAR
      months  : ARRAY [ 1 .. 13 ] OF char3;	    { names of months }
      c_hour, c_min, c_sec, c_csec : word;	       { current time }
      c_day, c_month, c_year : word;    	       { current date }
      c_dayOfWeek : word;
      day, year, h, m, s : STRING;
   BEGIN
      GetDate( c_year, c_month, c_day, c_dayOfWeek );
      GetTime( c_hour, c_min, c_sec, c_csec );
      months[ 1] := 'Jan'; months[ 2] := 'Feb'; months[ 3] := 'Mar';
      months[ 4] := 'Apr'; months[ 5] := 'May'; months[ 6] := 'Jun';
      months[ 7] := 'Jul'; months[ 8] := 'Aug'; months[ 9] := 'Sep';
      months[10] := 'Oct'; months[11] := 'Nov'; months[12] := 'Dec';
      Str( c_day:2, day ); Str( c_year:4, year );
      datum := Concat( months[c_month],'-',day,'-',year);
      Str( c_hour:2, h ); Str( c_min:2, m ); Str( c_sec:2, s );
      time  := Concat( h,':',m,':',s );
   END { date };

   PROCEDURE ClrLine( X, Y: integer );
   BEGIN
      GoToXY(1,Y); ClrEol; GoToXY(X,Y);
   END {ClrLine};

   PROCEDURE Error_cl( err: integer );
   BEGIN
      ClrScr; GoToXY(20,12); HighVideo;
      write( 'CLUSE - fatal error in ', prgtit ); NormVideo;
      GoToXY( 10,14); write( err_msg[err] );
      readln; Halt(1);
   END { Error_cl };

   PROCEDURE Select_File( head: heading; VAR name: Path_Name );
      VAR answ : STRING; ok : Boolean;
   BEGIN
      ClrScr; GoToXY(20,10); write('SELECT FILE' );
      Gotoxy(20,12); HighVideo; write(head); LowVideo;
      GoToXY(20,13); write('Press RETURN to confirm current value' );
      NormVideo;
      REPEAT
         ok := TRUE;
         ClrLine(20,15); write( 'Current file = ', name );
         ClrLine(20,16); write( 'New file     = ');
         readln( answ );
         IF answ <> '' THEN BEGIN name := answ; ok := FALSE END;
      UNTIL ok;
   END { File_select };

   PROCEDURE St_file( VAR f: text; rerwap: file_action;
      VAR fnam: Path_Name; head: heading; VAR ok: Boolean );
   VAR io_status: integer;
   BEGIN
      IF asklev > 0 THEN Select_File( head, fnam );
      Assign( f, fnam );
      IF rerwap = rew THEN Rewrite( f )
      ELSE BEGIN
{$i-}
         IF rerwap = res THEN Reset( f ) ELSE Append( f );
         io_status := IOResult;
{$i+}
         IF io_status <> 0 THEN IF ok THEN BEGIN
{$i-}
            Rewrite( f ); io_status := IOResult;
{$i+}
            IF io_status = 0 THEN BEGIN
               Close( f );
               IF rerwap = res THEN Reset( f ) ELSE Append( f );
               io_status := IOResult;
            END ELSE Error_Cl(1);
         END ELSE Error_Cl(1);
         ok := io_status = 0;
      END;
   END;

   PROCEDURE Ini_Env;
      VAR i: Word; st: String[4];
   BEGIN
   { *Cluse/PC control environment variables }
      c_tit  :=
         'STRAN/PC - (c)1990 V.Batagelj ';
      c_head := '*Cluse/PC ';
      GetDir( 0, c_fil );      { default path }
      c_fil  := c_fil+'\TEST';       { default file name }
      a_dat  := 'Dec-04-1990'; { date of last access }
      a_tim  := '12:00:00';    { time of last access }
      nunits := 0;  { number of units }
      nuvars := 0;  { number of variables }
      ncases := 1;  { number of cases/matrices }
      distyp := 0;  { type of dissimilarity }
      asstyp := 0;  { type of assoc. coeff. }
      statyp := 0;  { standardization type }
      source := 0;  { 0-RAW/1-DAT }
      linmax := 56; { lines/page }
      width  := 80; { page width }
      skip   := 0;  { 0-compact/1-spacing }
      debug  := 1;  { debug listing }
      asklev := 2;  { ask level }
      lstlev := 1;  { list level }
      f_raw  := Concat(c_fil,'.RAW');  { RAW file }
      f_dat  := Concat(c_fil,'.NET');  { DAT file }
      f_nam  := Concat(c_fil,'.NAM');  { NAM file }
      f_var  := Concat(c_fil,'.VAR');  { VAR file }
      f_lst  := Concat(c_fil,'.LST');  { LST file }
      f_dbg  := Concat(c_fil,'.DBG');  { DBG file }
      f_ana  := Concat(c_fil,'.ANA');  { ANA file }
      f_rel  := Concat(c_fil,'.REL');  { REL file }
      f_env  := Concat(c_fil,'.ENV');  { ENV file }
      FOR i := 1 TO camax DO BEGIN                    { DIS files }
         f_dis[i] := Concat(c_fil,Chr(Ord('0')+i),'.DIS');
      END;
   { CLUSE input/output environment }

      date(c_dat,c_tim);  { current date and time }
      dmin   := 0;
      dmax   := 1;        { scale }
      pagcnt := 0;        { page count }
      lincnt := 100;      { line count }
      lleng  := 70;       { line length }
      method := '      ';
      memode := '      ';
      prgtit := '      '; { program title }

      err_msg[ 1 ] := 'bad or unexistent file';
      err_msg[ 2 ] := 'error 2';
      err_msg[ 3 ] := 'bad CLUSE.INF file';
      err_msg[ 4 ] := 'bad STRAN environement file';
      Str(numax:4,st);
      err_msg[ 5 ] := 'STRAN'+clver+': # of units <= '+st;
      err_msg[ 6 ] := 'unsupported form of input data';
      err_msg[ 7 ] := 'error 7';
      err_msg[ 8 ] := 'error 8';
      err_msg[ 9 ] := 'error 9';
   END ;

   PROCEDURE Get_Env( ask: Boolean; VAR error: integer );
   { get CLUSE environement }
   VAR
      io_status: integer;
      head: heading;
      opened: Boolean;
      kontra: STRING[16];
      i: Word;
   BEGIN
      error := 0; opened := TRUE; Ini_Env;
      Assign( inf, 'CLUSE.INF' );
{$i-}
      Reset(inf); io_status := IOResult;
{$i+}
      IF io_status <> 0 THEN BEGIN
         Rewrite( inf );
         GetDir( 0, f_env ); f_env := f_env+'\TEST.ENV';
         head := 'STRAN environement file';
         Select_File( head, f_env );
         Assign( env, f_env );
{$i-}
         Reset( env ); io_status := IOResult;
{$i+}
         IF io_status <> 0 THEN BEGIN
{$i-}
            Rewrite( env ); io_status := IOResult;
{$i+}
            IF io_status <>0 THEN Error_Cl(1);
            Reset( env );
         END ELSE Error_Cl(1);
         writeln( inf, f_env );
         Close( inf ); ask := FALSE;
         opened := NOT Eof(env);
      END ELSE BEGIN
         IF Eof(inf) THEN Error_cl( 3 ) ELSE readln( inf, f_env );
      END;
      IF opened THEN BEGIN
         IF ask THEN BEGIN
            head := 'CLUSE environement file';
            Select_File( head, f_env );
            Assign( env, f_env );
{$i-}
            Reset( env ); io_status := IOResult;
{$i+}
            IF (io_status <> 0) OR Eof(env) THEN BEGIN
               Rewrite( env ); Reset( env );
            END;
         END ELSE BEGIN
            Assign( env, f_env );
{$i-}
            Reset( env ); io_status := IOResult;
{$i+}
            IF (io_status <> 0) OR Eof(env) THEN BEGIN
               Rewrite( env ); Reset( env );
            END;
         END;
         IF NOT Eof(env) THEN BEGIN
            readln( env, kontra );
            IF kontra <> c_head THEN Error_Cl(4);
         { *Cluse/IBM PC control environment variables }
            readln( env, c_tit );
            readln( env, a_dat ); readln( env, a_tim );
            readln( env, nunits, nuvars, ncases, distyp, asstyp, statyp,
               source, linmax, width, skip, debug, lstlev, asklev);
            readln( env, f_raw );
            readln( env, f_dat );
            readln( env, f_nam );
            readln( env, f_var );
            readln( env, f_rel );
            readln( env, f_lst );
            readln( env, f_dbg );
            readln( env, f_ana );
            FOR i := 1 TO camax DO readln( env, f_dis[i] );
         END;
      END;
   END { Get_Env };

   PROCEDURE Put_Env;
 { put CLUSE environement }
      VAR i: Word;
   BEGIN
      Rewrite( env );
      a_dat := c_dat; a_tim := c_tim;
    { *Cluse/PC control environment variables }
      writeln( env, c_head );
      writeln( env, c_tit );
      writeln( env, a_dat ); writeln( env, a_tim );
      writeln( env, nunits:7, nuvars:5, ncases:3, distyp:3, asstyp:3,
         statyp:3, source:2, linmax:5, width:5, skip:2, debug:2,
         lstlev:2, asklev:2 );
      writeln( env, f_raw );
      writeln( env, f_dat );
      writeln( env, f_nam );
      writeln( env, f_var );
      writeln( env, f_rel );
      writeln( env, f_lst );
      writeln( env, f_dbg );
      writeln( env, f_ana );
      FOR i := 1 TO camax DO writeln( env, f_dis[i] );
      writeln( env, '*EOD' );
      Rewrite( inf ); writeln( inf, f_env ); Close( inf );
   END { Put_Env };

   PROCEDURE newlin( num, kind: integer );
   {
     program package           C L U S E  / PC
     program unit              list control
     programming language      pascal
     programmed by             Vladimir Batagelj
     date                      feb 1981 / DEC 10 - fortran 77
     date                      mar 1987 / ATARI ST
     date                      apr 1988 / pascal

   Subroutine  NEWLCL  controls the format of the CLUSE-results
   listings. The parameter  num  is the intended number of lines
   to be printed on the current page in the following output
   action(s). If there is not enough space, the output is started
   on the new page.
   Each page of listing has a heading at the top. The content of
   the heading and the page format is defined by the contents of
   variables from the print-environment  which have to be
   initialized before use.
   }
      VAR numl, tl: integer; tit: STRING[60];
   BEGIN
      numl := num;
      IF numl < 0 THEN
         IF (lincnt - numl) > linmax THEN numl := 100 ELSE numl := 0;
      lincnt := lincnt + numl;
      IF lincnt > linmax THEN BEGIN
         pagcnt := pagcnt + 1; IF numl > linmax THEN numl := 0;
         lincnt := numl; tl := Length( c_tit );
         IF width > 120 THEN BEGIN
            write( lst, Chr(12) ); IF tl > 70 THEN tl := 0 ELSE tl := 70-tl;
            writeln( lst, '***  S T R A N  -  Network analysis  - ',
            prgtit:30, method:13, memode:12, ' ':17,
               'Page -', pagcnt:10 );
            write( lst, c_tit:70 );
            IF kind > 0 THEN write( lst, ' ':5, 'scale ( ',
               dmin:8:3, ' ,', dmax:9:3, ' )', ' ':5 )
               ELSE write( lst, ' ':41 );
            writeln( lst, 'Date: ', c_dat: 11 );
            writeln( lst ); writeln( lst );
         END ELSE BEGIN
            writeln( lst,Chr(12),'***  S T R A N  -  Network analysis  - ',
               Copy(method,1,7):7, memode:12, ' ':10,
               'Page -', pagcnt:5 );
            IF kind > 0 THEN write( lst, Copy(c_tit,1,20):20,
               '   scale   ( ', dmin:8:3, ' ,', dmax:9:3, ' )', ' ':5 )
            ELSE BEGIN
               IF tl <= 60 THEN BEGIN tit := c_tit; tl := 62-tl END
               ELSE BEGIN tit := Copy(c_tit,1,60); tl := 2 END;
               write( lst, tit, ' ':tl );
            END;
            writeln( lst, '      ', c_dat: 11 );
            writeln( lst ); writeln( lst );
         END;
      END;
   END { newlin };

   FUNCTION locecl( i, j, n: integer ): LongInt;
   BEGIN
      IF i < j THEN locecl := ( (2*n - i)*(i-1) ) DIV 2  +  j
         ELSE locecl := ( (2*n - j)*(j-1) ) DIV 2  +  i
   END { locecl };

END.