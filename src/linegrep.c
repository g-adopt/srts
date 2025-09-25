//-------------------------------------------------------------
//
//  linegrep -
//  grep's rows of data from the input file
//
//  if INPUT is:
//
//  a1 a2 a3 a4 a5 a6
//  b1 b2 b3 b4 b5 b6
//  c1 c2 c3 c4 c5 c6
//  ...
//  y1 y2 y3 y4 y5 y6
//  z1 z2 z3 z4 z5 z6
//
//  prompt> cat INPUT | linegrep -l 15 17
//
//  returns:
//  o1 o2 o3 o4 o5 o6
//  p1 p2 p3 p4 p5 p6
//  q1 q2 q3 q4 q5 q6
//
//  Jeroen Ritsema, May 2001, Caltech
//-------------------------------------------------------------
//  PK obtained from JR, 2013

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void usage();

main( argc, argv )
int   argc;
char       *argv[];
{
    int       linenumber = 0, index = 0, line1, line2, idum;
    char      ss[180];


   if (argc < 4) usage(-1);

   while ( ++index < argc && argv[index][0] == '-' ) {
        switch ( argv[index][1] ) {
            case 'l':
                if ( sscanf( argv[++index], "%d", &line1 ) != 1 ||
                     sscanf( argv[++index], "%d", &line2 ) != 1) usage(-1);
                break;
            default:
                usage(-1);
        }
   }
   if (line2 < line1 ) {
      idum = line1;
      line1 = line2;
      line2 = idum;
   }

   linenumber = 1;
   while ( fgets(ss,180,stdin) != NULL ) {
      if (linenumber >= line1 && linenumber <= line2) {
         fputs(ss,stdout);
      }
      linenumber++;
   }

exit( 0 );
}

void    usage( exitstatus )
int     exitstatus;
{
   fprintf(stderr,"Usage: linegrep -l line1 line2]\n"); 
   exit( exitstatus );
}

