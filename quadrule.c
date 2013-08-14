# include <stdlib.h>
# include <stdio.h>
# include <string.h>
# include <math.h>
# include <time.h>

# include "quadrule.h"

/******************************************************************************/

void bashforth_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    BASHFORTH_SET sets abscissas and weights for Adams-Bashforth quadrature.
  
  Discussion:
  
    Adams-Bashforth quadrature formulas are normally used in solving
    ordinary differential equations, and are not really suitable for
    general quadrature computations.  However, an Adams-Bashforth formula
    is equivalent to approximating the integral of F(Y(X)) between X(M)
    and X(M+1), using an explicit formula that relies only on known values
    of F(Y(X)) at X(M-N+1) through X(M).  For this reason, the formulas
    have been included here.
  
    Suppose the unknown function is denoted by Y(X), with derivative
    F(Y(X)), and that approximate values of the function are known at a
    series of X values, which we write as X(1), X(2), ..., X(M).  We write
    the value Y(X(1)) as Y(1) and so on.
  
    Then the solution of the ODE Y'=F(X,Y) at the next point X(M+1) is
    computed by:
  
      Y(M+1-1] = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dX
             = Y(M) + H * Sum ( 1 <= I <= N ) W(I) * F(Y(M+1-I))
               approximately.
  
    In the documentation that follows, we replace F(Y(X)) by F(X).
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w(x-1] = 1.0;
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) F(X) dX.
  
    The quadrature formula:
  
      Sum ( 1 <= I <= ORDER ) weight[I) * F ( 1 - I ),
  
    The Adams-Bashforth formulas require equally spaced data.
  
    Here is how the formula is applied in the case with non-unit spacing:
  
      Integral ( A <= X <= A+H ) F(X) dX =
      H * Sum ( 1 <= I <= ORDER ) weight[I) * F ( A - (I-1)*H ),
      approximately.
  
    The reference lists the second coefficient of the order 8 Adams-Bashforth
    formula as
      weight[2-1] =  -1162169.0 / 120960.0
    but this should be
      weight[2-1] =  -1152169.0 / 120960.0
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 April 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Vladimir Krylov,
    Approximate Calculation of Integrals,
    Dover, 2006,
    ISBN: 0486445798,
    LC: QA311.K713.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3,
    LC: QA47.M315.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.  ORDER should be
    between 1 and 10, or 12, 14, 16, 18 or 20.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    weight[0] is the weight at X = 0, weight[1] the weight at X = -1,
    and so on.  The weights are rational, and should sum to 1.  Some
    weights may be negative.
*/
{
  double d;
  int i;

  if ( order == 1 )
  {
    weight[1-1] = 1.0;
  }
  else if ( order == 2 )
  {
    d = 2.0;

    weight[1-1] =   3.0 / d;
    weight[2-1] = - 1.0 / d;
  }
  else if ( order == 3 )
  {
    d = 12.0;

    weight[1-1] =   23.0 / d;
    weight[2-1] = - 16.0 / d;
    weight[3-1] =    5.0 / d;
  }
  else if ( order == 4 )
  {
    d = 24.0;

    weight[1-1] =   55.0 / d;
    weight[2-1] = - 59.0 / d;
    weight[3-1] =   37.0 / d;
    weight[4-1] =  - 9.0 / d;
  }
  else if ( order == 5 )
  {
    d = 720.0;

    weight[1-1] =   1901.0 / d;
    weight[2-1] = - 2774.0 / d;
    weight[3-1] =   2616.0 / d;
    weight[4-1] = - 1274.0 / d;
    weight[5-1] =    251.0 / d;
  }
  else if ( order == 6 )
  {
    d = 1440.0;

    weight[1-1] =   4277.0 / d;
    weight[2-1] = - 7923.0 / d;
    weight[3-1] =   9982.0 / d;
    weight[4-1] = - 7298.0 / d;
    weight[5-1] =   2877.0 / d;
    weight[6-1] =  - 475.0 / d;
  }
  else if ( order == 7 )
  {
    d = 60480.0;

    weight[1-1] =    198721.0 / d;
    weight[2-1] =  - 447288.0 / d;
    weight[3-1] =    705549.0 / d;
    weight[4-1] =  - 688256.0 / d;
    weight[5-1] =    407139.0 / d;
    weight[6-1] =  - 134472.0 / d;
    weight[7-1] =     19087.0 / d;
  }
  else if ( order == 8 )
  {
    d = 120960.0;

    weight[1-1] =     434241.0 / d;
    weight[2-1] =  - 1152169.0 / d;
    weight[3-1] =    2183877.0 / d;
    weight[4-1] =  - 2664477.0 / d;
    weight[5-1] =    2102243.0 / d;
    weight[6-1] =  - 1041723.0 / d;
    weight[7-1] =     295767.0 / d;
    weight[8-1] =    - 36799.0 / d;
  }
  else if ( order == 9 )
  {
    d = 3628800.0;

    weight[1-1] =   14097247.0 / d;
    weight[2-1] =  -43125206.0 / d;
    weight[3-1] =   95476786.0 / d;
    weight[4-1] = -139855262.0 / d;
    weight[5-1] =  137968480.0 / d;
    weight[6-1] =  -91172642.0 / d;
    weight[7-1] =   38833486.0 / d;
    weight[8-1] =   -9664106.0 / d;
    weight[9-1] =    1070017.0 / d;
  }
  else if ( order == 10 )
  {
    d = 7257600.0;

    weight[ 1-1] =   30277247.0 / d;
    weight[ 2-1] = -104995189.0 / d;
    weight[ 3-1] =  265932680.0 / d;
    weight[ 4-1] = -454661776.0 / d;
    weight[ 5-1] =  538363838.0 / d;
    weight[ 6-1] = -444772162.0 / d;
    weight[ 7-1] =  252618224.0 / d;
    weight[ 8-1] =  -94307320.0 / d;
    weight[ 9-1] =   20884811.0 / d;
    weight[10-1] =   -2082753.0 / d;
  }
  else if ( order == 12 )
  {
    d = 958003200.0;

    weight[ 1-1] =    4527766399.0 / d;
    weight[ 2-1] =  -19433810163.0 / d;
    weight[ 3-1] =   61633227185.0 / d;
    weight[ 4-1] = -135579356757.0 / d;
    weight[ 5-1] =  214139355366.0 / d;
    weight[ 6-1] = -247741639374.0 / d;
    weight[ 7-1] =  211103573298.0 / d;
    weight[ 8-1] = -131365867290.0 / d;
    weight[ 9-1] =   58189107627.0 / d;
    weight[10-1] =  -17410248271.0 / d;
    weight[11-1] =    3158642445.0 / d;
    weight[12-1] =    -262747265.0 / d;
  }
  else if ( order == 14 )
  {
    d = 5230697472000.0;

    weight[ 1-1] =    27511554976875.0 / d;
    weight[ 2-1] =  -140970750679621.0 / d;
    weight[ 3-1] =   537247052515662.0 / d;
    weight[ 4-1] = -1445313351681906.0 / d;
    weight[ 5-1] =  2854429571790805.0 / d;
    weight[ 6-1] = -4246767353305755.0 / d;
    weight[ 7-1] =  4825671323488452.0 / d;
    weight[ 8-1] = -4204551925534524.0 / d;
    weight[ 9-1] =  2793869602879077.0 / d;
    weight[10-1] = -1393306307155755.0 / d;
    weight[11-1] =   505586141196430.0 / d;
    weight[12-1] =  -126174972681906.0 / d;
    weight[13-1] =    19382853593787.0 / d;
    weight[14-1] =    -1382741929621.0 / d;
  }
  else if ( order == 16 )
  {
    d = 62768369664000.0;

    weight[ 1-1] =     362555126427073.0 / d;
    weight[ 2-1] =   -2161567671248849.0 / d;
    weight[ 3-1] =    9622096909515337.0 / d;
    weight[ 4-1] =  -30607373860520569.0 / d;
    weight[ 5-1] =   72558117072259733.0 / d;
    weight[ 6-1] = -131963191940828581.0 / d;
    weight[ 7-1] =  187463140112902893.0 / d;
    weight[ 8-1] = -210020588912321949.0 / d;
    weight[ 9-1] =  186087544263596643.0 / d;
    weight[10-1] = -129930094104237331.0 / d;
    weight[11-1] =   70724351582843483.0 / d;
    weight[12-1] =  -29417910911251819.0 / d;
    weight[13-1] =    9038571752734087.0 / d;
    weight[14-1] =   -1934443196892599.0 / d;
    weight[15-1] =     257650275915823.0 / d;
    weight[16-1] =     -16088129229375.0 / d;
  }
  else if ( order == 18 )
  {
    d = 64023737057280000.0;

    weight[ 1-1] =     401972381695456831.0 / d;
    weight[ 2-1] =   -2735437642844079789.0 / d;
    weight[ 3-1] =   13930159965811142228.0 / d;
    weight[ 4-1] =  -51150187791975812900.0 / d;
    weight[ 5-1] =  141500575026572531760.0 / d;
    weight[ 6-1] = -304188128232928718008.0 / d;
    weight[ 7-1] =  518600355541383671092.0 / d;
    weight[ 8-1] = -710171024091234303204.0 / d;
    weight[ 9-1] =  786600875277595877750.0 / d;
    weight[10-1] = -706174326992944287370.0 / d;
    weight[11-1] =  512538584122114046748.0 / d;
    weight[12-1] = -298477260353977522892.0 / d;
    weight[13-1] =  137563142659866897224.0 / d;
    weight[14-1] =  -49070094880794267600.0 / d;
    weight[15-1] =   13071639236569712860.0 / d;
    weight[16-1] =   -2448689255584545196.0 / d;
    weight[17-1] =     287848942064256339.0 / d;
    weight[18-1] =     -15980174332775873.0 / d;
  }
  else if ( order == 20 )
  {
    d = 102181884343418880000.0;

    weight[ 1-1] =      691668239157222107697.0 / d;
    weight[ 2-1] =    -5292843584961252933125.0 / d;
    weight[ 3-1] =    30349492858024727686755.0 / d;
    weight[ 4-1] =  -126346544855927856134295.0 / d;
    weight[ 5-1] =   399537307669842150996468.0 / d;
    weight[ 6-1] =  -991168450545135070835076.0 / d;
    weight[ 7-1] =  1971629028083798845750380.0 / d;
    weight[ 8-1] = -3191065388846318679544380.0 / d;
    weight[ 9-1] =  4241614331208149947151790.0 / d;
    weight[10-1] = -4654326468801478894406214.0 / d;
    weight[11-1] =  4222756879776354065593786.0 / d;
    weight[12-1] = -3161821089800186539248210.0 / d;
    weight[13-1] =  1943018818982002395655620.0 / d;
    weight[14-1] =  -970350191086531368649620.0 / d;
    weight[15-1] =   387739787034699092364924.0 / d;
    weight[16-1] =  -121059601023985433003532.0 / d;
    weight[17-1] =    28462032496476316665705.0 / d;
    weight[18-1] =    -4740335757093710713245.0 / d;
    weight[19-1] =      498669220956647866875.0 / d;
    weight[20-1] =      -24919383499187492303.0 / d;
  }
  else
  {
    fprintf ( stderr, "\nBASHFORTH_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 through 10, 12, 14, 16, 18 or 20.\n");
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = ( double ) ( -i );
  }

  return;
}
/******************************************************************************/

void bdf_set ( int order, double alpha[], double *beta, double *gamma )

/******************************************************************************/
/*
  Purpose:
  
    BDF_SET sets weights for backward differentiation ODE weights.
  
  Discussion:
  
    GAMMA * Y(N+1) = Sum ( 1 <= I <= ORDER ) ALPHA(I) * Y(N+1-I)
                     + dX * BETA * Y'(X(N+1),Y(N+1))
  
    This is equivalent to the backward differentiation corrector formulas.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 April 2006
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int ORDER, the order of the rule, between 1 and 6.
  
    Output, double ALPHA[ORDER], *BETA, *GAMMA, the weights.
*/
{
  if ( order == 1 )
  {
    *beta =      1.0;
    *gamma =     1.0;
    alpha[1-1] = 1.0;
  }
  else if ( order == 2 )
  {
    *beta =        2.0;
    *gamma =       3.0;
    alpha[1-1] =   4.0;
    alpha[2-1] = - 1.0;
  }
  else if ( order == 3 )
  {
    *beta =        6.0;
    *gamma =      11.0;
    alpha[1-1] =  18.0;
    alpha[2-1] = - 9.0;
    alpha[3-1] =   2.0;
  }
  else if ( order == 4 )
  {
    *beta =        12.0;
    *gamma =       25.0;
    alpha[1-1] =   48.0;
    alpha[2-1] = - 36.0;
    alpha[3-1] =   16.0;
    alpha[4-1] =  - 3.0;
  }
  else if ( order == 5 )
  {
    *beta =         60.0;
    *gamma =       137.0;
    alpha[1-1] =   300.0;
    alpha[2-1] = - 300.0;
    alpha[3-1] =   200.0;
    alpha[4-1] =  - 75.0;
    alpha[5-1] =    12.0;
  }
  else if ( order == 6 )
  {
    *beta =         60.0;
    *gamma =       147.0;
    alpha[1-1] =   360.0;
    alpha[2-1] = - 450.0;
    alpha[3-1] =   400.0;
    alpha[4-1] = - 225.0;
    alpha[5-1] =    72.0;
    alpha[6-1] =  - 10.0;
  }
  else
  {
    fprintf(stderr,"\n" );
    fprintf ( stderr, "BDF_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal order requested = %d\n", order);
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void bdfc_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    BDFC_SET sets weights for backward differentiation corrector quadrature.
  
  Discussion:
  
    A backward differentiation corrector formula is defined for a set
    of evenly spaced abscissas X(I) with X(1-1] = 1 and X(2-1] = 0.  Assuming
    that the values of the function to be integrated are known at the
    abscissas, the formula is written in terms of the function value at
    X(1), and the backward differences at X(1) that approximate the
    derivatives there.
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w[x-1] = 1.0;
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * BD**(I-1) F ( 1 ).
  
    Here, "BD**(I-1) F ( 1 )" denotes the (I-1)st backward difference
    of F at X = 1, using a spacing of 1.  In particular,
  
    BD**0 F(1-1] = F(1)
    BD**1 F(1-1] = F(1) - F(0)
    BD**2 F(1-1] = F(1) - 2 * F(0) + F(-1 )
  
    The relationship between a backward difference corrector and the
    corresponding Adams-Moulton formula may be illustrated for the
    BDF corrector of order 4:
  
      BD**0 F(1) - 1/2 * BD**1 F(1) - 1/12 * BD**2 F(1) - 1/24 * BDF**3 F(1)
      =            F(1)
        -  1/2 * ( F(1) -         F(0) )
        - 1/12 * ( F(1) - 2     * F(0) +        F(-1) )
        - 1/24 * ( F(1) - 3     * F(0) + 3    * F(-1) -        F(-2) )
      =   9/24 *   F(1) + 19/24 * F(0) - 5/24 * F(-1) + 1/24 * F(-2)
  
    which is the Adams-Moulton formula of order 4.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 April 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Simeon Fatunla,
    Numerical Methods for Initial Value Problems in Ordinary Differential
    Equations,
    Academic Press, 1988.
  
  Parameters:
  
    Input, int ORDER, the order of the rule, which can be
    any value from 1 to 19.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
# define ORDER_MAX 19

  int i;
  static double w[ORDER_MAX] = {
                  1.0,
                 -1.0 /               2.0,
                 -1.0 /              12.0,
                 -1.0 /              24.0,
                -19.0 /             720.0,
                 -3.0 /             160.0,
               -863.0 /           60480.0,
               -275.0 /           24792.0,
             -33953.0 /         3628800.0,
              -8183.0 /         1036800.0,
           -3250433.0 /       479001600.0,
              -4671.0 /          788480.0,
       -13695779093.0 /   2615348736000.0,
        -2224234463.0 /    475517952000.0,
      -132282840127.0 /  31384184832000.0,
        -2639651053.0 /    689762304000.0,
    111956703448001.0 /   3201186852864.0,
           50188465.0 /     15613165568.0,
   2334028946344463.0 / 786014494949376.0 };

  if ( ORDER_MAX < order )
  {
    fprintf(stderr,"\n" );
    fprintf ( stderr, "BDFC_SET - Fatal error!\n" );
    fprintf ( stderr, "  Input order ORDER = %d\n", order );
    fprintf ( stderr, "  exceeds maximum order ORDER_MAX = %d\n", ORDER_MAX);
    exit ( 1 );
  }
  for ( i = 0; i < order; i++ )
  {
    weight[i] = w[i];
  }

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = ( double ) ( 1 - i );
  }

  return;
#  undef ORDER_MAX
}
/******************************************************************************/

void bdfp_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    BDFP_SET sets weights for backward differentiation predictor quadrature.
  
  Discussion:
  
    A backward differentiation predictor formula is defined for a set
    of evenly spaced abscissas X(I) with X(1) = 1 and X(2) = 0.  Assuming
    that the values of the function to be integrated are known at the
    abscissas, the formula is written in terms of the function value at
    X(2), and the backward differences at X(2) that approximate the
    derivatives there.  A backward differentiation predictor formula
    is equivalent to an Adams-Bashforth formula of the same order.
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w(x) = 1.0;
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * BD**(I-1) F ( 0 ),
  
    Here, "BD**(I-1) F ( 0 )" denotes the (I-1)st backward difference
    of F at X = 0, using a spacing of 1.  In particular,
  
    BD**0 F(0) = F(0)
    BD**1 F(0) = F(0) - F(-1)
    BD**2 F(0) = F(0) - 2 * F(-1) + F(-2 )
  
    The relationship between a backward difference predictor and the
    corresponding Adams-Bashforth formula may be illustrated for the
    BDF predictor of order 3:
  
      BD**0 F(0) + 0.5 * BD**1 F(0) + 5/12 * BD**2 F(0)
      =            F(0)
        + 1/2  * ( F(0) -         F(1) )
        + 5/12 * ( F(0) - 2     * F(-1) +      F(-2) )
      =  23/12 *   F(0) - 16/12 * F(-1) + 5/12 F(-2)
  
    which is the Adams-Bashforth formula of order 3.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 April 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Simeon Fatunla,
    Numerical Methods for Initial Value Problems in Ordinary Differential
    Equations,
    Academic Press, 1988.
  
  Parameters:
  
    Input, int ORDER, the order of the rule, which can be
    any value from 1 to 19.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weight of the rule.
*/
{
# define ORDER_MAX 19

  int i;
  double w[ORDER_MAX] = {
                        1.0,
                        1.0 /                2.0,
                        5.0 /               12.0,
                        3.0 /                8.0,
                      251.0 /              720.0,
                       95.0 /              288.0,
                    19087.0 /            60480.0,
                     5257.0 /            17280.0,
                  1070017.0 /          3628800.0,
                    25713.0 /            89600.0,
                 26842253.0 /         95800320.0,
                  4777223.0 /         17418240.0,
             703604254357.0 /    2615348736000.0,
             106364763817.0 /     402361344000.0,
            1166309819657.0 /    4483454976000.0,
                 25221445.0 /         98402304.0,
         8092989203533249.0 /    3201186852864.0,
           85455477715379.0 /      34237292544.0,
     12600467236042756559.0 / 5109094217170944.0 };

  if ( ORDER_MAX < order )
  {
    fprintf(stderr,"\n" );
    fprintf ( stderr, "BDFP_SET - Fatal error!\n" );
    fprintf ( stderr, "  Input order ORDER = %d\n", order );
    fprintf ( stderr, "  exceeds maximum order ORDER_MAX = %d\n", ORDER_MAX);
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    weight[i] = w[i];
  }
  for ( i = 0; i < order; i++ )
  {
    xtab[i] = ( double ) ( - i );
  }

  return;
# undef ORDER_MAX
}
/******************************************************************************/

double bdf_sum ( double func ( double x ), int order, double xtab[],
  double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    BDF_SUM carries out explicit backward difference quadrature on [0,1].
  
  Discussion:
  
    The integral to approximate is:
  
      Integral ( 0 <= X <= 1 ) F(X) dX
  
    The quadrature formula is:
  
      RESULT = Sum ( 1 <= I <= ORDER ) WEIGHT(I) * BDF**(I-1) FUNC ( 0 )
  
    The integral from 0 to 1 is approximated using data at X = 0,
    -1, -2, ..., -ORDER+1.  This is a form of extrapolation, and
    the approximation can become poor as ORDER increases.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 April 2006
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, double FUNC ( double X ), the name of the function which
    evaluates the integrand.
  
    Input, int ORDER, the order of the rule.
  
    Input, double XTAB[ORDER], the abscissas of the rule.
  
    Input, double WEIGHT[ORDER], the weights of the rule.
  
    Output, double BDF_SUM, the approximate value of the integral.
*/
{
  int i;
  int j;
  double result;

  double diftab[order];

  for ( i = 0; i < order; i++ )
  {
    diftab[i] = func ( xtab[i] );
  }

  for ( i = 2; i <= order; i++ )
  {
    for ( j = i; j <= order; j++ )
    {
      diftab[order+i-j-1] = ( diftab[order+i-j-2] - diftab[order+i-j-1] );
    }
  }

  result = r8vec_dot ( order, weight, diftab );

  return result;
}
/******************************************************************************/

char ch_cap ( char ch )

/******************************************************************************/
/*
  Purpose:
  
    CH_CAP capitalizes a single character.
  
  Discussion:
  
    This routine should be equivalent to the library "toupper" function.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    19 July 1998
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, char CH, the character to capitalize.
  
    Output, char CH_CAP, the capitalized character.
*/
{
  if ( 97 <= ch && ch <= 122 )
  {
    ch = ch - 32;
  }

  return ch;
}
/******************************************************************************/

void cheb_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    CHEB_SET sets abscissas and weights for Chebyshev quadrature.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x-1] = 1.0;
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( xtab[I) )
  
    The Chebyshev rule is distinguished by using equal weights.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    30 April 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Hermann Engels,
    Numerical Quadrature and Cubature,
    Academic Press, 1980.
  
    Zdenek Kopal,
    Numerical Analysis,
    John Wiley, 1955.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER may only have the values 1, 2, 3, 4, 5, 6, 7 or 9.
    There are NO other Chebyshev rules with real abscissas.
  
    Output, double XTAB[ORDER], the abscissas of the rule,
    which are symmetric in [-1,1].
  
    Output, double WEIGHT[ORDER], the weights of the rule,
    which should each equal 2 / ORDER.
*/
{
  int i;

  if ( order == 1 )
  {
    xtab[1-1] = 0.0;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = - 1.0 / sqrt ( 3.0 );
    xtab[2-1] =   1.0 / sqrt ( 3.0 );
  }
  else if ( order == 3 )
  {
    xtab[1-1] = - 1.0 / sqrt ( 2.0 );
    xtab[2-1] =   0.0;
    xtab[3-1] =   1.0 / sqrt ( 2.0 );
  }
  else if ( order == 4 )
  {
    xtab[1-1] =   - sqrt ( ( 1.0 + 2.0/ sqrt ( 5.0 ) ) / 3.0 );
    xtab[2-1] =   - sqrt ( ( 1.0 - 2.0/ sqrt ( 5.0 ) ) / 3.0 );
    xtab[3-1] =     sqrt ( ( 1.0 - 2.0/ sqrt ( 5.0 ) ) / 3.0 );
    xtab[4-1] =     sqrt ( ( 1.0 + 2.0/ sqrt ( 5.0 ) ) / 3.0 );
  }
  else if ( order == 5 )
  {
    xtab[1-1] = - sqrt ( ( 5.0 + sqrt ( 11.0 ) ) / 12.0 );
    xtab[2-1] = - sqrt ( ( 5.0 - sqrt ( 11.0 ) ) / 12.0 );
    xtab[3-1] =   0.0;
    xtab[4-1] =   sqrt ( ( 5.0 - sqrt ( 11.0 ) ) / 12.0 );
    xtab[5-1] =   sqrt ( ( 5.0 + sqrt ( 11.0 ) ) / 12.0 );
  }
  else if ( order == 6 )
  {
    xtab[1-1] = - 0.866246818107820591383598;
    xtab[2-1] = - 0.422518653761111529118546;
    xtab[3-1] = - 0.266635401516704720331534;
    xtab[4-1] =   0.266635401516704720331534;
    xtab[5-1] =   0.422518653761111529118546;
    xtab[6-1] =   0.866246818107820591383598;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = - 0.883861700758049035704224;
    xtab[2-1] = - 0.529656775285156811385048;
    xtab[3-1] = - 0.323911810519907637519673;
    xtab[4-1] =   0.0;
    xtab[5-1] =   0.323911810519907637519673;
    xtab[6-1] =   0.529656775285156811385048;
    xtab[7-1] =   0.883861700758049035704224;
  }
  else if ( order == 9 )
  {
    xtab[1-1] = - 0.911589307728434473664949;
    xtab[2-1] = - 0.601018655380238071428128;
    xtab[3-1] = - 0.528761783057879993260181;
    xtab[4-1] = - 0.167906184214803943068031;
    xtab[5-1] =   0.0;
    xtab[6-1] =   0.167906184214803943068031;
    xtab[7-1] =   0.528761783057879993260181;
    xtab[8-1] =   0.601018655380238071428128;
    xtab[9-1] =   0.911589307728434473664949;
  }
  else
  {
    fprintf(stderr,"\n" );
    fprintf ( stderr, "CHEB_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 through 7, and 9.\n" );
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    weight[i] = 2.0 / ( double ) ( order );
  }

  return;
}
/******************************************************************************/

void chebyshev1_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 quadrature rule.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x) = 1.0 / sqrt ( 1 - x^2 ).
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) / sqrt ( 1 - x^2 ) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 February 2008
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be greater than 0.
  
    Output, double X[ORDER], the abscissas.
  
    Output, double W[ORDER], the weights.
*/
{
  int i;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    fprintf(stderr,"\n" );
    fprintf ( stderr, "CHEBYSHEV1_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order);
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    w[i] = pi / ( double ) ( order );
  }
  for ( i = 0; i < order; i++ )
  {
    x[i] = cos ( pi * ( double ) ( 2 * order - 1 - 2 * i )
                    / ( double ) ( 2 * order ) );
  }

  return;
}
/******************************************************************************/

double chebyshev1_integral ( int expon )

/******************************************************************************/
/*
  Purpose:
  
    CHEBYSHEV1_INTEGRAL evaluates a monomial Chebyshev type 1 integral.
  
  Discussion:
  
    To test a Chebyshev type 1 quadrature rule, we use it to approximate the
    integral of a monomial:
  
      integral ( -1 <= x <= +1 ) x^n / sqrt ( 1 - x^2 ) dx
  
    This routine is given the value of the exponent, and returns the
    exact value of the integral.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 February 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int EXPON, the exponent.
  
    Output, double CHEBYSHEV1_INTEGRAL, the value of the exact integral.
*/
{
  double bot;
  double exact;
  int i;
  double pi = 3.141592653589793;
  double top;
/*
  Get the exact value of the integral.
*/
  if ( ( expon % 2 ) == 0 )
  {
    top = 1;
    bot = 1;
    for ( i = 2; i <= expon; i = i + 2 )
    {
      top = top * ( i - 1 );
      bot = bot *   i;
    }
	
    exact = pi * ( double ) ( top ) / ( double ) ( bot );
  }
  else
  {
    exact = 0.0;	
  }

  return exact;
}
/******************************************************************************/

void chebyshev2_compute ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 quadrature rule.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x) = sqrt ( 1 - x^2 ).
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X)  sqrt ( 1 - x^2 )  dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 February 2008
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be greater than 0.
  
    Output, double X[ORDER], the abscissas.
  
    Output, double W[ORDER], the weights.
*/
{
  double angle;
  int i;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    fprintf(stderr,"\n" );
    fprintf ( stderr, "CHEBYSHEV2_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order);
    exit ( 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    angle = pi * ( double ) ( order - i ) / ( double ) ( order + 1 );
    w[i] = pi / ( double ) ( order + 1 ) * pow ( sin ( angle ), 2 );
    x[i] = cos ( angle );
  }

  return;
}
/******************************************************************************/

double chebyshev2_integral ( int expon )

/******************************************************************************/
/*
  Purpose:
  
    CHEBYSHEV2_INTEGRAL evaluates a monomial Chebyshev type 2 integral.
  
  Discussion:
  
    To test a Chebyshev type 2 quadrature rule, we use it to approximate the
    integral of a monomial:
  
      integral ( -1 <= x <= +1 ) x^n * sqrt ( 1 - x^2 ) dx
  
    This routine is given the value of the exponent, and returns the
    exact value of the integral.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 February 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int EXPON, the exponent.
  
    Output, double CHEBYSHEV2_INTEGRAL, the value of the exact integral.
*/
{
  double bot;
  double exact;
  int i;
  double pi = 3.141592653589793;
  double top;
/*
  Get the exact value of the integral.
*/
  if ( ( expon % 2 ) == 0 )
  {
    top = 1;
    bot = 1;
    for ( i = 2; i <= expon; i = i + 2 )
    {
      top = top * ( i - 1 );
      bot = bot *   i;
    }

	bot = bot * ( double ) ( expon + 2 );

    exact = pi * ( double ) ( top ) / ( double ) ( bot );
  }
  else
  {
    exact = 0.0;
  }
  return exact;
}
/******************************************************************************/

void chebyshev3_compute ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    CHEBYSHEV3_COMPUTE computes a Gauss-Chebyshev type 3 quadrature rule.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x) = 1 / sqrt ( 1 - x**2 )
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) f(x) / sqrt ( 1 - x**2 ) dx
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
    The ORDER = 1 rule is exceptional.  It consists of a single
    point at 0, with weight PI.
  
    For rules with ORDER = 2 or greater, the following remarks apply:
  
    If ORDER points are used, then Gauss-Chebyshev quadrature
    will compute the integral exactly, whenever F(X) is a polynomial
    of degree 2*ORDER-3 or less.
  
    The abscissas include -1 and 1.
  
    The first and last weights are 0.5 * PI / ( ORDER - 1),
    and all other weights are PI / ( ORDER - 1 ).
  
    If the order is doubled, the abscissas of the new rule include
    all the points of the old rule.  This fact can be used to
    efficiently implement error estimation.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    30 April 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule, which must be at least 1.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  double angle;
  int i;
  double pi = 3.141592653589793;

  if ( order < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "CHEBYSHEV3_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  ORDER must be at least 1.\n" );
    fprintf ( stderr, "  The input value was ORDER = %d\n", order );
    exit ( 1 );
  }
/*
  Take care of the special case ORDER = 1.
*/
  if ( order == 1 )
  {
    xtab[0] = 0.0;
    weight[0] = pi;
    return;
  }

  for ( i = 0; i < order; i++ )
  {
    angle = ( double ) ( order - 1 - i ) * pi / ( double ) ( order - 1 );
    xtab[i] = cos ( angle );
  }

  weight[0] = pi / ( double ) ( 2 * ( order - 1 ) );
  for ( i = 1; i < order-1; i++ )
  {
    weight[i] = pi / ( double ) ( order - 1 );
  }
  weight[order-1] = pi / ( double ) ( 2 * ( order - 1 ) );

  return;
}
/******************************************************************************/

void clenshaw_curtis_compute ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis quadrature rule.
  
  Discussion:
  
    This method uses a direct approach.  The paper by Waldvogel
    exhibits a more efficient approach using Fourier transforms.
  
    The integration interval is [ -1, 1 ].
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
  
    The abscissas for the rule of order ORDER can be regarded
    as the cosines of equally spaced angles between 180 and 0 degrees:
  
      X(I) = cos ( ( I - 1 ) * PI / ( ORDER - 1 ) )
  
    except for the basic case ORDER = 1, when
  
      X(1) = 0.
  
    A Clenshaw-Curtis rule that uses ORDER points will integrate
    exactly all polynomials of degrees 0 through ORDER-1.  If ORDER
    is odd, then by symmetry the polynomial of degree ORDER will
    also be integrated exactly.
  
    If the value of ORDER is increased in a sensible way, then
    the new set of abscissas will include the old ones.  One such
    sequence would be ORDER(K) = 2*K+1 for K = 0, 1, 2, ...
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    18 October 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Charles Clenshaw, Alan Curtis,
    A Method for Numerical Integration on an Automatic Computer,
    Numerische Mathematik,
    Volume 2, Number 1, December 1960, pages 197-205.
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  
  Parameters:
  
    Input, int N, the order of the rule.
  
    Output, double X[N], W[N], the abscissas and weights of the rule.
*/
{
  double b;
  int i;
  int j;
  static double pi = 3.141592653589793;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "CLENSHAW_CURTIS_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
    return;
  }

  double theta[n];

  for ( i = 1; i <= n; i++ )
  {
    theta[i-1] = ( double ) ( i - 1 ) * pi
               / ( double ) ( n - 1 );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( theta[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 1.0;

    for ( j = 1; j <= ( ( n - 1 ) / 2 ); j++ )
    {
      if ( 2 * j == ( n - 1 ) )
      {
        b = 1.0;
      }
      else
      {
        b = 2.0;
      }

      w[i] = w[i] - b * cos ( 2.0 * ( double ) ( j ) * theta[i] )
           / ( double ) ( 4 * j * j - 1 );
    }
  }

  w[0] = w[0] / ( double ) ( n - 1 );
  for ( i = 2; i <= n-1; i++ )
  {
    w[i-1] = 2.0 * w[i-1] / ( double ) ( n - 1 );
  }
  w[n-1] = w[n-1] / ( double ) ( n - 1 );

  return;
}
/******************************************************************************/

void clenshaw_curtis_set ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    CLENSHAW_CURTIS_SET sets a Clenshaw-Curtis quadrature rule.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
  
    The abscissas for the rule of order ORDER can be regarded
    as the cosines of equally spaced angles between 180 and 0 degrees:
  
      X(I) = cos ( ( I - 1 ) * PI / ( ORDER - 1 ) )
  
    except for the basic case ORDER = 1, when
  
      X(1) = 0.
  
    A Clenshaw-Curtis rule that uses ORDER points will integrate
    exactly all polynomials of degrees 0 through ORDER-1.  If ORDER
    is odd, then by symmetry the polynomial of degree ORDER will
    also be integrated exactly.
  
    If the value of ORDER is increased in a sensible way, then
    the new set of abscissas will include the old ones.  One such
    sequence would be ORDER(K) = 2*K+1 for K = 0, 1, 2, ...
    Thus, in the table below, the abscissas for order 9 include
    those for order 5.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    14 May 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Charles Clenshaw, Alan Curtis,
    A Method for Numerical Integration on an Automatic Computer,
    Numerische Mathematik,
    Volume 2, Number 1, December 1960, pages 197-205.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 17, 33, 65 or 129.
  
    Output, double X[ORDER], the abscissas of the rule.
  
    Output, double W[ORDER], the weights of the rule.
    The weights are symmetric and sum to 2.
*/
{
  if ( order == 1 )
  {
    x[1-1] =  0.00000000000000000000;
    w[1-1] =  2.00000000000000000000;
  }
  else if ( order == 2 )
  {
    x[1-1] = -1.00000000000000000000;
    x[2-1] =  1.00000000000000000000;

    w[1-1] =  1.00000000000000000000;
    w[2-1] =  1.00000000000000000000;
  }
  else if ( order == 3 )
  {
    x[1-1] = -1.00000000000000000000;
    x[2-1] =  0.00000000000000000000;
    x[3-1] =  1.00000000000000000000;

    w[1-1] =  0.33333333333333333333;
    w[2-1] =  1.33333333333333333333;
    w[3-1] =  0.33333333333333333333;
  }
  else if ( order == 4 )
  {
    x[1-1] = -1.00000000000000000000;
    x[2-1] = -0.50000000000000000000;
    x[3-1] =  0.50000000000000000000;
    x[4-1] =  1.00000000000000000000;

    w[1-1] =  0.11111111111111111111;
    w[2-1] =  0.88888888888888888889;
    w[3-1] =  0.88888888888888888889;
    w[4-1] =  0.11111111111111111111;
  }
  else if ( order == 5 )
  {
    x[1-1] = -1.00000000000000000000;
    x[2-1] = -0.70710678118654752440;
    x[3-1] =  0.00000000000000000000;
    x[4-1] =  0.70710678118654752440;
    x[5-1] =  1.00000000000000000000;

    w[1-1] =  0.06666666666666666667;
    w[2-1] =  0.53333333333333333333;
    w[3-1] =  0.80000000000000000000;
    w[4-1] =  0.53333333333333333333;
    w[5-1] =  0.06666666666666666667;
  }
  else if ( order == 6 )
  {
    x[1-1] = -1.00000000000000000000;
    x[2-1] = -0.80901699437494742410;
    x[3-1] = -0.30901699437494742410;
    x[4-1] =  0.30901699437494742410;
    x[5-1] =  0.80901699437493732410;
    x[6-1] =  1.00000000000000000000;

    w[1-1] =  0.04000000000000000000;
    w[2-1] =  0.36074304120001121619;
    w[3-1] =  0.59925695879998878381;
    w[4-1] =  0.59925695879998878381;
    w[5-1] =  0.36074304120001121619;
    w[6-1] =  0.04000000000000000000;
  }
  else if ( order == 7 )
  {
    x[1-1] = -1.00000000000000000000;
    x[2-1] = -0.86602540378443864676;
    x[3-1] = -0.50000000000000000000;
    x[4-1] =  0.00000000000000000000;
    x[5-1] =  0.50000000000000000000;
    x[6-1] =  0.86602540378443864676;
    x[7-1] =  1.00000000000000000000;

    w[1-1] =  0.02857142857142857143;
    w[2-1] =  0.25396825396825396825;
    w[3-1] =  0.45714285714285714286;
    w[4-1] =  0.52063492063492063492;
    w[5-1] =  0.45714285714285714286;
    w[6-1] =  0.25396825396825396825;
    w[7-1] =  0.02857142857142857143;
  }
  else if ( order == 8 )
  {
    x[1-1] = -1.00000000000000000000;
    x[2-1] = -0.90096886790241912624;
    x[3-1] = -0.62348980185873353053;
    x[4-1] = -0.22252093395631440429;
    x[5-1] =  0.22252093395631440429;
    x[6-1] =  0.62348980185873353053;
    x[7-1] =  0.90096886790241910624;
    x[8-1] =  1.00000000000000000000;

    w[1-1] =  0.02040816326530612245;
    w[2-1] =  0.19014100721820835178;
    w[3-1] =  0.35224242371815911533;
    w[4-1] =  0.43720840579832641044;
    w[5-1] =  0.43720840579832641044;
    w[6-1] =  0.35224242371815911533;
    w[7-1] =  0.19014100721820835178;
    w[8-1] =  0.02040816326530612245;
  }
  else if ( order == 9 )
  {
    x[1-1] = -1.00000000000000000000;
    x[2-1] = -0.92387953251128675613;
    x[3-1] = -0.70710678118654752440;
    x[4-1] = -0.38268343236508977173;
    x[5-1] =  0.00000000000000000000;
    x[6-1] =  0.38268343236508977173;
    x[7-1] =  0.70710678118654752440;
    x[8-1] =  0.92387953251128675613;
    x[9-1] =  1.00000000000000000000;

    w[1-1] =  0.01587301587301587302;
    w[2-1] =  0.14621864921601815501;
    w[3-1] =  0.27936507936507936508;
    w[4-1] =  0.36171785872048978150;
    w[5-1] =  0.39365079365079365079;
    w[6-1] =  0.36171785872048978150;
    w[7-1] =  0.27936507936507936508;
    w[8-1] =  0.14621864921601815501;
    w[9-1] =  0.01587301587301587302;
  }
  else if ( order == 10 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.93969262078590838405;
    x[3-1]  = -0.76604444311897903520;
    x[4-1]  = -0.50000000000000000000;
    x[5-1]  = -0.17364817766693034885;
    x[6-1]  =  0.17364817766693034885;
    x[7-1]  =  0.50000000000000000000;
    x[8-1]  =  0.76604444311897903520;
    x[9-1]  =  0.93969262078590838405;
    x[10-1] =  1.00000000000000000000;

    w[1-1]  =  0.01234567901234567901;
    w[2-1]  =  0.11656745657203712296;
    w[3-1]  =  0.22528432333810440813;
    w[4-1]  =  0.30194003527336860670;
    w[5-1]  =  0.34386250580414418320;
    w[6-1]  =  0.34386250580414418320;
    w[7-1]  =  0.30194003527336860670;
    w[8-1]  =  0.22528432333810440813;
    w[9-1]  =  0.11656745657203712296;
    w[10-1] =  0.01234567901234567901;
  }
  else if ( order == 11 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.95105651629515357212;
    x[3-1]  = -0.80901699437494742410;
    x[4-1]  = -0.58778525229247312917;
    x[5-1]  = -0.30901699437494742410;
    x[6-1]  =  0.00000000000000000000;
    x[7-1]  =  0.30901699437494742410;
    x[8-1]  =  0.58778525229247312917;
    x[9-1]  =  0.80901699437494742410;
    x[10-1] =  0.95105651629515357212;
    x[11-1] =  1.00000000000000000000;

    w[1-1]  =  0.01010101010101010101;
    w[2-1]  =  0.09457905488370156116;
    w[3-1]  =  0.18563521442424776529;
    w[4-1]  =  0.25358833328368660623;
    w[5-1]  =  0.29921327042423708320;
    w[6-1]  =  0.31376623376623376623;
    w[7-1]  =  0.29921327042423708320;
    w[8-1]  =  0.25358833328368660623;
    w[9-1]  =  0.18563521442424776529;
    w[10-1] =  0.09457905488370156116;
    w[11-1] =  0.01010101010101010101;
  }
  else if ( order == 12 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.95949297361449738989;
    x[3-1]  = -0.84125353283118116886;
    x[4-1]  = -0.65486073394528506406;
    x[5-1]  = -0.41541501300188642553;
    x[6-1]  = -0.14231483827328514044;
    x[7-1]  =  0.14231483827328514044;
    x[8-1]  =  0.41541501300188642553;
    x[9-1]  =  0.65486073394528506406;
    x[10-1] =  0.84125353283118116886;
    x[11-1] =  0.95949297361449738989;
    x[12-1] =  1.00000000000000000000;

    w[1-1]  =  0.00826446280991735537;
    w[2-1]  =  0.07856015374620000543;
    w[3-1]  =  0.15504045508256136552;
    w[4-1]  =  0.21556254600086858099;
    w[5-1]  =  0.25991734106691617602;
    w[6-1]  =  0.28265504129353651666;
    w[7-1]  =  0.28265504129353651666;
    w[8-1]  =  0.25991734106691617602;
    w[9-1]  =  0.21556254600086858099;
    w[10-1] =  0.15504045508256136552;
    w[11-1] =  0.07856015374620000543;
    w[12-1] =  0.00826446280991735537;
  }
  else if ( order == 13 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.96592582628906828675;
    x[3-1]  = -0.86602540378443864676;
    x[4-1]  = -0.70710678118654752440;
    x[5-1]  = -0.50000000000000000000;
    x[6-1]  = -0.25881904510252076235;
    x[7-1]  =  0.00000000000000000000;
    x[8-1]  =  0.25881904510252076235;
    x[9-1]  =  0.50000000000000000000;
    x[10-1] =  0.70710678118654752440;
    x[11-1] =  0.86602540378443864676;
    x[12-1] =  0.96592582628906828675;
    x[13-1] =  1.00000000000000000000;

    w[1-1]  =  0.00699300699300699301;
    w[2-1]  =  0.06605742495207439452;
    w[3-1]  =  0.13154253154253154253;
    w[4-1]  =  0.18476338476338476338;
    w[5-1]  =  0.22697302697302697303;
    w[6-1]  =  0.25267569378104433860;
    w[7-1]  =  0.26198986198986198986;
    w[8-1]  =  0.25267569378104433860;
    w[9-1]  =  0.22697302697302697303;
    w[10-1] =  0.18476338476338476338;
    w[11-1] =  0.13154253154253154253;
    w[12-1] =  0.06605742495207439452;
    w[13-1] =  0.00699300699300699301;
  }
  else if ( order == 14 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.97094181742605202716;
    x[3-1]  = -0.88545602565320989590;
    x[4-1]  = -0.74851074817110109863;
    x[5-1]  = -0.56806474673115580251;
    x[6-1]  = -0.35460488704253562597;
    x[7-1]  = -0.12053668025532305335;
    x[8-1]  =  0.12053668025532305335;
    x[9-1]  =  0.35460488704253562597;
    x[10-1] =  0.56806474673115580251;
    x[11-1] =  0.74851074817110109863;
    x[12-1] =  0.88545602565320989590;
    x[13-1] =  0.97094181742605202716;
    x[14-1] =  1.00000000000000000000;

    w[1-1]  =  0.00591715976331360947;
    w[2-1]  =  0.05646531376341444627;
    w[3-1]  =  0.11276867248985655881;
    w[4-1]  =  0.16003802611671868523;
    w[5-1]  =  0.19899241036578321848;
    w[6-1]  =  0.22590304977856444935;
    w[7-1]  =  0.23991536772234903239;
    w[8-1]  =  0.23991536772234903239;
    w[9-1]  =  0.22590304977856444935;
    w[10-1] =  0.19899241036578321848;
    w[11-1] =  0.16003802611671868523;
    w[12-1] =  0.11276867248985655881;
    w[13-1] =  0.05646531376341444627;
    w[14-1] =  0.00591715976331360947;
  }
  else if ( order == 15 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.97492791218182360702;
    x[3-1]  = -0.90096886790241912624;
    x[4-1]  = -0.78183148246802980871;
    x[5-1]  = -0.62348980185873353053;
    x[6-1]  = -0.43388373911755812048;
    x[7-1]  = -0.22252093395631440429;
    x[8-1]  =  0.00000000000000000000;
    x[9-1]  =  0.22252093395631440429;
    x[10-1] =  0.43388373911755812048;
    x[11-1] =  0.62348980185873353053;
    x[12-1] =  0.78183148246802980871;
    x[13-1] =  0.90096886790241912624;
    x[14-1] =  0.97492791218182360702;
    x[15-1] =  1.00000000000000000000;

    w[1-1]  =  0.00512820512820512821;
    w[2-1]  =  0.04869938729508823855;
    w[3-1]  =  0.09782039167605215913;
    w[4-1]  =  0.13966507849560431803;
    w[5-1]  =  0.17560578900106674677;
    w[6-1]  =  0.20205146748238357364;
    w[7-1]  =  0.21888151163057340180;
    w[8-1]  =  0.22429633858205286777;
    w[9-1]  =  0.21888151163057340180;
    w[10-1] =  0.20205146748238357364;
    w[11-1] =  0.17560578900106674677;
    w[12-1] =  0.13966507849560431803;
    w[13-1] =  0.09782039167605215913;
    w[14-1] =  0.04869938729508823855;
    w[15-1] =  0.00512820512820512821;
  }
  else if ( order == 16 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.97814760073380563793;
    x[3-1]  = -0.91354545764260089550;
    x[4-1]  = -0.80901699437494742410;
    x[5-1]  = -0.66913060635885821383;
    x[6-1]  = -0.50000000000000000000;
    x[7-1]  = -0.30901699437494742410;
    x[8-1]  = -0.10452846326765347140;
    x[9-1]  =  0.10452846326765347140;
    x[10-1] =  0.30901699437494742410;
    x[11-1] =  0.50000000000000000000;
    x[12-1] =  0.66913060635885821383;
    x[13-1] =  0.80901699437494742410;
    x[14-1] =  0.91354545764260089550;
    x[15-1] =  0.97814760073380563793;
    x[16-1] =  1.00000000000000000000;

    w[1-1]  =  0.00444444444444444444;
    w[2-1]  =  0.04251476624752508988;
    w[3-1]  =  0.08553884025933288291;
    w[4-1]  =  0.12294010082849361533;
    w[5-1]  =  0.15573317603967369176;
    w[6-1]  =  0.18132978132978132978;
    w[7-1]  =  0.19921478132638853955;
    w[8-1]  =  0.20828410952436040635;
    w[9-1]  =  0.20828410952436040635;
    w[10-1] =  0.19921478132638853955;
    w[11-1] =  0.18132978132978132978;
    w[12-1] =  0.15573317603967369176;
    w[13-1] =  0.12294010082849361533;
    w[14-1] =  0.08553884025933288291;
    w[15-1] =  0.04251476624752508988;
    w[16-1] =  0.00444444444444444444;
  }
  else if ( order == 17 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.98078528040323044913;
    x[3-1]  = -0.92387953251128675613;
    x[4-1]  = -0.83146961230254523708;
    x[5-1]  = -0.70710678118654752440;
    x[6-1]  = -0.55557023301960222474;
    x[7-1]  = -0.38268343236508977173;
    x[8-1]  = -0.19509032201612826785;
    x[9-1]  =  0.00000000000000000000;
    x[10-1] =  0.19509032201612826785;
    x[11-1] =  0.38268343236508977173;
    x[12-1] =  0.55557023301960222474;
    x[13-1] =  0.70710678118654752440;
    x[14-1] =  0.83146961230254523708;
    x[15-1] =  0.92387953251128675613;
    x[16-1] =  0.98078528040323044913;
    x[17-1] =  1.00000000000000000000;

    w[1-1]  =  0.00392156862745098039;
    w[2-1]  =  0.03736870283720561032;
    w[3-1]  =  0.07548233154315183441;
    w[4-1]  =  0.10890555258189093044;
    w[5-1]  =  0.13895646836823307412;
    w[6-1]  =  0.16317266428170330256;
    w[7-1]  =  0.18147378423649335700;
    w[8-1]  =  0.19251386461292564687;
    w[9-1]  =  0.19641012582189052777;
    w[10-1] =  0.19251386461292564687;
    w[11-1] =  0.18147378423649335700;
    w[12-1] =  0.16317266428170330256;
    w[13-1] =  0.13895646836823307412;
    w[14-1] =  0.10890555258189093044;
    w[15-1] =  0.07548233154315183441;
    w[16-1] =  0.03736870283720561032;
    w[17-1] =  0.00392156862745098039;
  }
  else if ( order == 33 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.99518472667219688624;
    x[3-1]  = -0.98078528040323044913;
    x[4-1]  = -0.95694033573220886494;
    x[5-1]  = -0.92387953251128675613;
    x[6-1]  = -0.88192126434835502971;
    x[7-1]  = -0.83146961230254523708;
    x[8-1]  = -0.77301045336273696081;
    x[9-1]  = -0.70710678118654752440;
    x[10-1] = -0.63439328416364549822;
    x[11-1] = -0.55557023301960222474;
    x[12-1] = -0.47139673682599764856;
    x[13-1] = -0.38268343236508977173;
    x[14-1] = -0.29028467725446236764;
    x[15-1] = -0.19509032201612826785;
    x[16-1] = -0.098017140329560601994;
    x[17-1] =  0.000000000000000000000;
    x[18-1] =  0.098017140329560601994;
    x[19-1] =  0.19509032201612826785;
    x[20-1] =  0.29028467725446236764;
    x[21-1] =  0.38268343236508977173;
    x[22-1] =  0.47139673682599764856;
    x[23-1] =  0.55557023301960222474;
    x[24-1] =  0.63439328416364549822;
    x[25-1] =  0.70710678118654752440;
    x[26-1] =  0.77301045336273696081;
    x[27-1] =  0.83146961230254523708;
    x[28-1] =  0.88192126434835502971;
    x[29-1] =  0.92387953251128675613;
    x[30-1] =  0.95694033573220886494;
    x[31-1] =  0.98078528040323044913;
    x[32-1] =  0.99518472667219688624;
    x[33-1] =  1.00000000000000000000;

    w[1-1]  =  0.00097751710654936461;
    w[2-1]  =  0.00939319796295501470;
    w[3-1]  =  0.01923424513268114918;
    w[4-1]  =  0.02845791667723369009;
    w[5-1]  =  0.03759434191404720602;
    w[6-1]  =  0.04626276283775174949;
    w[7-1]  =  0.05455501630398031044;
    w[8-1]  =  0.06227210954529400455;
    w[9-1]  =  0.06942757563043545090;
    w[10-1] =  0.07588380044138847048;
    w[11-1] =  0.08163481765493851023;
    w[12-1] =  0.08657753844182743544;
    w[13-1] =  0.09070611286772099874;
    w[14-1] =  0.09394324443876873573;
    w[15-1] =  0.09629232594548817919;
    w[16-1] =  0.09769818820805558182;
    w[17-1] =  0.09817857778176829677;
    w[18-1] =  0.09769818820805558182;
    w[19-1] =  0.09629232594548817919;
    w[20-1] =  0.09394324443876873573;
    w[21-1] =  0.09070611286772099874;
    w[22-1] =  0.08657753844182743544;
    w[23-1] =  0.08163481765493851023;
    w[24-1] =  0.07588380044138847048;
    w[25-1] =  0.06942757563043545090;
    w[26-1] =  0.06227210954529400455;
    w[27-1] =  0.05455501630398031044;
    w[28-1] =  0.04626276283775174949;
    w[29-1] =  0.03759434191404720602;
    w[30-1] =  0.02845791667723369009;
    w[31-1] =  0.01923424513268114918;
    w[32-1] =  0.00939319796295501470;
    w[33-1] =  0.00097751710654936461;
  }
  else if ( order == 65 )
  {
    x[1-1]  = -1.00000000000000000000;
    x[2-1]  = -0.99879545620517239271;
    x[3-1]  = -0.99518472667219688624;
    x[4-1]  = -0.98917650996478097345;
    x[5-1]  = -0.98078528040323044913;
    x[6-1]  = -0.97003125319454399260;
    x[7-1]  = -0.95694033573220886494;
    x[8-1]  = -0.94154406518302077841;
    x[9-1]  = -0.92387953251128675613;
    x[10-1] = -0.90398929312344333159;
    x[11-1] = -0.88192126434835502971;
    x[12-1] = -0.85772861000027206990;
    x[13-1] = -0.83146961230254523708;
    x[14-1] = -0.80320753148064490981;
    x[15-1] = -0.77301045336273696081;
    x[16-1] = -0.74095112535495909118;
    x[17-1] = -0.70710678118654752440;
    x[18-1] = -0.67155895484701840063;
    x[19-1] = -0.63439328416364549822;
    x[20-1] = -0.59569930449243334347;
    x[21-1] = -0.55557023301960222474;
    x[22-1] = -0.51410274419322172659;
    x[23-1] = -0.47139673682599764856;
    x[24-1] = -0.42755509343028209432;
    x[25-1] = -0.38268343236508977173;
    x[26-1] = -0.33688985339222005069;
    x[27-1] = -0.29028467725446236764;
    x[28-1] = -0.24298017990326388995;
    x[29-1] = -0.19509032201612826785;
    x[30-1] = -0.14673047445536175166;
    x[31-1] = -0.098017140329560601994;
    x[32-1] = -0.049067674327418014255;
    x[33-1] =  0.000000000000000000000;
    x[34-1] =  0.049067674327418014255;
    x[35-1] =  0.098017140329560601994;
    x[36-1] =  0.14673047445536175166;
    x[37-1] =  0.19509032201612826785;
    x[38-1] =  0.24298017990326388995;
    x[39-1] =  0.29028467725446236764;
    x[40-1] =  0.33688985339222005069;
    x[41-1] =  0.38268343236508977173;
    x[42-1] =  0.42755509343028209432;
    x[43-1] =  0.47139673682599764856;
    x[44-1] =  0.51410274419322172659;
    x[45-1] =  0.55557023301960222474;
    x[46-1] =  0.59569930449243334347;
    x[47-1] =  0.63439328416364549822;
    x[48-1] =  0.67155895484701840063;
    x[49-1] =  0.70710678118654752440;
    x[50-1] =  0.74095112535495909118;
    x[51-1] =  0.77301045336273696081;
    x[52-1] =  0.80320753148064490981;
    x[53-1] =  0.83146961230254523708;
    x[54-1] =  0.85772861000027206990;
    x[55-1] =  0.88192126434835502971;
    x[56-1] =  0.90398929312344333159;
    x[57-1] =  0.92387953251128675613;
    x[58-1] =  0.94154406518302077841;
    x[59-1] =  0.95694033573220886494;
    x[60-1] =  0.97003125319454399260;
    x[61-1] =  0.98078528040323044913;
    x[62-1] =  0.98917650996478097345;
    x[63-1] =  0.99518472667219688624;
    x[64-1] =  0.99879545620517239271;
    x[65-1] =  1.00000000000000000000;

    w[1-1]  =  0.00024420024420024420;
    w[2-1]  =  0.00235149067531170332;
    w[3-1]  =  0.00483146544879091264;
    w[4-1]  =  0.00719269316173611402;
    w[5-1]  =  0.00958233879528379039;
    w[6-1]  =  0.01192339471421277160;
    w[7-1]  =  0.01425206043235199679;
    w[8-1]  =  0.01653498765728958965;
    w[9-1]  =  0.01878652974179578354;
    w[10-1] =  0.02098627442973743378;
    w[11-1] =  0.02314069493435819848;
    w[12-1] =  0.02523506498175476590;
    w[13-1] =  0.02727225714146838686;
    w[14-1] =  0.02924065319746833770;
    w[15-1] =  0.03114129710406762447;
    w[16-1] =  0.03296454656997632997;
    w[17-1] =  0.03471049818092511427;
    w[18-1] =  0.03637092028663918309;
    w[19-1] =  0.03794545992128481711;
    w[20-1] =  0.03942698871295609976;
    w[21-1] =  0.04081501340035783384;
    w[22-1] =  0.04210333111141810203;
    w[23-1] =  0.04329151496169082935;
    w[24-1] =  0.04437417923925731580;
    w[25-1] =  0.04535110955166067221;
    w[26-1] =  0.04621766751092557684;
    w[27-1] =  0.04697395904661414870;
    w[28-1] =  0.04761604458525019296;
    w[29-1] =  0.04814443257251220341;
    w[30-1] =  0.04855584485714105274;
    w[31-1] =  0.04885125664306609371;
    w[32-1] =  0.04902801843102555294;
    w[33-1] =  0.04908762351494245585;
    w[34-1] =  0.04902801843102555294;
    w[35-1] =  0.04885125664306609371;
    w[36-1] =  0.04855584485714105274;
    w[37-1] =  0.04814443257251220341;
    w[38-1] =  0.04761604458525019296;
    w[39-1] =  0.04697395904661414870;
    w[40-1] =  0.04621766751092557684;
    w[41-1] =  0.04535110955166067221;
    w[42-1] =  0.04437417923925731580;
    w[43-1] =  0.04329151496169082935;
    w[44-1] =  0.04210333111141810203;
    w[45-1] =  0.04081501340035783384;
    w[46-1] =  0.03942698871295609976;
    w[47-1] =  0.03794545992128481711;
    w[48-1] =  0.03637092028663918309;
    w[49-1] =  0.03471049818092511427;
    w[50-1] =  0.03296454656997632997;
    w[51-1] =  0.03114129710406762447;
    w[52-1] =  0.02924065319746833770;
    w[53-1] =  0.02727225714146838686;
    w[54-1] =  0.02523506498175476590;
    w[55-1] =  0.02314069493435819848;
    w[56-1] =  0.02098627442973743378;
    w[57-1] =  0.01878652974179578354;
    w[58-1] =  0.01653498765728958965;
    w[59-1] =  0.01425206043235199679;
    w[60-1] =  0.01192339471421277160;
    w[61-1] =  0.00958233879528379039;
    w[62-1] =  0.00719269316173611402;
    w[63-1] =  0.00483146544879091264;
    w[64-1] =  0.00235149067531170332;
    w[65-1] =  0.00024420024420024420;
  }
  else if ( order == 129 )
  {
    x[1-1]   = -1.00000000000000000000;
    x[2-1]   = -0.99969881869620422012;
    x[3-1]   = -0.99879545620517239271;
    x[4-1]   = -0.99729045667869021614;
    x[5-1]   = -0.99518472667219688624;
    x[6-1]   = -0.99247953459870999816;
    x[7-1]   = -0.98917650996478097345;
    x[8-1]   = -0.98527764238894124477;
    x[9-1]   = -0.98078528040323044913;
    x[10-1]  = -0.97570213003852854446;
    x[11-1]  = -0.97003125319454399260;
    x[12-1]  = -0.96377606579543986669;
    x[13-1]  = -0.95694033573220886494;
    x[14-1]  = -0.94952818059303666720;
    x[15-1]  = -0.94154406518302077841;
    x[16-1]  = -0.93299279883473888771;
    x[17-1]  = -0.92387953251128675613;
    x[18-1]  = -0.91420975570353065464;
    x[19-1]  = -0.90398929312344333159;
    x[20-1]  = -0.89322430119551532034;
    x[21-1]  = -0.88192126434835502971;
    x[22-1]  = -0.87008699110871141865;
    x[23-1]  = -0.85772861000027206990;
    x[24-1]  = -0.84485356524970707326;
    x[25-1]  = -0.83146961230254523708;
    x[26-1]  = -0.81758481315158369650;
    x[27-1]  = -0.80320753148064490981;
    x[28-1]  = -0.78834642762660626201;
    x[29-1]  = -0.77301045336273696081;
    x[30-1]  = -0.75720884650648454758;
    x[31-1]  = -0.74095112535495909118;
    x[32-1]  = -0.72424708295146692094;
    x[33-1]  = -0.70710678118654752440;
    x[34-1]  = -0.68954054473706692462;
    x[35-1]  = -0.67155895484701840063;
    x[36-1]  = -0.65317284295377676408;
    x[37-1]  = -0.63439328416364549822;
    x[38-1]  = -0.61523159058062684548;
    x[39-1]  = -0.59569930449243334347;
    x[40-1]  = -0.57580819141784530075;
    x[41-1]  = -0.55557023301960222474;
    x[42-1]  = -0.53499761988709721066;
    x[43-1]  = -0.51410274419322172659;
    x[44-1]  = -0.49289819222978403687;
    x[45-1]  = -0.47139673682599764856;
    x[46-1]  = -0.44961132965460660005;
    x[47-1]  = -0.42755509343028209432;
    x[48-1]  = -0.40524131400498987091;
    x[49-1]  = -0.38268343236508977173;
    x[50-1]  = -0.35989503653498814878;
    x[51-1]  = -0.33688985339222005069;
    x[52-1]  = -0.31368174039889147666;
    x[53-1]  = -0.29028467725446236764;
    x[54-1]  = -0.26671275747489838633;
    x[55-1]  = -0.24298017990326388995;
    x[56-1]  = -0.21910124015686979723;
    x[57-1]  = -0.19509032201612826785;
    x[58-1]  = -0.17096188876030122636;
    x[59-1]  = -0.14673047445536175166;
    x[60-1]  = -0.12241067519921619850;
    x[61-1]  = -0.098017140329560601994;
    x[62-1]  = -0.073564563599667423529;
    x[63-1]  = -0.049067674327418014255;
    x[64-1]  = -0.024541228522912288032;
    x[65-1]  =  0.00000000000000000000;
    x[66-1]  =  0.024541228522912288032;
    x[67-1]  =  0.049067674327418014255;
    x[68-1]  =  0.073564563599667423529;
    x[69-1]  =  0.098017140329560601994;
    x[70-1]  =  0.12241067519921619850;
    x[71-1]  =  0.14673047445536175166;
    x[72-1]  =  0.17096188876030122636;
    x[73-1]  =  0.19509032201612826785;
    x[74-1]  =  0.21910124015686979723;
    x[75-1]  =  0.24298017990326388995;
    x[76-1]  =  0.26671275747489838633;
    x[77-1]  =  0.29028467725446236764;
    x[78-1]  =  0.31368174039889147666;
    x[79-1]  =  0.33688985339222005069;
    x[80-1]  =  0.35989503653498814878;
    x[81-1]  =  0.38268343236508977173;
    x[82-1]  =  0.40524131400498987091;
    x[83-1]  =  0.42755509343028209432;
    x[84-1]  =  0.44961132965460660005;
    x[85-1]  =  0.47139673682599764856;
    x[86-1]  =  0.49289819222978403687;
    x[87-1]  =  0.51410274419322172659;
    x[88-1]  =  0.53499761988709721066;
    x[89-1]  =  0.55557023301960222474;
    x[90-1]  =  0.57580819141784530075;
    x[91-1]  =  0.59569930449243334347;
    x[92-1]  =  0.61523159058062684548;
    x[93-1]  =  0.63439328416364549822;
    x[94-1]  =  0.65317284295377676408;
    x[95-1]  =  0.67155895484701840063;
    x[96-1]  =  0.68954054473706692462;
    x[97-1]  =  0.70710678118654752440;
    x[98-1]  =  0.72424708295146692094;
    x[99-1]  =  0.74095112535495909118;
    x[100-1] =  0.75720884650648454758;
    x[101-1] =  0.77301045336273696081;
    x[102-1] =  0.78834642762660626201;
    x[103-1] =  0.80320753148064490981;
    x[104-1] =  0.81758481315158369650;
    x[105-1] =  0.83146961230254523708;
    x[106-1] =  0.84485356524970707326;
    x[107-1] =  0.85772861000027206990;
    x[108-1] =  0.87008699110871141865;
    x[109-1] =  0.88192126434835502971;
    x[110-1] =  0.89322430119551532034;
    x[111-1] =  0.90398929312344333159;
    x[112-1] =  0.91420975570353065464;
    x[113-1] =  0.92387953251128675613;
    x[114-1] =  0.93299279883473888771;
    x[115-1] =  0.94154406518302077841;
    x[116-1] =  0.94952818059303666720;
    x[117-1] =  0.95694033573220886494;
    x[118-1] =  0.96377606579543986669;
    x[119-1] =  0.97003125319454399260;
    x[120-1] =  0.97570213003852854446;
    x[121-1] =  0.98078528040323044913;
    x[122-1] =  0.98527764238894124477;
    x[123-1] =  0.98917650996478097345;
    x[124-1] =  0.99247953459870999816;
    x[125-1] =  0.99518472667219688624;
    x[126-1] =  0.99729045667869021614;
    x[127-1] =  0.99879545620517239271;
    x[128-1] =  0.99969881869620422012;
    x[129-1] =  1.00000000000000000000;

    w[1-1]   =  0.00006103888176768602;
    w[2-1]   =  0.00058807215382869754;
    w[3-1]   =  0.00120930061875273991;
    w[4-1]   =  0.00180308126695362360;
    w[5-1]   =  0.00240715327877140915;
    w[6-1]   =  0.00300345869904497128;
    w[7-1]   =  0.00360197835812614147;
    w[8-1]   =  0.00419553798718534675;
    w[9-1]   =  0.00478862143341336763;
    w[10-1]  =  0.00537724746840184621;
    w[11-1]  =  0.00596388034730799521;
    w[12-1]  =  0.00654590843862298928;
    w[13-1]  =  0.00712483332325489785;
    w[14-1]  =  0.00769875778896082811;
    w[15-1]  =  0.00826865154203087108;
    w[16-1]  =  0.00883303867470133581;
    w[17-1]  =  0.00939256583934814871;
    w[18-1]  =  0.00994602784923457905;
    w[19-1]  =  0.01049386202576892125;
    w[20-1]  =  0.01103504877427254184;
    w[21-1]  =  0.01156988348290849967;
    w[22-1]  =  0.01209748052807164113;
    w[23-1]  =  0.01261803597977743271;
    w[24-1]  =  0.01313076516693974630;
    w[25-1]  =  0.01363579321293772047;
    w[26-1]  =  0.01413241437853094133;
    w[27-1]  =  0.01462070254634350205;
    w[28-1]  =  0.01510001572479266783;
    w[29-1]  =  0.01557039073899425960;
    w[30-1]  =  0.01603123858745057916;
    w[31-1]  =  0.01648256956220377909;
    w[32-1]  =  0.01692383985846499368;
    w[33-1]  =  0.01735504125411394958;
    w[34-1]  =  0.01777566938875279997;
    w[35-1]  =  0.01818570377926339481;
    w[36-1]  =  0.01858467519566908661;
    w[37-1]  =  0.01897255587067948426;
    w[38-1]  =  0.01934890842392451844;
    w[39-1]  =  0.01971370183700155725;
    w[40-1]  =  0.02006652805198357604;
    w[41-1]  =  0.02040735612003867863;
    w[42-1]  =  0.02073580533490147816;
    w[43-1]  =  0.02105184759002011131;
    w[44-1]  =  0.02135512797425970725;
    w[45-1]  =  0.02164562356712882440;
    w[46-1]  =  0.02192300400598756892;
    w[47-1]  =  0.02218725355897195088;
    w[48-1]  =  0.02243806539722630184;
    w[49-1]  =  0.02267543270456671718;
    w[50-1]  =  0.02289907134390605882;
    w[51-1]  =  0.02310898491627407168;
    w[52-1]  =  0.02330491126131143273;
    w[53-1]  =  0.02348686571193163505;
    w[54-1]  =  0.02365460746057766523;
    w[55-1]  =  0.02380816473024258975;
    w[56-1]  =  0.02394731750476901502;
    w[57-1]  =  0.02407210792327850000;
    w[58-1]  =  0.02418233623893147567;
    w[59-1]  =  0.02427805942075745923;
    w[60-1]  =  0.02435909748927643184;
    w[61-1]  =  0.02442552306156708690;
    w[62-1]  =  0.02447717542743444284;
    w[63-1]  =  0.02451414358881568292;
    w[64-1]  =  0.02453628559651495473;
    w[65-1]  =  0.02454370750551418263;
    w[66-1]  =  0.02453628559651495473;
    w[67-1]  =  0.02451414358881568292;
    w[68-1]  =  0.02447717542743444284;
    w[69-1]  =  0.02442552306156708690;
    w[70-1]  =  0.02435909748927643184;
    w[71-1]  =  0.02427805942075745923;
    w[72-1]  =  0.02418233623893147567;
    w[73-1]  =  0.02407210792327850000;
    w[74-1]  =  0.02394731750476901502;
    w[75-1]  =  0.02380816473024258975;
    w[76-1]  =  0.02365460746057766523;
    w[77-1]  =  0.02348686571193163505;
    w[78-1]  =  0.02330491126131143273;
    w[79-1]  =  0.02310898491627407168;
    w[80-1]  =  0.02289907134390605882;
    w[81-1]  =  0.02267543270456671718;
    w[82-1]  =  0.02243806539722630184;
    w[83-1]  =  0.02218725355897195088;
    w[84-1]  =  0.02192300400598756892;
    w[85-1]  =  0.02164562356712882440;
    w[86-1]  =  0.02135512797425970725;
    w[87-1]  =  0.02105184759002011131;
    w[88-1]  =  0.02073580533490147816;
    w[89-1]  =  0.02040735612003867863;
    w[90-1]  =  0.02006652805198357604;
    w[91-1]  =  0.01971370183700155725;
    w[92-1]  =  0.01934890842392451844;
    w[93-1]  =  0.01897255587067948426;
    w[94-1]  =  0.01858467519566908661;
    w[95-1]  =  0.01818570377926339481;
    w[96-1]  =  0.01777566938875279997;
    w[97-1]  =  0.01735504125411394958;
    w[98-1]  =  0.01692383985846499368;
    w[99-1]  =  0.01648256956220377909;
    w[100-1] =  0.01603123858745057916;
    w[101-1] =  0.01557039073899425960;
    w[102-1] =  0.01510001572479266783;
    w[103-1] =  0.01462070254634350205;
    w[104-1] =  0.01413241437853094133;
    w[105-1] =  0.01363579321293772047;
    w[106-1] =  0.01313076516693974630;
    w[107-1] =  0.01261803597977743271;
    w[108-1] =  0.01209748052807164113;
    w[109-1] =  0.01156988348290849967;
    w[110-1] =  0.01103504877427254184;
    w[111-1] =  0.01049386202576892125;
    w[112-1] =  0.00994602784923457905;
    w[113-1] =  0.00939256583934814871;
    w[114-1] =  0.00883303867470133581;
    w[115-1] =  0.00826865154203087108;
    w[116-1] =  0.00769875778896082811;
    w[117-1] =  0.00712483332325489785;
    w[118-1] =  0.00654590843862298928;
    w[119-1] =  0.00596388034730799521;
    w[120-1] =  0.00537724746840184621;
    w[121-1] =  0.00478862143341336763;
    w[122-1] =  0.00419553798718534675;
    w[123-1] =  0.00360197835812614147;
    w[124-1] =  0.00300345869904497128;
    w[125-1] =  0.00240715327877140915;
    w[126-1] =  0.00180308126695362360;
    w[127-1] =  0.00120930061875273991;
    w[128-1] =  0.00058807215382869754;
    w[129-1] =  0.00006103888176768602;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "CLENSHAW_CURTIS_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 to 17, 33, 65 or 129.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void fejer1_compute ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    FEJER1_COMPUTE computes a Fejer type 1 quadrature rule.
  
  Discussion:
  
    This method uses a direct approach.  The paper by Waldvogel
    exhibits a more efficient approach using Fourier transforms.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    05 March 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.
  
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  
  Parameters:
  
    Input, int N, the order of the rule.
  
    Output, double X[N], W[N], the abscissas and weights of the rule.
*/
{
  int i;
  int j;
  static double pi = 3.141592653589793;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FEJER1_COMPUTE - Fatal error!\n" );
    fprintf ( stderr,"  N < 1.\n" );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
    return;
  }

  double theta[n];

  for ( i = 1; i <= n; i++ )
  {
    theta[i-1] = ( double ) ( 2 * ( n - i ) + 1 ) * pi
               / ( double ) ( 2 * n     );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( theta[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 1.0;
    for ( j = 1; j <= ( n / 2 ); j++ )
    {
      w[i] = w[i] - 2.0
        * cos ( 2.0 * ( double ) ( j ) * theta[i] )
        / ( double ) ( 4 * j * j - 1 );
    }
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 * w[i] / ( double ) ( n );
  }

  return;
}
/******************************************************************************/

void fejer1_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    FEJER1_SET sets abscissas and weights for Fejer type 1 quadrature.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    05 March 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.
  
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.  ORDER should be
    between 1 and 10.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 1 )
  {
    xtab[0]   =  0.000000000000000;
    weight[0] =  2.000000000000000;
  }
  else if ( order == 2 )
  {
    xtab[0] =   -0.7071067811865475;
    xtab[1] =    0.7071067811865475;

    weight[0] =  1.000000000000000;
    weight[1] =  1.000000000000000;
  }
  else if ( order == 3 )
  {
    xtab[0] =  -0.8660254037844387;
    xtab[1] =   0.0000000000000000;
    xtab[2] =   0.8660254037844387;

    weight[0] =  0.4444444444444444;
    weight[1] =  1.111111111111111;
    weight[2] =  0.4444444444444444;
  }
  else if ( order == 4 )
  {
    xtab[0] =  -0.9238795325112867;
    xtab[1] =  -0.3826834323650897;
    xtab[2] =   0.3826834323650898;
    xtab[3] =   0.9238795325112867;

    weight[0] = 0.2642977396044841;
    weight[1] = 0.7357022603955158;
    weight[2] = 0.7357022603955158;
    weight[3] = 0.2642977396044841;
  }
  else if ( order == 5 )
  {
    xtab[0] =  -0.9510565162951535;
    xtab[1] =  -0.5877852522924730;
    xtab[2] =   0.0000000000000000;
    xtab[3] =   0.5877852522924731;
    xtab[4] =   0.9510565162951535;

    weight[0] = 0.1677812284666835;
    weight[1] = 0.5255521048666498;
    weight[2] = 0.6133333333333333;
    weight[3] = 0.5255521048666498;
    weight[4] = 0.1677812284666835;
  }
  else if ( order == 6 )
  {
    xtab[0] =  -0.9659258262890682;
    xtab[1] =  -0.7071067811865475;
    xtab[2] =  -0.2588190451025206;
    xtab[3] =   0.2588190451025207;
    xtab[4] =   0.7071067811865476;
    xtab[5] =   0.9659258262890683;

    weight[0] = 0.1186610213812358;
    weight[1] = 0.3777777777777778;
    weight[2] = 0.5035612008409863;
    weight[3] = 0.5035612008409863;
    weight[4] = 0.3777777777777778;
    weight[5] = 0.1186610213812358;
  }
  else if ( order == 7 )
  {
    xtab[0] =  -0.9749279121818237;
    xtab[1] =  -0.7818314824680295;
    xtab[2] =  -0.4338837391175581;
    xtab[3] =   0.0000000000000000;
    xtab[4] =   0.4338837391175582;
    xtab[5] =   0.7818314824680298;
    xtab[6] =   0.9749279121818236;

    weight[0] = 0.08671618072672234;
    weight[1] = 0.2878313947886919;
    weight[2] = 0.3982415401308441;
    weight[3] = 0.4544217687074830;
    weight[4] = 0.3982415401308441;
    weight[5] = 0.2878313947886919;
    weight[6] = 0.08671618072672234;
  }
  else if ( order == 8 )
  {
    xtab[0] =  -0.9807852804032304;
    xtab[1] =  -0.8314696123025453;
    xtab[2] =  -0.5555702330196020;
    xtab[3] =  -0.1950903220161282;
    xtab[4] =   0.1950903220161283;
    xtab[5] =   0.5555702330196023;
    xtab[6] =   0.8314696123025452;
    xtab[7] =   0.9807852804032304;

    weight[0] = 0.06698294569858981;
    weight[1] = 0.2229879330145788;
    weight[2] = 0.3241525190645244;
    weight[3] = 0.3858766022223071;
    weight[4] = 0.3858766022223071;
    weight[5] = 0.3241525190645244;
    weight[6] = 0.2229879330145788;
    weight[7] = 0.06698294569858981;
 }
 else if ( order == 9 )
 {
    xtab[0] =  -0.9848077530122080;
    xtab[1] =  -0.8660254037844385;
    xtab[2] =  -0.6427876096865394;
    xtab[3] =  -0.3420201433256685;
    xtab[4] =   0.0000000000000000;
    xtab[5] =   0.3420201433256688;
    xtab[6] =   0.6427876096865394;
    xtab[7] =   0.8660254037844387;
    xtab[8] =   0.9848077530122080;

    weight[0] = 0.05273664990990676;
    weight[1] = 0.1791887125220458;
    weight[2] = 0.2640372225410044;
    weight[3] = 0.3308451751681364;
    weight[4] = 0.3463844797178130;
    weight[5] = 0.3308451751681364;
    weight[6] = 0.2640372225410044;
    weight[7] = 0.1791887125220458;
    weight[8] = 0.05273664990990676;
  }
  else if ( order == 10 )
  {
    xtab[0] =  -0.9876883405951377;
    xtab[1] =  -0.8910065241883678;
    xtab[2] =  -0.7071067811865475;
    xtab[3] =  -0.4539904997395467;
    xtab[4] =  -0.1564344650402306;
    xtab[5] =   0.1564344650402309;
    xtab[6] =   0.4539904997395468;
    xtab[7] =   0.7071067811865476;
    xtab[8] =   0.8910065241883679;
    xtab[9] =   0.9876883405951378;

    weight[0] = 0.04293911957413078;
    weight[1] = 0.1458749193773909;
    weight[2] = 0.2203174603174603;
    weight[3] = 0.2808792186638755;
    weight[4] = 0.3099892820671425;
    weight[5] = 0.3099892820671425;
    weight[6] = 0.2808792186638755;
    weight[7] = 0.2203174603174603;
    weight[8] = 0.1458749193773909;
    weight[9] = 0.04293911957413078;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FEJER1_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 through 10.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void fejer2_compute ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    FEJER2_COMPUTE computes a Fejer type 2 quadrature rule.
  
  Discussion:
  
    This method uses a direct approach.  The paper by Waldvogel
    exhibits a more efficient approach using Fourier transforms.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    05 March 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.
  
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  
  Parameters:
  
    Input, int N, the order of the rule.
  
    Output, double X[N], W[N], the abscissas and weights of the rule.
*/
{
  int i;
  int j;
  double p;
  static double pi = 3.141592653589793;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FEJER2_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  N < 1.\n" );
    exit ( 1 );
  }

  if ( n == 1 )
  {
    x[0] = 0.0;
    w[0] = 2.0;
    return;
  }

  double theta[n];

  for ( i = 1; i <= n; i++ )
  {
    theta[i-1] = ( double ) ( n + 1 - i ) * pi
               / ( double ) ( n + 1     );
  }

  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( theta[i] );
  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 1.0;

    for ( j = 1; j <= ( ( n - 1 ) / 2 ); j++ )
    {
      w[i] = w[i] - 2.0 * cos ( 2.0 * ( double ) ( j ) * theta[i] )
        / ( double ) ( 4 * j * j - 1 );
    }

    if ( 2 < n )
    {
      p = 2.0 * ( double ) ( ( n + 1 ) / 2 ) - 1.0;
      w[i] = w[i] - cos ( ( p + 1.0 ) * theta[i] ) / p;
    }

  }

  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 * w[i] / ( double ) ( n + 1 );
  }

  return;
}
/******************************************************************************/

void fejer2_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    FEJER2_SET sets abscissas and weights for Fejer type 2 quadrature.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    05 March 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
    Walter Gautschi,
    Numerical Quadrature in the Presence of a Singularity,
    SIAM Journal on Numerical Analysis,
    Volume 4, Number 3, 1967, pages 357-362.
  
    Joerg Waldvogel,
    Fast Construction of the Fejer and Clenshaw-Curtis Quadrature Rules,
    BIT Numerical Mathematics,
    Volume 43, Number 1, 2003, pages 1-18.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.  ORDER should be
    between 1 and 10.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 1 )
  {
    xtab[0]   =  0.000000000000000;
    weight[0] =  2.000000000000000;
  }
  else if ( order == 2 )
  {
    xtab[0] =   -0.5000000000000000;
    xtab[1] =    0.5000000000000000;

    weight[0] =  1.0000000000000000;
    weight[1] =  1.0000000000000000;
  }
  else if ( order == 3 )
  {
    xtab[0] =  -0.7071067811865476;
    xtab[1] =   0.0000000000000000;
    xtab[2] =   0.7071067811865476;

    weight[0] =  0.6666666666666666;
    weight[1] =  0.6666666666666666;
    weight[2] =  0.6666666666666666;
  }
  else if ( order == 4 )
  {
    xtab[0] =  -0.8090169943749475;
    xtab[1] =  -0.3090169943749475;
    xtab[2] =   0.3090169943749475;
    xtab[3] =   0.8090169943749475;

    weight[0] = 0.4254644007500070;
    weight[1] = 0.5745355992499930;
    weight[2] = 0.5745355992499930;
    weight[3] = 0.4254644007500070;
  }
  else if ( order == 5 )
  {
    xtab[0] =  -0.8660254037844387;
    xtab[1] =  -0.5000000000000000;
    xtab[2] =   0.0000000000000000;
    xtab[3] =   0.5000000000000000;
    xtab[4] =   0.8660254037844387;

    weight[0] = 0.3111111111111111;
    weight[1] = 0.4000000000000000;
    weight[2] = 0.5777777777777777;
    weight[3] = 0.4000000000000000;
    weight[4] = 0.3111111111111111;
  }
  else if ( order == 6 )
  {
    xtab[0] =  -0.9009688679024191;
    xtab[1] =  -0.6234898018587336;
    xtab[2] =  -0.2225209339563144;
    xtab[3] =   0.2225209339563144;
    xtab[4] =   0.6234898018587336;
    xtab[5] =   0.9009688679024191;

    weight[0] = 0.2269152467244296;
    weight[1] = 0.3267938603769863;
    weight[2] = 0.4462908928985841;
    weight[3] = 0.4462908928985841;
    weight[4] = 0.3267938603769863;
    weight[5] = 0.2269152467244296;
  }
  else if ( order == 7 )
  {
    xtab[0] =  -0.9238795325112867;
    xtab[1] =  -0.7071067811865476;
    xtab[2] =  -0.3826834323650898;
    xtab[3] =   0.0000000000000000;
    xtab[4] =   0.3826834323650898;
    xtab[5] =   0.7071067811865476;
    xtab[6] =   0.9238795325112867;

    weight[0] = 0.1779646809620499;
    weight[1] = 0.2476190476190476;
    weight[2] = 0.3934638904665215;
    weight[3] = 0.3619047619047619;
    weight[4] = 0.3934638904665215;
    weight[5] = 0.2476190476190476;
    weight[6] = 0.1779646809620499;
  }
  else if ( order == 8 )
  {
    xtab[0] =  -0.9396926207859084;
    xtab[1] =  -0.7660444431189780;
    xtab[2] =  -0.5000000000000000;
    xtab[3] =  -0.1736481776669304;
    xtab[4] =   0.1736481776669304;
    xtab[5] =   0.5000000000000000;
    xtab[6] =   0.7660444431189780;
    xtab[7] =   0.9396926207859084;

    weight[0] = 0.1397697435050225;
    weight[1] = 0.2063696457302284;
    weight[2] = 0.3142857142857143;
    weight[3] = 0.3395748964790348;
    weight[4] = 0.3395748964790348;
    weight[5] = 0.3142857142857143;
    weight[6] = 0.2063696457302284;
    weight[7] = 0.1397697435050225;
  }
  else if ( order == 9 )
  {
    xtab[0] =  -0.9510565162951535;
    xtab[1] =  -0.8090169943749475;
    xtab[2] =  -0.5877852522924731;
    xtab[3] =  -0.3090169943749475;
    xtab[4] =   0.0000000000000000;
    xtab[5] =   0.3090169943749475;
    xtab[6] =   0.5877852522924731;
    xtab[7] =   0.8090169943749475;
    xtab[8] =   0.9510565162951535;

    weight[0] = 0.1147810750857217;
    weight[1] = 0.1654331942222276;
    weight[2] = 0.2737903534857068;
    weight[3] = 0.2790112502222170;
    weight[4] = 0.3339682539682539;
    weight[5] = 0.2790112502222170;
    weight[6] = 0.2737903534857068;
    weight[7] = 0.1654331942222276;
    weight[8] = 0.1147810750857217;
  }
  else if ( order == 10 )
  {
    xtab[0] =  -0.9594929736144974;
    xtab[1] =  -0.8412535328311812;
    xtab[2] =  -0.6548607339452851;
    xtab[3] =  -0.4154150130018864;
    xtab[4] =  -0.1423148382732851;
    xtab[5] =   0.1423148382732851;
    xtab[6] =   0.4154150130018864;
    xtab[7] =   0.6548607339452851;
    xtab[8] =   0.8412535328311812;
    xtab[9] =   0.9594929736144974;

    weight[0] = 0.09441954173982806;
    weight[1] = 0.1411354380109716;
    weight[2] = 0.2263866903636005;
    weight[3] = 0.2530509772156453;
    weight[4] = 0.2850073526699544;
    weight[5] = 0.2850073526699544;
    weight[6] = 0.2530509772156453;
    weight[7] = 0.2263866903636005;
    weight[8] = 0.1411354380109716;
    weight[9] = 0.09441954173982806;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "FEJER2_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 through 10.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void gegenbauer_compute ( int order, double alpha, double xtab[],
  double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    GEGENBAUER_COMPUTE computes a Gauss-Gegenbauer quadrature rule.
  
  Discussion:
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
    Thanks to Janiki Raman for pointing out a problem in an earlier
    version of the code that occurred when ALPHA was -0.5.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    24 June 2008
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input, int ORDER, the order of the quadrature rule.
  
    Input, double ALPHA, the exponent of (1-X^2) in the weight.
    -1.0 < ALPHA is required.
  
    Output, double XTAB[ORDER], the abscissas.
  
    Output, double WEIGHT[ORDER], the weights.
*/
{
  double an;
  double cc;
  double delta;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double r3;
  double temp;
  double x;
/*
  Check ORDER.
*/
  if ( order < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GEGENBAUER_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  1 <= ORDER is required.\n" );
    exit ( 1 );
  }

  double c[order];
/*
  Check ALPHA.
*/
  if ( alpha <= -1.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "GEGENBAUER_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  -1.0 < ALPHA is required.\n" );
    exit ( 1 );
  }
/*
  Set the recursion coefficients.
*/
  c[0] = 0.0;
  if ( 2 <= order )
  {
    c[1] = 1.0 / ( 2.0 * alpha + 3.0 );
  }

  for ( i = 3; i <= order; i++ )
  {
    c[i-1] = ( double ) ( i - 1 )
          * ( alpha + alpha + ( double ) ( i - 1 ) ) /
          ( ( alpha + alpha + ( double ) ( 2 * i - 1 ) )
          * ( alpha + alpha + ( double ) ( 2 * i - 3 ) ) );
  }

  delta = r8_gamma ( alpha         + 1.0 )
        * r8_gamma (         alpha + 1.0 )
        / r8_gamma ( alpha + alpha + 2.0 );

  prod = 1.0;
  for ( i = 2; i <= order; i++ )
  {
    prod = prod * c[i-1];
  }
  cc = delta * pow ( 2.0, alpha + alpha + 1.0 ) * prod;

  for ( i = 1; i <= order; i++ )
  {
    if ( i == 1 )
    {
      an = alpha / ( double ) ( order );

      r1 = ( 1.0 + alpha )
        * ( 2.78 / ( 4.0 + ( double ) ( order * order ) )
        + 0.768 * an / ( double ) ( order ) );

      r2 = 1.0 + 2.44 * an + 1.282 * an * an;

      x = ( r2 - r1 ) / r2;
    }
    else if ( i == 2 )
    {
      r1 = ( 4.1 + alpha ) /
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( ( double ) ( order ) - 8.0 ) *
        ( 1.0 + 0.12 * alpha ) / ( double ) ( order );

      r3 = 1.0 + 0.012 * alpha *
        ( 1.0 + 0.25 * r8_abs ( alpha ) ) / ( double ) ( order );

      x = x - r1 * r2 * r3 * ( 1.0 - x );
    }
    else if ( i == 3 )
    {
      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 )
        / ( double ) ( order );

      r3 = 1.0 + 8.0 * alpha /
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) );

      x = x - r1 * r2 * r3 * ( xtab[0] - x );
    }
    else if ( i < order - 1 )
    {
      x = 3.0 * xtab[i-2] - 3.0 * xtab[i-3] + xtab[i-4];
    }
    else if ( i == order - 1 )
    {
      r1 = ( 1.0 + 0.235 * alpha ) / ( 0.766 + 0.119 * alpha );

      r2 = 1.0 / ( 1.0 + 0.639
        * ( ( double ) ( order ) - 4.0 )
        / ( 1.0 + 0.71 * ( ( double ) ( order ) - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) *
        ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }
    else if ( i == order )
    {
      r1 = ( 1.0 + 0.37 * alpha ) / ( 1.67 + 0.28 * alpha );

      r2 = 1.0 /
        ( 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 )
        / ( double ) ( order ) );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha /
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }

    gegenbauer_root ( &x, order, alpha, &dp2, &p1, c );

    xtab[i-1] = x;
    weight[i-1] = cc / ( dp2 * p1 );
  }
/*
  Reverse the order of the values.
*/
  for ( i = 1; i <= order/2; i++ )
  {
    temp          = xtab[i-1];
    xtab[i-1]     = xtab[order-i];
    xtab[order-i] = temp;
  }

  for ( i = 1; i <=order/2; i++ )
  {
    temp            = weight[i-1];
    weight[i-1]     = weight[order-i];
    weight[order-i] = temp;
  }

  return;
}
/******************************************************************************/

double gegenbauer_integral ( int expon, double alpha )

/******************************************************************************/
/*
  Purpose:
  
    GEGENBAUER_INTEGRAL evaluates the integral of a monomial with Gegenbauer weight.
  
  Discussion:
  
    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x^2)^ALPHA dx
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 February 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int EXPON, the exponent.
  
    Input, double ALPHA, the exponent of (1-X^2) in the weight factor.
  
    Output, double GEGENBAUER_INTEGRAL, the value of the integral.
*/
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double s;
  double value;
  double value1;

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
    return value;
  }

  c = ( double ) ( expon );

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = r8_gamma ( 1.0 + c ) * 2.0
    * r8_gamma ( 1.0 + alpha  ) * value1
    / r8_gamma ( 2.0 + alpha  + c );

  return value;
}
/******************************************************************************/

void gegenbauer_recur ( double *p2, double *dp2, double *p1, double x, int order,
  double alpha, double c[] )

/******************************************************************************/
/*
  Purpose:
  
    GEGENBAUER_RECUR finds the value and derivative of a Gegenbauer polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 February 2008
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Output, double *P2, the value of J(ORDER)(X).
  
    Output, double *DP2, the value of J'(ORDER)(X).
  
    Output, double *P1, the value of J(ORDER-1)(X).
  
    Input, double X, the point at which polynomials are evaluated.
  
    Input, int ORDER, the order of the polynomial to be computed.
  
    Input, double ALPHA, the exponents of (1-X^2).
  
    Input, double C[ORDER], the recursion coefficients.
*/
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x;
  *dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2 = x *  ( *p1 ) - c[i-1] * p0;
    *dp2 = x * dp1 + ( *p1 ) - c[i-1] * dp0;
  }
  return;
}
/******************************************************************************/

void gegenbauer_root ( double *x, int order, double alpha,  double *dp2,
  double *p1, double c[] )

/******************************************************************************/
/*
  Purpose:
  
    GEGENBAUER_ROOT improves an approximate root of a Gegenbauer polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 February 2008
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input/output, double *X, the approximate root, which
    should be improved on output.
  
    Input, int ORDER, the order of the polynomial to be computed.
  
    Input, double ALPHA, the exponents of (1-X^2).
  
    Output, double *DP2, the value of J'(ORDER)(X).
  
    Output, double *P1, the value of J(ORDER-1)(X).
  
    Input, double C[ORDER], the recursion coefficients.
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    gegenbauer_recur ( &p2, dp2, p1, *x, order, alpha, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }
  return;
}
/******************************************************************************/

void gen_hermite_compute ( int order, double alpha, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule.
  
  Discussion:
  
    The integral to be approximated has the form:
  
      Integral ( -oo < x < +oo ) x^ALPHA exp(-x^2) f(x) dx
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    02 March 2008
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
  
    Input, double ALPHA, the parameter.
  
    Output, double X[ORDER], W[ORDER], the abscissas and weights
    for the requested generalized Gauss-Hermite rule.
*/
{
  double alpha_laguerre;
  double arg;
  int i;
  int order_laguerre;

  if ( order == 1 )
  {
    arg = ( alpha + 1.0 ) / 2.0;
    x[0] = 0.0;
    w[0] = r8_gamma ( arg );
    return;
  }

  if ( ( order % 2 ) == 0 )
  {
    order_laguerre = order / 2;
    alpha_laguerre = ( alpha - 1.0 ) / 2.0;
  }
  else
  {
    order_laguerre = ( order - 1 ) / 2;
    alpha_laguerre = ( alpha + 1.0 ) / 2.0;
  }

   double w_laguerre[order_laguerre];
   double x_laguerre[order_laguerre];

  gen_laguerre_compute ( order_laguerre, alpha_laguerre,
    x_laguerre, w_laguerre );

  if ( ( order % 2 ) == 0 )
  {
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
    }
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[order_laguerre+i] = sqrt ( x_laguerre[i] );
	}
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[i] = 0.5 * w_laguerre[order_laguerre-1-i];
    }
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre+i] = 0.5 * w_laguerre[i];
    }
  }
  else if ( ( order % 2 ) == 1 )
  {
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[i] = - sqrt ( x_laguerre[order_laguerre-1-i] );
    }
    x[order_laguerre] = 0.0;
    for ( i = 0; i < order_laguerre; i++ )
    {
      x[order_laguerre+1+i] = sqrt ( x_laguerre[i] );
	}
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[i] = 0.5 * w_laguerre[order_laguerre-1-i] / x_laguerre[order_laguerre-1-i];
    }

    arg = ( alpha + 1.0 ) / 2.0;
    w[order_laguerre] = r8_gamma ( arg );
    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre] = w[order_laguerre] - w_laguerre[i] / x_laguerre[i];
    }

    for ( i = 0; i < order_laguerre; i++ )
    {
      w[order_laguerre+1+i] = 0.5 * w_laguerre[i] / x_laguerre[i];
    }
  }

  return;
}
/******************************************************************************/

double gen_hermite_integral ( int expon, double alpha )

/******************************************************************************/
/*
  Purpose:
  
    GEN_HERMITE_INTEGRAL evaluates a monomial generalized Hermite integral.
  
  Discussion:
  
    H(n,alpha) = Integral ( -oo < x < +oo ) x^n |x|^alpha exp(-x^2) dx
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    19 February 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int EXPON, the exponent of the monomial.
    0 <= EXPON.
  
    Input, double ALPHA, the exponent of |X| in the weight function.
    -1.0 < ALPHA.
  
    Output, double GEN_HERMITE_INTEGRAL, the value of the integral.
*/
{
  double a;
  double arg;
  double value;

  if ( ( expon % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    a = alpha + ( double ) ( expon );
    if ( a <= - 1.0 )
    {
      value = - r8_huge ( );
    }
    else
    {
      arg = ( a + 1.0 ) / 2.0;
      value = r8_gamma ( arg );
    }
  }
  return value;
}
/******************************************************************************/

void gen_laguerre_compute ( int order, double alpha, double xtab[],
  double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre quadrature rule.
  
  Discussion:
  
    In the simplest case, ALPHA is 0, and we are approximating the
    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
    it is easy to modify the rule to approximate the integral from
    A to +oo as well.
  
    If ALPHA is nonzero, then there is no simple way to extend the
    rule to approximate the integral from A to +oo.  The simplest
    procedures would be to approximate the integral from 0 to A.
  
    The integration interval is [ A, +oo ) or [ 0, +oo ).
  
    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x**alpha.
  
  
    If the integral to approximate is:
  
        Integral ( A <= X < +oo ) EXP ( - X ) * F(X) dX
      or
        Integral ( 0 <= X < +oo ) EXP ( - X ) * X**ALPHA * F(X) dX
  
    then the quadrature rule is:
  
      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( A+XTAB(I) )
    or
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  
    If the integral to approximate is:
  
        Integral ( A <= X < +oo ) F(X) dX
      or
        Integral ( 0 <= X < +oo ) X**ALPHA * F(X) dX
  
    then the quadrature rule is:
  
      EXP ( - A ) * Sum ( 1 <= I <= ORDER )
        WEIGHT(I) * EXP(A+XTAB(I)) * F ( A+XTAB(I) )
    or
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    18 February 2008
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input, int ORDER, the order of the quadrature rule to be computed.
    ORDER must be at least 1.
  
    Input, double ALPHA, the exponent of the X factor.
    Set ALPHA = 0.0 for the simplest rule.
    ALPHA must be nonnegative.
  
    Output, double XTAB[ORDER], the Gauss-Laguerre abscissas.
  
    Output, double WEIGHT[ORDER], the Gauss-Laguerre weights.
*/
{
  double cc;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double ratio;
  double x;

  double b[order];
  double c[order];
/*
  Set the recursion coefficients.
*/
  for ( i = 0; i < order; i++ )
  {
    b[i] = ( alpha + ( double ) ( 2 * i + 1 ) );
  }

  for ( i = 0; i < order; i++ )
  {
    c[i] = ( double ) ( i ) * ( alpha + ( double ) ( i ) );
  }
  prod = 1.0;
  for ( i = 1; i < order; i++ )
  {
    prod = prod * c[i];
  }
  cc = r8_gamma ( alpha + 1.0 ) * prod;

  for ( i = 0; i < order; i++ )
  {
/*
  Compute an estimate for the root.
*/
    if ( i == 0 )
    {
      x = ( 1.0 + alpha ) * ( 3.0+ 0.92 * alpha ) /
        ( 1.0 + 2.4 * ( double ) ( order ) + 1.8 * alpha );
    }
    else if ( i == 1 )
    {
      x = x + ( 15.0 + 6.25 * alpha ) /
        ( 1.0 + 0.9 * alpha + 2.5 * ( double ) ( order ) );
    }
    else
    {
      r1 = ( 1.0 + 2.55 * ( double ) ( i - 1 ) )
        / ( 1.9 * ( double ) ( i - 1 ) );

      r2 = 1.26 * ( double ) ( i - 1 ) * alpha /
        ( 1.0 + 3.5 * ( double ) ( i - 1 ) );

      ratio = ( r1 + r2 ) / ( 1.0 + 0.3 * alpha );

      x = x + ratio * ( x - xtab[i-2] );
    }
/*
  Use iteration to find the root.
*/
    gen_laguerre_root ( &x, order, alpha, &dp2, &p1, b, c );
/*
  Set the abscissa and weight.
*/
    xtab[i] = x;
    weight[i] = ( cc / dp2 ) / p1;
  }

  return;
}
/******************************************************************************/

double gen_laguerre_integral ( int expon, double alpha )

/******************************************************************************/
/*
  Purpose:
  
    GEN_LAGUERRE_INTEGRAL evaluates a monomial generalized Laguerre integral.
  
  Discussion:
  
    L(n,alpha) = Integral ( 0 <= x < +oo ) x^n * x^alpha exp(-x) dx
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    20 February 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int EXPON, the exponent of the monomial.
    0 <= EXPON.
  
    Input, double ALPHA, the exponent of X in the weight function.
    -1.0 < ALPHA.
  
    Output, double GEN_LAGUERRE_INTEGRAL, the value of the integral.
*/
{
  double arg;
  double value;

  arg = alpha + ( double ) ( expon + 1.0 );
  value = r8_gamma ( arg );

  return value;
}
/******************************************************************************/

void gen_laguerre_recur ( double *p2, double *dp2, double *p1, double x,
  int order, double alpha, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:
  
    GEN_LAGUERRE_RECUR evaluates a generalized Laguerre polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    18 February 2008
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Output, double *P2, the value of L(ORDER)(X).
  
    Output, double *DP2, the value of L'(ORDER)(X).
  
    Output, double *P1, the value of L(ORDER-1)(X).
  
    Input, double X, the point at which polynomials are evaluated.
  
    Input, int ORDER, the order of the polynomial to be computed.
  
    Input, double ALPHA, the exponent of the X factor in the
    integrand.
  
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x - alpha - 1.0;
  *dp2 = 1.0;

  for ( i = 1; i < order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
    *dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
  }

  return;
}
/******************************************************************************/

void gen_laguerre_root ( double *x, int order, double alpha, double *dp2,
  double *p1, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:
  
    GEN_LAGUERRE_ROOT improves a root of a generalized Laguerre polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    18 February 2008
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input/output, double *X, the approximate root, which
    should be improved on output.
  
    Input, int ORDER, the order of the polynomial to be computed.
  
    Input, double ALPHA, the exponent of the X factor.
  
    Output, double *DP2, the value of L'(ORDER)(X).
  
    Output, double *P1, the value of L(ORDER-1)(X).
  
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    gen_laguerre_recur ( &p2, dp2, p1, *x, order, alpha, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

void hermite_compute ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    HERMITE_COMPUTE computes a Gauss-Hermite quadrature rule.
  
  Discussion:
  
    The abscissas are the zeros of the N-th order Hermite polynomial.
  
    The integration interval is ( -oo, +oo ).
  
    The weight function is w(x) = exp ( - x*x ).
  
    The integral to approximate:
  
      Integral ( -oo < X < +oo ) exp ( - X*X ) * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    30 April 2006
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input, int ORDER, the order of the formula to be computed.
  
    Output, double XTAB[ORDER], the Gauss-Hermite abscissas.
  
    Output, double WEIGHT[ORDER], the Gauss-Hermite weights.
*/
{
  double cc;
  double dp2;
  int i;
  double p1;
  double pi = 3.141592653589793;
  double s;
  double x;

  cc = sqrt ( pi ) * r8_gamma ( ( double ) ( order ) )
    / pow ( 2.0, order - 1 );

  s = pow ( 2.0 * ( double ) ( order ) + 1.0, 1.0 / 6.0 );

  for ( i = 0; i < ( order + 1 ) / 2; i++ )
  {
    if ( i == 0 )
    {
      x = s * s * s - 1.85575 / s;
    }
    else if ( i == 1 )
    {
      x = x - 1.14 * pow ( ( double ) ( order ), 0.426 ) / x;
    }
    else if ( i == 2 )
    {
      x = 1.86 * x - 0.86 * xtab[0];
    }
    else if ( i == 3 )
    {
      x = 1.91 * x - 0.91 * xtab[1];
    }
    else
    {
      x = 2.0 * x - xtab[i-2];
    }

    hermite_root ( &x, order, &dp2, &p1 );

    xtab[i] = x;
    weight[i] = ( cc / dp2 ) / p1;

    xtab[order-i-1] = -x;
    weight[order-i-1] = weight[i];
  }
/*
  Reverse the order of the values.
*/
  for ( i = 0; i < ( order / 2 ); i++ )
  {
    x               = xtab[i];
    xtab[i]         = xtab[order-1-i];
    xtab[order-1-i] = x;
  }

  return;
}
/******************************************************************************/

double hermite_integral ( int n )

/******************************************************************************/
/*
  Purpose:
  
    HERMITE_INTEGRAL evaluates a monomial Hermite integral.
  
  Discussion:
  
    H(n) = Integral ( -oo < x < +oo ) x^n exp(-x^2) dx
  
    H(n) is 0 for n odd.
  
    H(n) = (n-1)!! * sqrt(pi) / 2^(n/2) for n even.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    19 February 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the order of the integral.
    0 <= N.
  
    Output, double VALUE, the value of the integral.
*/
{
  double pi = 3.141592653589793;
  double value;

  if ( n < 0 )
  {
    value = - r8_huge ( );
  }
  else if ( ( n % 2 ) == 1 )
  {
    value = 0.0;
  }
  else
  {
    value = r8_factorial2 ( n - 1 ) * sqrt ( pi ) / pow ( 2.0, n / 2 );
  }

  return value;
}
/******************************************************************************/

void hermite_recur ( double *p2, double *dp2, double *p1, double x, int order )

/******************************************************************************/
/*
  Purpose:
  
    HERMITE_RECUR finds the value and derivative of a Hermite polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    30 April 2006
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Output, double *P2, the value of H(ORDER)(X).
  
    Output, double *DP2, the value of H'(ORDER)(X).
  
    Output, double *P1, the value of H(ORDER-1)(X).
  
    Input, double X, the point at which polynomials are evaluated.
  
    Input, int ORDER, the order of the polynomial to be computed.
*/
{
  int i;
  double dq0;
  double dq1;
  double dq2;
  double q0;
  double q1;
  double q2;

  q1 = 1.0;
  dq1 = 0.0;

  q2 = x;
  dq2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    q0 = q1;
    dq0 = dq1;

    q1 = q2;
    dq1 = dq2;

    q2  = x * q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * q0;
    dq2 = x * dq1 + q1 - 0.5 * ( ( double ) ( i ) - 1.0 ) * dq0;
  }

  *p2 = q2;
  *dp2 = dq2;
  *p1 = q1;

  return;
}
/******************************************************************************/

void hermite_root ( double *x, int order, double *dp2, double *p1 )

/******************************************************************************/
/*
  Purpose:
  
    HERMITE_ROOT improves an approximate root of a Hermite polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    19 February 2008
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input/output, double *X, the approximate root, which
    should be improved on output.
  
    Input, int ORDER, the order of the Hermite polynomial.
  
    Output, double *DP2, the value of H'(ORDER)(X).
  
    Output, double *P1, the value of H(ORDER-1)(X).
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    hermite_recur ( &p2, dp2, p1, *x, order );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }

  return;
}
/******************************************************************************/

void hermite_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    HERMITE_SET sets abscissas and weights for Hermite quadrature.
  
  Discussion:
  
    The integration interval is ( -oo, +oo ).
  
    The weight function is w(x) = exp ( - x**2 ).
  
    The integral to approximate:
  
      Integral ( -oo < X < +oo ) exp ( - X**2 ) * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) ).
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    03 October 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Vladimir Krylov,
    Approximate Calculation of Integrals,
    Dover, 2006,
    ISBN: 0486445798,
    LC: QA311.K713.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3,
    LC: QA47.M315.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 20, or one of the values
    30, 31, 32, 40, 50, 60, 63, 64 or 127.
  
    Output, double XTAB[ORDER], the abscissas of the rule,
    which are symmetrically placed around 0.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are positive and symmetric, and should sum
    to SQRT(PI).
*/
{
  if ( order == 1 )
  {
    xtab[1-1] = 0.0;

    weight[1-1] = 1.77245385090551602729816748334;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = - 0.707106781186547524400844362105E+00;
    xtab[2-1] =   0.707106781186547524400844362105E+00;

    weight[1-1] = 0.886226925452758013649083741671E+00;
    weight[2-1] = 0.886226925452758013649083741671E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = - 0.122474487139158904909864203735E+01;
    xtab[2-1] =   0.0E+00;
    xtab[3-1] =   0.122474487139158904909864203735E+01;

    weight[1-1] = 0.295408975150919337883027913890E+00;
    weight[2-1] = 0.118163590060367735153211165556E+01;
    weight[3-1] = 0.295408975150919337883027913890E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = - 0.165068012388578455588334111112E+01;
    xtab[2-1] = - 0.524647623275290317884060253835E+00;
    xtab[3-1] =   0.524647623275290317884060253835E+00;
    xtab[4-1] =   0.165068012388578455588334111112E+01;

    weight[1-1] = 0.813128354472451771430345571899E-01;
    weight[2-1] = 0.804914090005512836506049184481E+00;
    weight[3-1] = 0.804914090005512836506049184481E+00;
    weight[4-1] = 0.813128354472451771430345571899E-01;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = - 0.202018287045608563292872408814E+01;
    xtab[2-1] = - 0.958572464613818507112770593893E+00;
    xtab[3-1] =   0.0E+00;
    xtab[4-1] =   0.958572464613818507112770593893E+00;
    xtab[5-1] =   0.202018287045608563292872408814E+01;

    weight[1-1] = 0.199532420590459132077434585942E-01;
    weight[2-1] = 0.393619323152241159828495620852E+00;
    weight[3-1] = 0.945308720482941881225689324449E+00;
    weight[4-1] = 0.393619323152241159828495620852E+00;
    weight[5-1] = 0.199532420590459132077434585942E-01;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = - 0.235060497367449222283392198706E+01;
    xtab[2-1] = - 0.133584907401369694971489528297E+01;
    xtab[3-1] = - 0.436077411927616508679215948251E+00;
    xtab[4-1] =   0.436077411927616508679215948251E+00;
    xtab[5-1] =   0.133584907401369694971489528297E+01;
    xtab[6-1] =   0.235060497367449222283392198706E+01;

    weight[1-1] = 0.453000990550884564085747256463E-02;
    weight[2-1] = 0.157067320322856643916311563508E+00;
    weight[3-1] = 0.724629595224392524091914705598E+00;
    weight[4-1] = 0.724629595224392524091914705598E+00;
    weight[5-1] = 0.157067320322856643916311563508E+00;
    weight[6-1] = 0.453000990550884564085747256463E-02;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = - 0.265196135683523349244708200652E+01;
    xtab[2-1] = - 0.167355162876747144503180139830E+01;
    xtab[3-1] = - 0.816287882858964663038710959027E+00;
    xtab[4-1] =   0.0E+00;
    xtab[5-1] =   0.816287882858964663038710959027E+00;
    xtab[6-1] =   0.167355162876747144503180139830E+01;
    xtab[7-1] =   0.265196135683523349244708200652E+01;

    weight[1-1] = 0.971781245099519154149424255939E-03;
    weight[2-1] = 0.545155828191270305921785688417E-01;
    weight[3-1] = 0.425607252610127800520317466666E+00;
    weight[4-1] = 0.810264617556807326764876563813E+00;
    weight[5-1] = 0.425607252610127800520317466666E+00;
    weight[6-1] = 0.545155828191270305921785688417E-01;
    weight[7-1] = 0.971781245099519154149424255939E-03;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = - 0.293063742025724401922350270524E+01;
    xtab[2-1] = - 0.198165675669584292585463063977E+01;
    xtab[3-1] = - 0.115719371244678019472076577906E+01;
    xtab[4-1] = - 0.381186990207322116854718885584E+00;
    xtab[5-1] =   0.381186990207322116854718885584E+00;
    xtab[6-1] =   0.115719371244678019472076577906E+01;
    xtab[7-1] =   0.198165675669584292585463063977E+01;
    xtab[8-1] =   0.293063742025724401922350270524E+01;

    weight[1-1] = 0.199604072211367619206090452544E-03;
    weight[2-1] = 0.170779830074134754562030564364E-01;
    weight[3-1] = 0.207802325814891879543258620286E+00;
    weight[4-1] = 0.661147012558241291030415974496E+00;
    weight[5-1] = 0.661147012558241291030415974496E+00;
    weight[6-1] = 0.207802325814891879543258620286E+00;
    weight[7-1] = 0.170779830074134754562030564364E-01;
    weight[8-1] = 0.199604072211367619206090452544E-03;
  }
  else if ( order == 9 )
  {
    xtab[1-1] = - 0.319099320178152760723004779538E+01;
    xtab[2-1] = - 0.226658058453184311180209693284E+01;
    xtab[3-1] = - 0.146855328921666793166701573925E+01;
    xtab[4-1] = - 0.723551018752837573322639864579E+00;
    xtab[5-1] =   0.0E+00;
    xtab[6-1] =   0.723551018752837573322639864579E+00;
    xtab[7-1] =   0.146855328921666793166701573925E+01;
    xtab[8-1] =   0.226658058453184311180209693284E+01;
    xtab[9-1] =   0.319099320178152760723004779538E+01;

    weight[1-1] = 0.396069772632643819045862946425E-04;
    weight[2-1] = 0.494362427553694721722456597763E-02;
    weight[3-1] = 0.884745273943765732879751147476E-01;
    weight[4-1] = 0.432651559002555750199812112956E+00;
    weight[5-1] = 0.720235215606050957124334723389E+00;
    weight[6-1] = 0.432651559002555750199812112956E+00;
    weight[7-1] = 0.884745273943765732879751147476E-01;
    weight[8-1] = 0.494362427553694721722456597763E-02;
    weight[9-1] = 0.396069772632643819045862946425E-04;
  }
  else if ( order == 10 )
  {
    xtab[1-1] =  - 0.343615911883773760332672549432E+01;
    xtab[2-1] =  - 0.253273167423278979640896079775E+01;
    xtab[3-1] =  - 0.175668364929988177345140122011E+01;
    xtab[4-1] =  - 0.103661082978951365417749191676E+01;
    xtab[5-1] =  - 0.342901327223704608789165025557E+00;
    xtab[6-1] =    0.342901327223704608789165025557E+00;
    xtab[7-1] =    0.103661082978951365417749191676E+01;
    xtab[8-1] =    0.175668364929988177345140122011E+01;
    xtab[9-1] =    0.253273167423278979640896079775E+01;
    xtab[10-1] =   0.343615911883773760332672549432E+01;

    weight[1-1] =  0.764043285523262062915936785960E-05;
    weight[2-1] =  0.134364574678123269220156558585E-02;
    weight[3-1] =  0.338743944554810631361647312776E-01;
    weight[4-1] =  0.240138611082314686416523295006E+00;
    weight[5-1] =  0.610862633735325798783564990433E+00;
    weight[6-1] =  0.610862633735325798783564990433E+00;
    weight[7-1] =  0.240138611082314686416523295006E+00;
    weight[8-1] =  0.338743944554810631361647312776E-01;
    weight[9-1] =  0.134364574678123269220156558585E-02;
    weight[10-1] = 0.764043285523262062915936785960E-05;
  }
  else if ( order == 11 )
  {
    xtab[1-1] =  - 0.366847084655958251845837146485E+01;
    xtab[2-1] =  - 0.278329009978165177083671870152E+01;
    xtab[3-1] =  - 0.202594801582575533516591283121E+01;
    xtab[4-1] =  - 0.132655708449493285594973473558E+01;
    xtab[5-1] =  - 0.656809566882099765024611575383E+00;
    xtab[6-1] =    0.0E+00;
    xtab[7-1] =    0.656809566882099765024611575383E+00;
    xtab[8-1] =    0.132655708449493285594973473558E+01;
    xtab[9-1] =    0.202594801582575533516591283121E+01;
    xtab[10-1] =   0.278329009978165177083671870152E+01;
    xtab[11-1] =   0.366847084655958251845837146485E+01;

    weight[1-1] =  0.143956039371425822033088366032E-05;
    weight[2-1] =  0.346819466323345510643413772940E-03;
    weight[3-1] =  0.119113954449115324503874202916E-01;
    weight[4-1] =  0.117227875167708503381788649308E+00;
    weight[5-1] =  0.429359752356125028446073598601E+00;
    weight[6-1] =  0.654759286914591779203940657627E+00;
    weight[7-1] =  0.429359752356125028446073598601E+00;
    weight[8-1] =  0.117227875167708503381788649308E+00;
    weight[9-1] =  0.119113954449115324503874202916E-01;
    weight[10-1] = 0.346819466323345510643413772940E-03;
    weight[11-1] = 0.143956039371425822033088366032E-05;
  }
  else if ( order == 12 )
  {
    xtab[1-1] =  - 0.388972489786978191927164274724E+01;
    xtab[2-1] =  - 0.302063702512088977171067937518E+01;
    xtab[3-1] =  - 0.227950708050105990018772856942E+01;
    xtab[4-1] =  - 0.159768263515260479670966277090E+01;
    xtab[5-1] =  - 0.947788391240163743704578131060E+00;
    xtab[6-1] =  - 0.314240376254359111276611634095E+00;
    xtab[7-1] =    0.314240376254359111276611634095E+00;
    xtab[8-1] =    0.947788391240163743704578131060E+00;
    xtab[9-1] =    0.159768263515260479670966277090E+01;
    xtab[10-1] =   0.227950708050105990018772856942E+01;
    xtab[11-1] =   0.302063702512088977171067937518E+01;
    xtab[12-1] =   0.388972489786978191927164274724E+01;

    weight[1-1] =  0.265855168435630160602311400877E-06;
    weight[2-1] =  0.857368704358785865456906323153E-04;
    weight[3-1] =  0.390539058462906185999438432620E-02;
    weight[4-1] =  0.516079856158839299918734423606E-01;
    weight[5-1] =  0.260492310264161129233396139765E+00;
    weight[6-1] =  0.570135236262479578347113482275E+00;
    weight[7-1] =  0.570135236262479578347113482275E+00;
    weight[8-1] =  0.260492310264161129233396139765E+00;
    weight[9-1] =  0.516079856158839299918734423606E-01;
    weight[10-1] = 0.390539058462906185999438432620E-02;
    weight[11-1] = 0.857368704358785865456906323153E-04;
    weight[12-1] = 0.265855168435630160602311400877E-06;
  }
  else if ( order == 13 )
  {
    xtab[1-1] =  - 0.410133759617863964117891508007E+01;
    xtab[2-1] =  - 0.324660897837240998812205115236E+01;
    xtab[3-1] =  - 0.251973568567823788343040913628E+01;
    xtab[4-1] =  - 0.185310765160151214200350644316E+01;
    xtab[5-1] =  - 0.122005503659074842622205526637E+01;
    xtab[6-1] =  - 0.605763879171060113080537108602E+00;
    xtab[7-1] =    0.0E+00;
    xtab[8-1] =    0.605763879171060113080537108602E+00;
    xtab[9-1] =    0.122005503659074842622205526637E+01;
    xtab[10-1] =   0.185310765160151214200350644316E+01;
    xtab[11-1] =   0.251973568567823788343040913628E+01;
    xtab[12-1] =   0.324660897837240998812205115236E+01;
    xtab[13-1] =   0.410133759617863964117891508007E+01;

    weight[1-1] =  0.482573185007313108834997332342E-07;
    weight[2-1] =  0.204303604027070731248669432937E-04;
    weight[3-1] =  0.120745999271938594730924899224E-02;
    weight[4-1] =  0.208627752961699392166033805050E-01;
    weight[5-1] =  0.140323320687023437762792268873E+00;
    weight[6-1] =  0.421616296898543221746893558568E+00;
    weight[7-1] =  0.604393187921161642342099068579E+00;
    weight[8-1] =  0.421616296898543221746893558568E+00;
    weight[9-1] =  0.140323320687023437762792268873E+00;
    weight[10-1] = 0.208627752961699392166033805050E-01;
    weight[11-1] = 0.120745999271938594730924899224E-02;
    weight[12-1] = 0.204303604027070731248669432937E-04;
    weight[13-1] = 0.482573185007313108834997332342E-07;
  }
  else if ( order == 14 )
  {
    xtab[1-1] =  - 0.430444857047363181262129810037E+01;
    xtab[2-1] =  - 0.346265693360227055020891736115E+01;
    xtab[3-1] =  - 0.274847072498540256862499852415E+01;
    xtab[4-1] =  - 0.209518325850771681573497272630E+01;
    xtab[5-1] =  - 0.147668273114114087058350654421E+01;
    xtab[6-1] =  - 0.878713787329399416114679311861E+00;
    xtab[7-1] =  - 0.291745510672562078446113075799E+00;
    xtab[8-1] =    0.291745510672562078446113075799E+00;
    xtab[9-1] =    0.878713787329399416114679311861E+00;
    xtab[10-1] =   0.147668273114114087058350654421E+01;
    xtab[11-1] =   0.209518325850771681573497272630E+01;
    xtab[12-1] =   0.274847072498540256862499852415E+01;
    xtab[13-1] =   0.346265693360227055020891736115E+01;
    xtab[14-1] =   0.430444857047363181262129810037E+01;

    weight[1-1] =  0.862859116812515794532041783429E-08;
    weight[2-1] =  0.471648435501891674887688950105E-05;
    weight[3-1] =  0.355092613551923610483661076691E-03;
    weight[4-1] =  0.785005472645794431048644334608E-02;
    weight[5-1] =  0.685055342234652055387163312367E-01;
    weight[6-1] =  0.273105609064246603352569187026E+00;
    weight[7-1] =  0.536405909712090149794921296776E+00;
    weight[8-1] =  0.536405909712090149794921296776E+00;
    weight[9-1] =  0.273105609064246603352569187026E+00;
    weight[10-1] = 0.685055342234652055387163312367E-01;
    weight[11-1] = 0.785005472645794431048644334608E-02;
    weight[12-1] = 0.355092613551923610483661076691E-03;
    weight[13-1] = 0.471648435501891674887688950105E-05;
    weight[14-1] = 0.862859116812515794532041783429E-08;
  }
  else if ( order == 15 )
  {
    xtab[1-1] =  - 0.449999070730939155366438053053E+01;
    xtab[2-1] =  - 0.366995037340445253472922383312E+01;
    xtab[3-1] =  - 0.296716692790560324848896036355E+01;
    xtab[4-1] =  - 0.232573248617385774545404479449E+01;
    xtab[5-1] =  - 0.171999257518648893241583152515E+01;
    xtab[6-1] =  - 0.113611558521092066631913490556E+01;
    xtab[7-1] =  - 0.565069583255575748526020337198E+00;
    xtab[8-1] =    0.0E+00;
    xtab[9-1] =    0.565069583255575748526020337198E+00;
    xtab[10-1] =   0.113611558521092066631913490556E+01;
    xtab[11-1] =   0.171999257518648893241583152515E+01;
    xtab[12-1] =   0.232573248617385774545404479449E+01;
    xtab[13-1] =   0.296716692790560324848896036355E+01;
    xtab[14-1] =   0.366995037340445253472922383312E+01;
    xtab[15-1] =   0.449999070730939155366438053053E+01;

    weight[1-1] =  0.152247580425351702016062666965E-08;
    weight[2-1] =  0.105911554771106663577520791055E-05;
    weight[3-1] =  0.100004441232499868127296736177E-03;
    weight[4-1] =  0.277806884291277589607887049229E-02;
    weight[5-1] =  0.307800338725460822286814158758E-01;
    weight[6-1] =  0.158488915795935746883839384960E+00;
    weight[7-1] =  0.412028687498898627025891079568E+00;
    weight[8-1] =  0.564100308726417532852625797340E+00;
    weight[9-1] =  0.412028687498898627025891079568E+00;
    weight[10-1] = 0.158488915795935746883839384960E+00;
    weight[11-1] = 0.307800338725460822286814158758E-01;
    weight[12-1] = 0.277806884291277589607887049229E-02;
    weight[13-1] = 0.100004441232499868127296736177E-03;
    weight[14-1] = 0.105911554771106663577520791055E-05;
    weight[15-1] = 0.152247580425351702016062666965E-08;
  }
  else if ( order == 16 )
  {
    xtab[1-1] =  - 0.468873893930581836468849864875E+01;
    xtab[2-1] =  - 0.386944790486012269871942409801E+01;
    xtab[3-1] =  - 0.317699916197995602681399455926E+01;
    xtab[4-1] =  - 0.254620215784748136215932870545E+01;
    xtab[5-1] =  - 0.195178799091625397743465541496E+01;
    xtab[6-1] =  - 0.138025853919888079637208966969E+01;
    xtab[7-1] =  - 0.822951449144655892582454496734E+00;
    xtab[8-1] =  - 0.273481046138152452158280401965E+00;
    xtab[9-1] =    0.273481046138152452158280401965E+00;
    xtab[10-1] =   0.822951449144655892582454496734E+00;
    xtab[11-1] =   0.138025853919888079637208966969E+01;
    xtab[12-1] =   0.195178799091625397743465541496E+01;
    xtab[13-1] =   0.254620215784748136215932870545E+01;
    xtab[14-1] =   0.317699916197995602681399455926E+01;
    xtab[15-1] =   0.386944790486012269871942409801E+01;
    xtab[16-1] =   0.468873893930581836468849864875E+01;

    weight[1-1] =  0.265480747401118224470926366050E-09;
    weight[2-1] =  0.232098084486521065338749423185E-06;
    weight[3-1] =  0.271186009253788151201891432244E-04;
    weight[4-1] =  0.932284008624180529914277305537E-03;
    weight[5-1] =  0.128803115355099736834642999312E-01;
    weight[6-1] =  0.838100413989858294154207349001E-01;
    weight[7-1] =  0.280647458528533675369463335380E+00;
    weight[8-1] =  0.507929479016613741913517341791E+00;
    weight[9-1] =  0.507929479016613741913517341791E+00;
    weight[10-1] = 0.280647458528533675369463335380E+00;
    weight[11-1] = 0.838100413989858294154207349001E-01;
    weight[12-1] = 0.128803115355099736834642999312E-01;
    weight[13-1] = 0.932284008624180529914277305537E-03;
    weight[14-1] = 0.271186009253788151201891432244E-04;
    weight[15-1] = 0.232098084486521065338749423185E-06;
    weight[16-1] = 0.265480747401118224470926366050E-09;
  }
  else if ( order == 17 )
  {
    xtab[1-1] =  - 0.487134519367440308834927655662E+01;
    xtab[2-1] =  - 0.406194667587547430689245559698E+01;
    xtab[3-1] =  - 0.337893209114149408338327069289E+01;
    xtab[4-1] =  - 0.275776291570388873092640349574E+01;
    xtab[5-1] =  - 0.217350282666662081927537907149E+01;
    xtab[6-1] =  - 0.161292431422123133311288254454E+01;
    xtab[7-1] =  - 0.106764872574345055363045773799E+01;
    xtab[8-1] =  - 0.531633001342654731349086553718E+00;
    xtab[9-1] =    0.0E+00;
    xtab[10-1] =   0.531633001342654731349086553718E+00;
    xtab[11-1] =   0.106764872574345055363045773799E+01;
    xtab[12-1] =   0.161292431422123133311288254454E+01;
    xtab[13-1] =   0.217350282666662081927537907149E+01;
    xtab[14-1] =   0.275776291570388873092640349574E+01;
    xtab[15-1] =   0.337893209114149408338327069289E+01;
    xtab[16-1] =   0.406194667587547430689245559698E+01;
    xtab[17-1] =   0.487134519367440308834927655662E+01;

    weight[1-1] =  0.458057893079863330580889281222E-10;
    weight[2-1] =  0.497707898163079405227863353715E-07;
    weight[3-1] =  0.711228914002130958353327376218E-05;
    weight[4-1] =  0.298643286697753041151336643059E-03;
    weight[5-1] =  0.506734995762753791170069495879E-02;
    weight[6-1] =  0.409200341495762798094994877854E-01;
    weight[7-1] =  0.172648297670097079217645196219E+00;
    weight[8-1] =  0.401826469470411956577635085257E+00;
    weight[9-1] =  0.530917937624863560331883103379E+00;
    weight[10-1] = 0.401826469470411956577635085257E+00;
    weight[11-1] = 0.172648297670097079217645196219E+00;
    weight[12-1] = 0.409200341495762798094994877854E-01;
    weight[13-1] = 0.506734995762753791170069495879E-02;
    weight[14-1] = 0.298643286697753041151336643059E-03;
    weight[15-1] = 0.711228914002130958353327376218E-05;
    weight[16-1] = 0.497707898163079405227863353715E-07;
    weight[17-1] = 0.458057893079863330580889281222E-10;
  }
  else if ( order == 18 )
  {
    xtab[1-1] =  - 0.504836400887446676837203757885E+01;
    xtab[2-1] =  - 0.424811787356812646302342016090E+01;
    xtab[3-1] =  - 0.357376906848626607950067599377E+01;
    xtab[4-1] =  - 0.296137750553160684477863254906E+01;
    xtab[5-1] =  - 0.238629908916668600026459301424E+01;
    xtab[6-1] =  - 0.183553160426162889225383944409E+01;
    xtab[7-1] =  - 0.130092085838961736566626555439E+01;
    xtab[8-1] =  - 0.776682919267411661316659462284E+00;
    xtab[9-1] =  - 0.258267750519096759258116098711E+00;
    xtab[10-1] =   0.258267750519096759258116098711E+00;
    xtab[11-1] =   0.776682919267411661316659462284E+00;
    xtab[12-1] =   0.130092085838961736566626555439E+01;
    xtab[13-1] =   0.183553160426162889225383944409E+01;
    xtab[14-1] =   0.238629908916668600026459301424E+01;
    xtab[15-1] =   0.296137750553160684477863254906E+01;
    xtab[16-1] =   0.357376906848626607950067599377E+01;
    xtab[17-1] =   0.424811787356812646302342016090E+01;
    xtab[18-1] =   0.504836400887446676837203757885E+01;

    weight[1-1] =  0.782819977211589102925147471012E-11;
    weight[2-1] =  0.104672057957920824443559608435E-07;
    weight[3-1] =  0.181065448109343040959702385911E-05;
    weight[4-1] =  0.918112686792940352914675407371E-04;
    weight[5-1] =  0.188852263026841789438175325426E-02;
    weight[6-1] =  0.186400423875446519219315221973E-01;
    weight[7-1] =  0.973017476413154293308537234155E-01;
    weight[8-1] =  0.284807285669979578595606820713E+00;
    weight[9-1] =  0.483495694725455552876410522141E+00;
    weight[10-1] = 0.483495694725455552876410522141E+00;
    weight[11-1] = 0.284807285669979578595606820713E+00;
    weight[12-1] = 0.973017476413154293308537234155E-01;
    weight[13-1] = 0.186400423875446519219315221973E-01;
    weight[14-1] = 0.188852263026841789438175325426E-02;
    weight[15-1] = 0.918112686792940352914675407371E-04;
    weight[16-1] = 0.181065448109343040959702385911E-05;
    weight[17-1] = 0.104672057957920824443559608435E-07;
    weight[18-1] = 0.782819977211589102925147471012E-11;
  }
  else if ( order == 19 )
  {
    xtab[1-1] =  - 0.522027169053748216460967142500E+01;
    xtab[2-1] =  - 0.442853280660377943723498532226E+01;
    xtab[3-1] =  - 0.376218735196402009751489394104E+01;
    xtab[4-1] =  - 0.315784881834760228184318034120E+01;
    xtab[5-1] =  - 0.259113378979454256492128084112E+01;
    xtab[6-1] =  - 0.204923170985061937575050838669E+01;
    xtab[7-1] =  - 0.152417061939353303183354859367E+01;
    xtab[8-1] =  - 0.101036838713431135136859873726E+01;
    xtab[9-1] =  - 0.503520163423888209373811765050E+00;
    xtab[10-1] =   0.0E+00;
    xtab[11-1] =   0.503520163423888209373811765050E+00;
    xtab[12-1] =   0.101036838713431135136859873726E+01;
    xtab[13-1] =   0.152417061939353303183354859367E+01;
    xtab[14-1] =   0.204923170985061937575050838669E+01;
    xtab[15-1] =   0.259113378979454256492128084112E+01;
    xtab[16-1] =   0.315784881834760228184318034120E+01;
    xtab[17-1] =   0.376218735196402009751489394104E+01;
    xtab[18-1] =   0.442853280660377943723498532226E+01;
    xtab[19-1] =   0.522027169053748216460967142500E+01;

    weight[1-1] =  0.132629709449851575185289154385E-11;
    weight[2-1] =  0.216305100986355475019693077221E-08;
    weight[3-1] =  0.448824314722312295179447915594E-06;
    weight[4-1] =  0.272091977631616257711941025214E-04;
    weight[5-1] =  0.670877521407181106194696282100E-03;
    weight[6-1] =  0.798886677772299020922211491861E-02;
    weight[7-1] =  0.508103869090520673569908110358E-01;
    weight[8-1] =  0.183632701306997074156148485766E+00;
    weight[9-1] =  0.391608988613030244504042313621E+00;
    weight[10-1] = 0.502974888276186530840731361096E+00;
    weight[11-1] = 0.391608988613030244504042313621E+00;
    weight[12-1] = 0.183632701306997074156148485766E+00;
    weight[13-1] = 0.508103869090520673569908110358E-01;
    weight[14-1] = 0.798886677772299020922211491861E-02;
    weight[15-1] = 0.670877521407181106194696282100E-03;
    weight[16-1] = 0.272091977631616257711941025214E-04;
    weight[17-1] = 0.448824314722312295179447915594E-06;
    weight[18-1] = 0.216305100986355475019693077221E-08;
    weight[19-1] = 0.132629709449851575185289154385E-11;
  }
  else if ( order == 20 )
  {
    xtab[1-1] =  - 0.538748089001123286201690041068E+01;
    xtab[2-1] =  - 0.460368244955074427307767524898E+01;
    xtab[3-1] =  - 0.394476404011562521037562880052E+01;
    xtab[4-1] =  - 0.334785456738321632691492452300E+01;
    xtab[5-1] =  - 0.278880605842813048052503375640E+01;
    xtab[6-1] =  - 0.225497400208927552308233334473E+01;
    xtab[7-1] =  - 0.173853771211658620678086566214E+01;
    xtab[8-1] =  - 0.123407621539532300788581834696E+01;
    xtab[9-1] =  - 0.737473728545394358705605144252E+00;
    xtab[10-1] = - 0.245340708300901249903836530634E+00;
    xtab[11-1] =   0.245340708300901249903836530634E+00;
    xtab[12-1] =   0.737473728545394358705605144252E+00;
    xtab[13-1] =   0.123407621539532300788581834696E+01;
    xtab[14-1] =   0.173853771211658620678086566214E+01;
    xtab[15-1] =   0.225497400208927552308233334473E+01;
    xtab[16-1] =   0.278880605842813048052503375640E+01;
    xtab[17-1] =   0.334785456738321632691492452300E+01;
    xtab[18-1] =   0.394476404011562521037562880052E+01;
    xtab[19-1] =   0.460368244955074427307767524898E+01;
    xtab[20-1] =   0.538748089001123286201690041068E+01;

    weight[1-1] =  0.222939364553415129252250061603E-12;
    weight[2-1] =  0.439934099227318055362885145547E-09;
    weight[3-1] =  0.108606937076928169399952456345E-06;
    weight[4-1] =  0.780255647853206369414599199965E-05;
    weight[5-1] =  0.228338636016353967257145917963E-03;
    weight[6-1] =  0.324377334223786183218324713235E-02;
    weight[7-1] =  0.248105208874636108821649525589E-01;
    weight[8-1] =  0.109017206020023320013755033535E+00;
    weight[9-1] =  0.286675505362834129719659706228E+00;
    weight[10-1] = 0.462243669600610089650328639861E+00;
    weight[11-1] = 0.462243669600610089650328639861E+00;
    weight[12-1] = 0.286675505362834129719659706228E+00;
    weight[13-1] = 0.109017206020023320013755033535E+00;
    weight[14-1] = 0.248105208874636108821649525589E-01;
    weight[15-1] = 0.324377334223786183218324713235E-02;
    weight[16-1] = 0.228338636016353967257145917963E-03;
    weight[17-1] = 0.780255647853206369414599199965E-05;
    weight[18-1] = 0.108606937076928169399952456345E-06;
    weight[19-1] = 0.439934099227318055362885145547E-09;
    weight[20-1] = 0.222939364553415129252250061603E-12;
  }
  else if ( order == 30 )
  {
    xtab[ 1-1] =   -6.86334529352989158106110835756E+00;
    xtab[ 2-1] =   -6.13827922012393462039499237854E+00;
    xtab[ 3-1] =   -5.53314715156749572511833355558E+00;
    xtab[ 4-1] =   -4.98891896858994394448649710633E+00;
    xtab[ 5-1] =   -4.48305535709251834188703761971E+00;
    xtab[ 6-1] =   -4.00390860386122881522787601332E+00;
    xtab[ 7-1] =   -3.54444387315534988692540090217E+00;
    xtab[ 8-1] =   -3.09997052958644174868873332237E+00;
    xtab[ 9-1] =   -2.66713212453561720057110646422E+00;
    xtab[10-1] =   -2.24339146776150407247297999483E+00;
    xtab[11-1] =   -1.82674114360368803883588048351E+00;
    xtab[12-1] =   -1.41552780019818851194072510555E+00;
    xtab[13-1] =   -1.00833827104672346180498960870E+00;
    xtab[14-1] =   -0.603921058625552307778155678757E+00;
    xtab[15-1] =   -0.201128576548871485545763013244E+00;
    xtab[16-1] =    0.201128576548871485545763013244E+00;
    xtab[17-1] =    0.603921058625552307778155678757E+00;
    xtab[18-1] =    1.00833827104672346180498960870E+00;
    xtab[19-1] =    1.41552780019818851194072510555E+00;
    xtab[20-1] =    1.82674114360368803883588048351E+00;
    xtab[21-1] =    2.24339146776150407247297999483E+00;
    xtab[22-1] =    2.66713212453561720057110646422E+00;
    xtab[23-1] =    3.09997052958644174868873332237E+00;
    xtab[24-1] =    3.54444387315534988692540090217E+00;
    xtab[25-1] =    4.00390860386122881522787601332E+00;
    xtab[26-1] =    4.48305535709251834188703761971E+00;
    xtab[27-1] =    4.98891896858994394448649710633E+00;
    xtab[28-1] =    5.53314715156749572511833355558E+00;
    xtab[29-1] =    6.13827922012393462039499237854E+00;
    xtab[30-1] =    6.86334529352989158106110835756E+00;

    weight[ 1-1] =   0.290825470013122622941102747365E-20;
    weight[ 2-1] =   0.281033360275090370876277491534E-16;
    weight[ 3-1] =   0.287860708054870606219239791142E-13;
    weight[ 4-1] =   0.810618629746304420399344796173E-11;
    weight[ 5-1] =   0.917858042437852820850075742492E-09;
    weight[ 6-1] =   0.510852245077594627738963204403E-07;
    weight[ 7-1] =   0.157909488732471028834638794022E-05;
    weight[ 8-1] =   0.293872522892298764150118423412E-04;
    weight[ 9-1] =   0.348310124318685523420995323183E-03;
    weight[10-1] =   0.273792247306765846298942568953E-02;
    weight[11-1] =   0.147038297048266835152773557787E-01;
    weight[12-1] =   0.551441768702342511680754948183E-01;
    weight[13-1] =   0.146735847540890099751693643152E+00;
    weight[14-1] =   0.280130930839212667413493211293E+00;
    weight[15-1] =   0.386394889541813862555601849165E+00;
    weight[16-1] =   0.386394889541813862555601849165E+00;
    weight[17-1] =   0.280130930839212667413493211293E+00;
    weight[18-1] =   0.146735847540890099751693643152E+00;
    weight[19-1] =   0.551441768702342511680754948183E-01;
    weight[20-1] =   0.147038297048266835152773557787E-01;
    weight[21-1] =   0.273792247306765846298942568953E-02;
    weight[22-1] =   0.348310124318685523420995323183E-03;
    weight[23-1] =   0.293872522892298764150118423412E-04;
    weight[24-1] =   0.157909488732471028834638794022E-05;
    weight[25-1] =   0.510852245077594627738963204403E-07;
    weight[26-1] =   0.917858042437852820850075742492E-09;
    weight[27-1] =   0.810618629746304420399344796173E-11;
    weight[28-1] =   0.287860708054870606219239791142E-13;
    weight[29-1] =   0.281033360275090370876277491534E-16;
    weight[30-1] =   0.290825470013122622941102747365E-20;
  }
  else if ( order == 31 )
  {
    xtab[  1-1] =   -6.9956801237185402753248521473232E+00;
    xtab[  2-1] =   -6.2750787049428601427036567812530E+00;
    xtab[  3-1] =   -5.6739614446185883296332558789276E+00;
    xtab[  4-1] =   -5.1335955771123807045862968913996E+00;
    xtab[  5-1] =   -4.6315595063128599420667997654336E+00;
    xtab[  6-1] =   -4.1562717558181451724831352315314E+00;
    xtab[  7-1] =   -3.7007434032314694224497164589673E+00;
    xtab[  8-1] =   -3.2603207323135408104645401509648E+00;
    xtab[  9-1] =   -2.8316804533902054557015640151425E+00;
    xtab[ 10-1] =   -2.4123177054804201051740184582119E+00;
    xtab[ 11-1] =   -2.0002585489356389657975562598571E+00;
    xtab[ 12-1] =   -1.5938858604721398261388419455550E+00;
    xtab[ 13-1] =   -1.1918269983500464260821358649242E+00;
    xtab[ 14-1] =  -0.79287697691530893968593032998830E+00;
    xtab[ 15-1] =  -0.39594273647142311094670041663436E+00;
    xtab[ 16-1] =    0.0000000000000000000000000000000E+00;
    xtab[ 17-1] =   0.39594273647142311094670041663436E+00;
    xtab[ 18-1] =   0.79287697691530893968593032998830E+00;
    xtab[ 19-1] =    1.1918269983500464260821358649242E+00;
    xtab[ 20-1] =    1.5938858604721398261388419455550E+00;
    xtab[ 21-1] =    2.0002585489356389657975562598571E+00;
    xtab[ 22-1] =    2.4123177054804201051740184582119E+00;
    xtab[ 23-1] =    2.8316804533902054557015640151425E+00;
    xtab[ 24-1] =    3.2603207323135408104645401509648E+00;
    xtab[ 25-1] =    3.7007434032314694224497164589673E+00;
    xtab[ 26-1] =    4.1562717558181451724831352315314E+00;
    xtab[ 27-1] =    4.6315595063128599420667997654336E+00;
    xtab[ 28-1] =    5.1335955771123807045862968913996E+00;
    xtab[ 29-1] =    5.6739614446185883296332558789276E+00;
    xtab[ 30-1] =    6.2750787049428601427036567812530E+00;
    xtab[ 31-1] =    6.9956801237185402753248521473232E+00;

    weight[  1-1] =   0.46189683944498305857470556847735E-21;
    weight[  2-1] =   0.51106090079112519643027197715274E-17;
    weight[  3-1] =   0.58995564987355133075257722133966E-14;
    weight[  4-1] =   0.18603735214463569590294465062239E-11;
    weight[  5-1] =   0.23524920032013205739850619940094E-09;
    weight[  6-1] =   0.14611988344865057576066495091513E-07;
    weight[  7-1] =   0.50437125589241034841778074689627E-06;
    weight[  8-1] =   0.10498602757642934202945441341697E-04;
    weight[  9-1] =   0.13952090395003623854995664958146E-03;
    weight[ 10-1] =   0.12336833073030489880608311394968E-02;
    weight[ 11-1] =   0.74827999140119116765002499116934E-02;
    weight[ 12-1] =   0.31847230731201222775249585776902E-01;
    weight[ 13-1] =   0.96717948160569462991143316029341E-01;
    weight[ 14-1] =   0.21213278866810461318136114862419E+00;
    weight[ 15-1] =   0.33877265789305344906000174083214E+00;
    weight[ 16-1] =   0.39577855609737786462923720809676E+00;
    weight[ 17-1] =   0.33877265789305344906000174083214E+00;
    weight[ 18-1] =   0.21213278866810461318136114862419E+00;
    weight[ 19-1] =   0.96717948160569462991143316029341E-01;
    weight[ 20-1] =   0.31847230731201222775249585776902E-01;
    weight[ 21-1] =   0.74827999140119116765002499116934E-02;
    weight[ 22-1] =   0.12336833073030489880608311394968E-02;
    weight[ 23-1] =   0.13952090395003623854995664958146E-03;
    weight[ 24-1] =   0.10498602757642934202945441341697E-04;
    weight[ 25-1] =   0.50437125589241034841778074689627E-06;
    weight[ 26-1] =   0.14611988344865057576066495091513E-07;
    weight[ 27-1] =   0.23524920032013205739850619940094E-09;
    weight[ 28-1] =   0.18603735214463569590294465062239E-11;
    weight[ 29-1] =   0.58995564987355133075257722133966E-14;
    weight[ 30-1] =   0.51106090079112519643027197715274E-17;
    weight[ 31-1] =   0.46189683944498305857470556847735E-21;
  }
  else if ( order == 32 )
  {
    xtab[ 1-1] =   -7.12581390983E+00;
    xtab[ 2-1] =   -6.40949814927E+00;
    xtab[ 3-1] =   -5.81222594952E+00;
    xtab[ 4-1] =   -5.27555098652E+00;
    xtab[ 5-1] =   -4.77716450350E+00;
    xtab[ 6-1] =   -4.30554795335E+00;
    xtab[ 7-1] =   -3.85375548547E+00;
    xtab[ 8-1] =   -3.41716749282E+00;
    xtab[ 9-1] =   -2.99249082500E+00;
    xtab[10-1] =   -2.57724953773E+00;
    xtab[11-1] =   -2.16949918361E+00;
    xtab[12-1] =   -1.76765410946E+00;
    xtab[13-1] =   -1.37037641095E+00;
    xtab[14-1] =  -0.976500463590E+00;
    xtab[15-1] =  -0.584978765436E+00;
    xtab[16-1] =  -0.194840741569E+00;
    xtab[17-1] =   0.194840741569E+00;
    xtab[18-1] =   0.584978765436E+00;
    xtab[19-1] =   0.976500463590E+00;
    xtab[20-1] =    1.37037641095E+00;
    xtab[21-1] =    1.76765410946E+00;
    xtab[22-1] =    2.16949918361E+00;
    xtab[23-1] =    2.57724953773E+00;
    xtab[24-1] =    2.99249082500E+00;
    xtab[25-1] =    3.41716749282E+00;
    xtab[26-1] =    3.85375548547E+00;
    xtab[27-1] =    4.30554795335E+00;
    xtab[28-1] =    4.77716450350E+00;
    xtab[29-1] =    5.27555098652E+00;
    xtab[30-1] =    5.81222594952E+00;
    xtab[31-1] =    6.40949814927E+00;
    xtab[32-1] =    7.12581390983E+00;

    weight[ 1-1] =   0.731067642736E-22;
    weight[ 2-1] =   0.923173653649E-18;
    weight[ 3-1] =   0.119734401709E-14;
    weight[ 4-1] =   0.421501021125E-12;
    weight[ 5-1] =   0.593329146300E-10;
    weight[ 6-1] =   0.409883216476E-08;
    weight[ 7-1] =   0.157416779254E-06;
    weight[ 8-1] =   0.365058512955E-05;
    weight[ 9-1] =   0.541658406172E-04;
    weight[10-1] =   0.536268365526E-03;
    weight[11-1] =   0.365489032664E-02;
    weight[12-1] =   0.175534288315E-01;
    weight[13-1] =   0.604581309557E-01;
    weight[14-1] =   0.151269734076E+00;
    weight[15-1] =   0.277458142302E+00;
    weight[16-1] =   0.375238352592E+00;
    weight[17-1] =   0.375238352592E+00;
    weight[18-1] =   0.277458142302E+00;
    weight[19-1] =   0.151269734076E+00;
    weight[20-1] =   0.604581309557E-01;
    weight[21-1] =   0.175534288315E-01;
    weight[22-1] =   0.365489032664E-02;
    weight[23-1] =   0.536268365526E-03;
    weight[24-1] =   0.541658406172E-04;
    weight[25-1] =   0.365058512955E-05;
    weight[26-1] =   0.157416779254E-06;
    weight[27-1] =   0.409883216476E-08;
    weight[28-1] =   0.593329146300E-10;
    weight[29-1] =   0.421501021125E-12;
    weight[30-1] =   0.119734401709E-14;
    weight[31-1] =   0.923173653649E-18;
    weight[32-1] =   0.731067642736E-22;
  }
  else if ( order == 40 )
  {
    xtab[ 1-1] =   -8.09876113925E+00;
    xtab[ 2-1] =   -7.41158253149E+00;
    xtab[ 3-1] =   -6.84023730525E+00;
    xtab[ 4-1] =   -6.32825535122E+00;
    xtab[ 5-1] =   -5.85409505603E+00;
    xtab[ 6-1] =   -5.40665424797E+00;
    xtab[ 7-1] =   -4.97926097855E+00;
    xtab[ 8-1] =   -4.56750207284E+00;
    xtab[ 9-1] =   -4.16825706683E+00;
    xtab[10-1] =   -3.77920675344E+00;
    xtab[11-1] =   -3.39855826586E+00;
    xtab[12-1] =   -3.02487988390E+00;
    xtab[13-1] =   -2.65699599844E+00;
    xtab[14-1] =   -2.29391714188E+00;
    xtab[15-1] =   -1.93479147228E+00;
    xtab[16-1] =   -1.57886989493E+00;
    xtab[17-1] =   -1.22548010905E+00;
    xtab[18-1] =  -0.874006612357E+00;
    xtab[19-1] =  -0.523874713832E+00;
    xtab[20-1] =  -0.174537214598E+00;
    xtab[21-1] =   0.174537214598E+00;
    xtab[22-1] =   0.523874713832E+00;
    xtab[23-1] =   0.874006612357E+00;
    xtab[24-1] =    1.22548010905E+00;
    xtab[25-1] =    1.57886989493E+00;
    xtab[26-1] =    1.93479147228E+00;
    xtab[27-1] =    2.29391714188E+00;
    xtab[28-1] =    2.65699599844E+00;
    xtab[29-1] =    3.02487988390E+00;
    xtab[30-1] =    3.39855826586E+00;
    xtab[31-1] =    3.77920675344E+00;
    xtab[32-1] =    4.16825706683E+00;
    xtab[33-1] =    4.56750207284E+00;
    xtab[34-1] =    4.97926097855E+00;
    xtab[35-1] =    5.40665424797E+00;
    xtab[36-1] =    5.85409505603E+00;
    xtab[37-1] =    6.32825535122E+00;
    xtab[38-1] =    6.84023730525E+00;
    xtab[39-1] =    7.41158253149E+00;
    xtab[40-1] =    8.09876113925E+00;

    weight[ 1-1] =   0.259104371384E-28;
    weight[ 2-1] =   0.854405696375E-24;
    weight[ 3-1] =   0.256759336540E-20;
    weight[ 4-1] =   0.198918101211E-17;
    weight[ 5-1] =   0.600835878947E-15;
    weight[ 6-1] =   0.880570764518E-13;
    weight[ 7-1] =   0.715652805267E-11;
    weight[ 8-1] =   0.352562079135E-09;
    weight[ 9-1] =   0.112123608322E-07;
    weight[10-1] =   0.241114416359E-06;
    weight[11-1] =   0.363157615067E-05;
    weight[12-1] =   0.393693398108E-04;
    weight[13-1] =   0.313853594540E-03;
    weight[14-1] =   0.187149682959E-02;
    weight[15-1] =   0.846088800823E-02;
    weight[16-1] =   0.293125655361E-01;
    weight[17-1] =   0.784746058652E-01;
    weight[18-1] =   0.163378732713E+00;
    weight[19-1] =   0.265728251876E+00;
    weight[20-1] =   0.338643277425E+00;
    weight[21-1] =   0.338643277425E+00;
    weight[22-1] =   0.265728251876E+00;
    weight[23-1] =   0.163378732713E+00;
    weight[24-1] =   0.784746058652E-01;
    weight[25-1] =   0.293125655361E-01;
    weight[26-1] =   0.846088800823E-02;
    weight[27-1] =   0.187149682959E-02;
    weight[28-1] =   0.313853594540E-03;
    weight[29-1] =   0.393693398108E-04;
    weight[30-1] =   0.363157615067E-05;
    weight[31-1] =   0.241114416359E-06;
    weight[32-1] =   0.112123608322E-07;
    weight[33-1] =   0.352562079135E-09;
    weight[34-1] =   0.715652805267E-11;
    weight[35-1] =   0.880570764518E-13;
    weight[36-1] =   0.600835878947E-15;
    weight[37-1] =   0.198918101211E-17;
    weight[38-1] =   0.256759336540E-20;
    weight[39-1] =   0.854405696375E-24;
    weight[40-1] =   0.259104371384E-28;
  }
  else if ( order == 50 )
  {
    xtab[ 1-1] =   -9.18240695813E+00;
    xtab[ 2-1] =   -8.52277103092E+00;
    xtab[ 3-1] =   -7.97562236821E+00;
    xtab[ 4-1] =   -7.48640942986E+00;
    xtab[ 5-1] =   -7.03432350977E+00;
    xtab[ 6-1] =   -6.60864797386E+00;
    xtab[ 7-1] =   -6.20295251927E+00;
    xtab[ 8-1] =   -5.81299467542E+00;
    xtab[ 9-1] =   -5.43578608722E+00;
    xtab[10-1] =   -5.06911758492E+00;
    xtab[11-1] =   -4.71129366617E+00;
    xtab[12-1] =   -4.36097316045E+00;
    xtab[13-1] =   -4.01706817286E+00;
    xtab[14-1] =   -3.67867706252E+00;
    xtab[15-1] =   -3.34503831394E+00;
    xtab[16-1] =   -3.01549776957E+00;
    xtab[17-1] =   -2.68948470227E+00;
    xtab[18-1] =   -2.36649390430E+00;
    xtab[19-1] =   -2.04607196869E+00;
    xtab[20-1] =   -1.72780654752E+00;
    xtab[21-1] =   -1.41131775490E+00;
    xtab[22-1] =   -1.09625112896E+00;
    xtab[23-1] =  -0.782271729555E+00;
    xtab[24-1] =  -0.469059056678E+00;
    xtab[25-1] =  -0.156302546889E+00;
    xtab[26-1] =   0.156302546889E+00;
    xtab[27-1] =   0.469059056678E+00;
    xtab[28-1] =   0.782271729555E+00;
    xtab[29-1] =    1.09625112896E+00;
    xtab[30-1] =    1.41131775490E+00;
    xtab[31-1] =    1.72780654752E+00;
    xtab[32-1] =    2.04607196869E+00;
    xtab[33-1] =    2.36649390430E+00;
    xtab[34-1] =    2.68948470227E+00;
    xtab[35-1] =    3.01549776957E+00;
    xtab[36-1] =    3.34503831394E+00;
    xtab[37-1] =    3.67867706252E+00;
    xtab[38-1] =    4.01706817286E+00;
    xtab[39-1] =    4.36097316045E+00;
    xtab[40-1] =    4.71129366617E+00;
    xtab[41-1] =    5.06911758492E+00;
    xtab[42-1] =    5.43578608722E+00;
    xtab[43-1] =    5.81299467542E+00;
    xtab[44-1] =    6.20295251927E+00;
    xtab[45-1] =    6.60864797386E+00;
    xtab[46-1] =    7.03432350977E+00;
    xtab[47-1] =    7.48640942986E+00;
    xtab[48-1] =    7.97562236821E+00;
    xtab[49-1] =    8.52277103092E+00;
    xtab[50-1] =    9.18240695813E+00;

    weight[ 1-1] =   0.183379404857E-36;
    weight[ 2-1] =   0.167380166790E-31;
    weight[ 3-1] =   0.121524412340E-27;
    weight[ 4-1] =   0.213765830835E-24;
    weight[ 5-1] =   0.141709359957E-21;
    weight[ 6-1] =   0.447098436530E-19;
    weight[ 7-1] =   0.774238295702E-17;
    weight[ 8-1] =   0.809426189344E-15;
    weight[ 9-1] =   0.546594403180E-13;
    weight[10-1] =   0.250665552389E-11;
    weight[11-1] =   0.811187736448E-10;
    weight[12-1] =   0.190904054379E-08;
    weight[13-1] =   0.334679340401E-07;
    weight[14-1] =   0.445702996680E-06;
    weight[15-1] =   0.458168270794E-05;
    weight[16-1] =   0.368401905377E-04;
    weight[17-1] =   0.234269892109E-03;
    weight[18-1] =   0.118901178175E-02;
    weight[19-1] =   0.485326382616E-02;
    weight[20-1] =   0.160319410684E-01;
    weight[21-1] =   0.430791591566E-01;
    weight[22-1] =   0.945489354768E-01;
    weight[23-1] =   0.170032455676E+00;
    weight[24-1] =   0.251130856331E+00;
    weight[25-1] =   0.305085129203E+00;
    weight[26-1] =   0.305085129203E+00;
    weight[27-1] =   0.251130856331E+00;
    weight[28-1] =   0.170032455676E+00;
    weight[29-1] =   0.945489354768E-01;
    weight[30-1] =   0.430791591566E-01;
    weight[31-1] =   0.160319410684E-01;
    weight[32-1] =   0.485326382616E-02;
    weight[33-1] =   0.118901178175E-02;
    weight[34-1] =   0.234269892109E-03;
    weight[35-1] =   0.368401905377E-04;
    weight[36-1] =   0.458168270794E-05;
    weight[37-1] =   0.445702996680E-06;
    weight[38-1] =   0.334679340401E-07;
    weight[39-1] =   0.190904054379E-08;
    weight[40-1] =   0.811187736448E-10;
    weight[41-1] =   0.250665552389E-11;
    weight[42-1] =   0.546594403180E-13;
    weight[43-1] =   0.809426189344E-15;
    weight[44-1] =   0.774238295702E-17;
    weight[45-1] =   0.447098436530E-19;
    weight[46-1] =   0.141709359957E-21;
    weight[47-1] =   0.213765830835E-24;
    weight[48-1] =   0.121524412340E-27;
    weight[49-1] =   0.167380166790E-31;
    weight[50-1] =   0.183379404857E-36;
  }
  else if ( order == 60 )
  {
    xtab[ 1-1] =   -10.1591092462E+00;
    xtab[ 2-1] =   -9.52090367701E+00;
    xtab[ 3-1] =   -8.99239800140E+00;
    xtab[ 4-1] =   -8.52056928412E+00;
    xtab[ 5-1] =   -8.08518865425E+00;
    xtab[ 6-1] =   -7.67583993750E+00;
    xtab[ 7-1] =   -7.28627659440E+00;
    xtab[ 8-1] =   -6.91238153219E+00;
    xtab[ 9-1] =   -6.55125916706E+00;
    xtab[10-1] =   -6.20077355799E+00;
    xtab[11-1] =   -5.85929019639E+00;
    xtab[12-1] =   -5.52552108614E+00;
    xtab[13-1] =   -5.19842653458E+00;
    xtab[14-1] =   -4.87715007747E+00;
    xtab[15-1] =   -4.56097375794E+00;
    xtab[16-1] =   -4.24928643596E+00;
    xtab[17-1] =   -3.94156073393E+00;
    xtab[18-1] =   -3.63733587617E+00;
    xtab[19-1] =   -3.33620465355E+00;
    xtab[20-1] =   -3.03780333823E+00;
    xtab[21-1] =   -2.74180374807E+00;
    xtab[22-1] =   -2.44790690231E+00;
    xtab[23-1] =   -2.15583787123E+00;
    xtab[24-1] =   -1.86534153123E+00;
    xtab[25-1] =   -1.57617901198E+00;
    xtab[26-1] =   -1.28812467487E+00;
    xtab[27-1] =   -1.00096349956E+00;
    xtab[28-1] =  -0.714488781673E+00;
    xtab[29-1] =  -0.428500064221E+00;
    xtab[30-1] =  -0.142801238703E+00;
    xtab[31-1] =   0.142801238703E+00;
    xtab[32-1] =   0.428500064221E+00;
    xtab[33-1] =   0.714488781673E+00;
    xtab[34-1] =    1.00096349956E+00;
    xtab[35-1] =    1.28812467487E+00;
    xtab[36-1] =    1.57617901198E+00;
    xtab[37-1] =    1.86534153123E+00;
    xtab[38-1] =    2.15583787123E+00;
    xtab[39-1] =    2.44790690231E+00;
    xtab[40-1] =    2.74180374807E+00;
    xtab[41-1] =    3.03780333823E+00;
    xtab[42-1] =    3.33620465355E+00;
    xtab[43-1] =    3.63733587617E+00;
    xtab[44-1] =    3.94156073393E+00;
    xtab[45-1] =    4.24928643596E+00;
    xtab[46-1] =    4.56097375794E+00;
    xtab[47-1] =    4.87715007747E+00;
    xtab[48-1] =    5.19842653458E+00;
    xtab[49-1] =    5.52552108614E+00;
    xtab[50-1] =    5.85929019639E+00;
    xtab[51-1] =    6.20077355799E+00;
    xtab[52-1] =    6.55125916706E+00;
    xtab[53-1] =    6.91238153219E+00;
    xtab[54-1] =    7.28627659440E+00;
    xtab[55-1] =    7.67583993750E+00;
    xtab[56-1] =    8.08518865425E+00;
    xtab[57-1] =    8.52056928412E+00;
    xtab[58-1] =    8.99239800140E+00;
    xtab[59-1] =    9.52090367701E+00;
    xtab[60-1] =    10.1591092462E+00;

    weight[ 1-1] =   0.110958724796E-44;
    weight[ 2-1] =   0.243974758810E-39;
    weight[ 3-1] =   0.377162672698E-35;
    weight[ 4-1] =   0.133255961176E-31;
    weight[ 5-1] =   0.171557314767E-28;
    weight[ 6-1] =   0.102940599693E-25;
    weight[ 7-1] =   0.334575695574E-23;
    weight[ 8-1] =   0.651256725748E-21;
    weight[ 9-1] =   0.815364047300E-19;
    weight[10-1] =   0.692324790956E-17;
    weight[11-1] =   0.415244410968E-15;
    weight[12-1] =   0.181662457614E-13;
    weight[13-1] =   0.594843051597E-12;
    weight[14-1] =   0.148895734905E-10;
    weight[15-1] =   0.289935901280E-09;
    weight[16-1] =   0.445682277521E-08;
    weight[17-1] =   0.547555461926E-07;
    weight[18-1] =   0.543351613419E-06;
    weight[19-1] =   0.439428693625E-05;
    weight[20-1] =   0.291874190415E-04;
    weight[21-1] =   0.160277334681E-03;
    weight[22-1] =   0.731773556963E-03;
    weight[23-1] =   0.279132482894E-02;
    weight[24-1] =   0.893217836028E-02;
    weight[25-1] =   0.240612727660E-01;
    weight[26-1] =   0.547189709320E-01;
    weight[27-1] =   0.105298763697E+00;
    weight[28-1] =   0.171776156918E+00;
    weight[29-1] =   0.237868904958E+00;
    weight[30-1] =   0.279853117522E+00;
    weight[31-1] =   0.279853117522E+00;
    weight[32-1] =   0.237868904958E+00;
    weight[33-1] =   0.171776156918E+00;
    weight[34-1] =   0.105298763697E+00;
    weight[35-1] =   0.547189709320E-01;
    weight[36-1] =   0.240612727660E-01;
    weight[37-1] =   0.893217836028E-02;
    weight[38-1] =   0.279132482894E-02;
    weight[39-1] =   0.731773556963E-03;
    weight[40-1] =   0.160277334681E-03;
    weight[41-1] =   0.291874190415E-04;
    weight[42-1] =   0.439428693625E-05;
    weight[43-1] =   0.543351613419E-06;
    weight[44-1] =   0.547555461926E-07;
    weight[45-1] =   0.445682277521E-08;
    weight[46-1] =   0.289935901280E-09;
    weight[47-1] =   0.148895734905E-10;
    weight[48-1] =   0.594843051597E-12;
    weight[49-1] =   0.181662457614E-13;
    weight[50-1] =   0.415244410968E-15;
    weight[51-1] =   0.692324790956E-17;
    weight[52-1] =   0.815364047300E-19;
    weight[53-1] =   0.651256725748E-21;
    weight[54-1] =   0.334575695574E-23;
    weight[55-1] =   0.102940599693E-25;
    weight[56-1] =   0.171557314767E-28;
    weight[57-1] =   0.133255961176E-31;
    weight[58-1] =   0.377162672698E-35;
    weight[59-1] =   0.243974758810E-39;
    weight[60-1] =   0.110958724796E-44;
  }
  else if ( order == 63 )
  {
    xtab[  1-1] =   -10.435499877854168053468115427285E+00;
    xtab[  2-1] =   -9.8028759912974963635223935286507E+00;
    xtab[  3-1] =   -9.2792019543050391319404745506496E+00;
    xtab[  4-1] =   -8.8118581437284546442526628275570E+00;
    xtab[  5-1] =   -8.3807683451863219343010651043788E+00;
    xtab[  6-1] =   -7.9755950801420373181541806298501E+00;
    xtab[  7-1] =   -7.5901395198641066762479783194468E+00;
    xtab[  8-1] =   -7.2203167078889678461161324222529E+00;
    xtab[  9-1] =   -6.8632544331795368527353285876066E+00;
    xtab[ 10-1] =   -6.5168348106821160605273395854042E+00;
    xtab[ 11-1] =   -6.1794379922705969862418461787263E+00;
    xtab[ 12-1] =   -5.8497884000810673462526582961482E+00;
    xtab[ 13-1] =   -5.5268572526403031425047575122840E+00;
    xtab[ 14-1] =   -5.2097979830408354861575136416263E+00;
    xtab[ 15-1] =   -4.8979018644975742350745099214868E+00;
    xtab[ 16-1] =   -4.5905665744435190229271294569091E+00;
    xtab[ 17-1] =   -4.2872733352824404031727616199454E+00;
    xtab[ 18-1] =   -3.9875699104197157485227052068068E+00;
    xtab[ 19-1] =   -3.6910577000963465117322810559754E+00;
    xtab[ 20-1] =   -3.3973817713303911852755941806287E+00;
    xtab[ 21-1] =   -3.1062230279282566329138616746036E+00;
    xtab[ 22-1] =   -2.8172919672837977750747135657355E+00;
    xtab[ 23-1] =   -2.5303236304712010926855221718499E+00;
    xtab[ 24-1] =   -2.2450734604812066298995918179330E+00;
    xtab[ 25-1] =   -1.9613138583081485293922008411321E+00;
    xtab[ 26-1] =   -1.6788312791720137520802800622638E+00;
    xtab[ 27-1] =   -1.3974237486049625107570752063702E+00;
    xtab[ 28-1] =   -1.1168987050996462690510970277840E+00;
    xtab[ 29-1] =  -0.83707109558947615977737795461293E+00;
    xtab[ 30-1] =  -0.55776166427908221668763665253822E+00;
    xtab[ 31-1] =  -0.27879538567115223986687628627202E+00;
    xtab[ 32-1] =   0.00000000000000000000000000000000E+00;
    xtab[ 33-1] =   0.27879538567115223986687628627202E+00;
    xtab[ 34-1] =   0.55776166427908221668763665253822E+00;
    xtab[ 35-1] =   0.83707109558947615977737795461293E+00;
    xtab[ 36-1] =    1.1168987050996462690510970277840E+00;
    xtab[ 37-1] =    1.3974237486049625107570752063702E+00;
    xtab[ 38-1] =    1.6788312791720137520802800622638E+00;
    xtab[ 39-1] =    1.9613138583081485293922008411321E+00;
    xtab[ 40-1] =    2.2450734604812066298995918179330E+00;
    xtab[ 41-1] =    2.5303236304712010926855221718499E+00;
    xtab[ 42-1] =    2.8172919672837977750747135657355E+00;
    xtab[ 43-1] =    3.1062230279282566329138616746036E+00;
    xtab[ 44-1] =    3.3973817713303911852755941806287E+00;
    xtab[ 45-1] =    3.6910577000963465117322810559754E+00;
    xtab[ 46-1] =    3.9875699104197157485227052068068E+00;
    xtab[ 47-1] =    4.2872733352824404031727616199454E+00;
    xtab[ 48-1] =    4.5905665744435190229271294569091E+00;
    xtab[ 49-1] =    4.8979018644975742350745099214868E+00;
    xtab[ 50-1] =    5.2097979830408354861575136416263E+00;
    xtab[ 51-1] =    5.5268572526403031425047575122840E+00;
    xtab[ 52-1] =    5.8497884000810673462526582961482E+00;
    xtab[ 53-1] =    6.1794379922705969862418461787263E+00;
    xtab[ 54-1] =    6.5168348106821160605273395854042E+00;
    xtab[ 55-1] =    6.8632544331795368527353285876066E+00;
    xtab[ 56-1] =    7.2203167078889678461161324222529E+00;
    xtab[ 57-1] =    7.5901395198641066762479783194468E+00;
    xtab[ 58-1] =    7.9755950801420373181541806298501E+00;
    xtab[ 59-1] =    8.3807683451863219343010651043788E+00;
    xtab[ 60-1] =    8.8118581437284546442526628275570E+00;
    xtab[ 61-1] =    9.2792019543050391319404745506496E+00;
    xtab[ 62-1] =    9.8028759912974963635223935286507E+00;
    xtab[ 63-1] =    10.435499877854168053468115427285E+00;

    weight[  1-1] =   0.37099206434787551197827130470031E-47;
    weight[  2-1] =   0.10400778615192299534481914814892E-41;
    weight[  3-1] =   0.19796804708258311251124226474396E-37;
    weight[  4-1] =   0.84687478191640015120141181138947E-34;
    weight[  5-1] =   0.13071305930779945903630127634063E-30;
    weight[  6-1] =   0.93437837175367456929765381518998E-28;
    weight[  7-1] =   0.36027426635173044862245783257252E-25;
    weight[  8-1] =   0.82963863115951789374753323156164E-23;
    weight[  9-1] =   0.12266629909105281472971700203949E-20;
    weight[ 10-1] =   0.12288435628797061539461585325494E-18;
    weight[ 11-1] =   0.86925536958188009075932426691516E-17;
    weight[ 12-1] =   0.44857058689176221240330804981619E-15;
    weight[ 13-1] =   0.17335817955735154599902643794700E-13;
    weight[ 14-1] =   0.51265062385038307838565047455223E-12;
    weight[ 15-1] =   0.11808921844532942490513037158404E-10;
    weight[ 16-1] =   0.21508698297808025739828859845140E-09;
    weight[ 17-1] =   0.31371929535285447801497640621672E-08;
    weight[ 18-1] =   0.37041625984781705796752840204084E-07;
    weight[ 19-1] =   0.35734732949879669663960738150956E-06;
    weight[ 20-1] =   0.28393114498380927832990899215541E-05;
    weight[ 21-1] =   0.18709113003730498008961134765721E-04;
    weight[ 22-1] =   0.10284880800653635546698378640623E-03;
    weight[ 23-1] =   0.47411702610173128107201781718693E-03;
    weight[ 24-1] =   0.18409222622384813438539657470055E-02;
    weight[ 25-1] =   0.60436044551187631655712178246467E-02;
    weight[ 26-1] =   0.16829299199599730926458559757600E-01;
    weight[ 27-1] =   0.39858264027692992170237391875317E-01;
    weight[ 28-1] =   0.80467087993950415219587554532823E-01;
    weight[ 29-1] =   0.13871950817615293377792092082674E+00;
    weight[ 30-1] =   0.20448695346833761570957197160475E+00;
    weight[ 31-1] =   0.25799889943058042204920467417642E+00;
    weight[ 32-1] =   0.27876694884838411919175686949858E+00;
    weight[ 33-1] =   0.25799889943058042204920467417642E+00;
    weight[ 34-1] =   0.20448695346833761570957197160475E+00;
    weight[ 35-1] =   0.13871950817615293377792092082674E+00;
    weight[ 36-1] =   0.80467087993950415219587554532823E-01;
    weight[ 37-1] =   0.39858264027692992170237391875317E-01;
    weight[ 38-1] =   0.16829299199599730926458559757600E-01;
    weight[ 39-1] =   0.60436044551187631655712178246467E-02;
    weight[ 40-1] =   0.18409222622384813438539657470055E-02;
    weight[ 41-1] =   0.47411702610173128107201781718693E-03;
    weight[ 42-1] =   0.10284880800653635546698378640623E-03;
    weight[ 43-1] =   0.18709113003730498008961134765721E-04;
    weight[ 44-1] =   0.28393114498380927832990899215541E-05;
    weight[ 45-1] =   0.35734732949879669663960738150956E-06;
    weight[ 46-1] =   0.37041625984781705796752840204084E-07;
    weight[ 47-1] =   0.31371929535285447801497640621672E-08;
    weight[ 48-1] =   0.21508698297808025739828859845140E-09;
    weight[ 49-1] =   0.11808921844532942490513037158404E-10;
    weight[ 50-1] =   0.51265062385038307838565047455223E-12;
    weight[ 51-1] =   0.17335817955735154599902643794700E-13;
    weight[ 52-1] =   0.44857058689176221240330804981619E-15;
    weight[ 53-1] =   0.86925536958188009075932426691516E-17;
    weight[ 54-1] =   0.12288435628797061539461585325494E-18;
    weight[ 55-1] =   0.12266629909105281472971700203949E-20;
    weight[ 56-1] =   0.82963863115951789374753323156164E-23;
    weight[ 57-1] =   0.36027426635173044862245783257252E-25;
    weight[ 58-1] =   0.93437837175367456929765381518998E-28;
    weight[ 59-1] =   0.13071305930779945903630127634063E-30;
    weight[ 60-1] =   0.84687478191640015120141181138947E-34;
    weight[ 61-1] =   0.19796804708258311251124226474396E-37;
    weight[ 62-1] =   0.10400778615192299534481914814892E-41;
    weight[ 63-1] =   0.37099206434787551197827130470031E-47;
  }
  else if ( order == 64 )
  {
    xtab[ 1-1] =   -10.5261231680E+00;
    xtab[ 2-1] =   -9.89528758683E+00;
    xtab[ 3-1] =   -9.37315954965E+00;
    xtab[ 4-1] =   -8.90724909996E+00;
    xtab[ 5-1] =   -8.47752908338E+00;
    xtab[ 6-1] =   -8.07368728501E+00;
    xtab[ 7-1] =   -7.68954016404E+00;
    xtab[ 8-1] =   -7.32101303278E+00;
    xtab[ 9-1] =   -6.96524112055E+00;
    xtab[10-1] =   -6.62011226264E+00;
    xtab[11-1] =   -6.28401122877E+00;
    xtab[12-1] =   -5.95566632680E+00;
    xtab[13-1] =   -5.63405216435E+00;
    xtab[14-1] =   -5.31832522463E+00;
    xtab[15-1] =   -5.00777960220E+00;
    xtab[16-1] =   -4.70181564741E+00;
    xtab[17-1] =   -4.39991716823E+00;
    xtab[18-1] =   -4.10163447457E+00;
    xtab[19-1] =   -3.80657151395E+00;
    xtab[20-1] =   -3.51437593574E+00;
    xtab[21-1] =   -3.22473129199E+00;
    xtab[22-1] =   -2.93735082300E+00;
    xtab[23-1] =   -2.65197243543E+00;
    xtab[24-1] =   -2.36835458863E+00;
    xtab[25-1] =   -2.08627287988E+00;
    xtab[26-1] =   -1.80551717147E+00;
    xtab[27-1] =   -1.52588914021E+00;
    xtab[28-1] =   -1.24720015694E+00;
    xtab[29-1] =  -0.969269423071E+00;
    xtab[30-1] =  -0.691922305810E+00;
    xtab[31-1] =  -0.414988824121E+00;
    xtab[32-1] =  -0.138302244987E+00;
    xtab[33-1] =   0.138302244987E+00;
    xtab[34-1] =   0.414988824121E+00;
    xtab[35-1] =   0.691922305810E+00;
    xtab[36-1] =   0.969269423071E+00;
    xtab[37-1] =    1.24720015694E+00;
    xtab[38-1] =    1.52588914021E+00;
    xtab[39-1] =    1.80551717147E+00;
    xtab[40-1] =    2.08627287988E+00;
    xtab[41-1] =    2.36835458863E+00;
    xtab[42-1] =    2.65197243543E+00;
    xtab[43-1] =    2.93735082300E+00;
    xtab[44-1] =    3.22473129199E+00;
    xtab[45-1] =    3.51437593574E+00;
    xtab[46-1] =    3.80657151395E+00;
    xtab[47-1] =    4.10163447457E+00;
    xtab[48-1] =    4.39991716823E+00;
    xtab[49-1] =    4.70181564741E+00;
    xtab[50-1] =    5.00777960220E+00;
    xtab[51-1] =    5.31832522463E+00;
    xtab[52-1] =    5.63405216435E+00;
    xtab[53-1] =    5.95566632680E+00;
    xtab[54-1] =    6.28401122877E+00;
    xtab[55-1] =    6.62011226264E+00;
    xtab[56-1] =    6.96524112055E+00;
    xtab[57-1] =    7.32101303278E+00;
    xtab[58-1] =    7.68954016404E+00;
    xtab[59-1] =    8.07368728501E+00;
    xtab[60-1] =    8.47752908338E+00;
    xtab[61-1] =    8.90724909996E+00;
    xtab[62-1] =    9.37315954965E+00;
    xtab[63-1] =    9.89528758683E+00;
    xtab[64-1] =    10.5261231680E+00;

    weight[ 1-1] =   0.553570653584E-48;
    weight[ 2-1] =   0.167974799010E-42;
    weight[ 3-1] =   0.342113801099E-38;
    weight[ 4-1] =   0.155739062462E-34;
    weight[ 5-1] =   0.254966089910E-31;
    weight[ 6-1] =   0.192910359546E-28;
    weight[ 7-1] =   0.786179778889E-26;
    weight[ 8-1] =   0.191170688329E-23;
    weight[ 9-1] =   0.298286278427E-21;
    weight[10-1] =   0.315225456649E-19;
    weight[11-1] =   0.235188471067E-17;
    weight[12-1] =   0.128009339117E-15;
    weight[13-1] =   0.521862372645E-14;
    weight[14-1] =   0.162834073070E-12;
    weight[15-1] =   0.395917776693E-11;
    weight[16-1] =   0.761521725012E-10;
    weight[17-1] =   0.117361674232E-08;
    weight[18-1] =   0.146512531647E-07;
    weight[19-1] =   0.149553293672E-06;
    weight[20-1] =   0.125834025103E-05;
    weight[21-1] =   0.878849923082E-05;
    weight[22-1] =   0.512592913577E-04;
    weight[23-1] =   0.250983698512E-03;
    weight[24-1] =   0.103632909950E-02;
    weight[25-1] =   0.362258697852E-02;
    weight[26-1] =   0.107560405098E-01;
    weight[27-1] =   0.272031289536E-01;
    weight[28-1] =   0.587399819634E-01;
    weight[29-1] =   0.108498349306E+00;
    weight[30-1] =   0.171685842349E+00;
    weight[31-1] =   0.232994786062E+00;
    weight[32-1] =   0.271377424940E+00;
    weight[33-1] =   0.271377424940E+00;
    weight[34-1] =   0.232994786062E+00;
    weight[35-1] =   0.171685842349E+00;
    weight[36-1] =   0.108498349306E+00;
    weight[37-1] =   0.587399819634E-01;
    weight[38-1] =   0.272031289536E-01;
    weight[39-1] =   0.107560405098E-01;
    weight[40-1] =   0.362258697852E-02;
    weight[41-1] =   0.103632909950E-02;
    weight[42-1] =   0.250983698512E-03;
    weight[43-1] =   0.512592913577E-04;
    weight[44-1] =   0.878849923082E-05;
    weight[45-1] =   0.125834025103E-05;
    weight[46-1] =   0.149553293672E-06;
    weight[47-1] =   0.146512531647E-07;
    weight[48-1] =   0.117361674232E-08;
    weight[49-1] =   0.761521725012E-10;
    weight[50-1] =   0.395917776693E-11;
    weight[51-1] =   0.162834073070E-12;
    weight[52-1] =   0.521862372645E-14;
    weight[53-1] =   0.128009339117E-15;
    weight[54-1] =   0.235188471067E-17;
    weight[55-1] =   0.315225456649E-19;
    weight[56-1] =   0.298286278427E-21;
    weight[57-1] =   0.191170688329E-23;
    weight[58-1] =   0.786179778889E-26;
    weight[59-1] =   0.192910359546E-28;
    weight[60-1] =   0.254966089910E-31;
    weight[61-1] =   0.155739062462E-34;
    weight[62-1] =   0.342113801099E-38;
    weight[63-1] =   0.167974799010E-42;
    weight[64-1] =   0.553570653584E-48;
  }
  else if ( order == 127 )
  {
    xtab[  1-1] =   -15.228338148167350978246954433464E+00;
    xtab[  2-1] =   -14.669595158833972632746354112896E+00;
    xtab[  3-1] =   -14.209085995284870755168244250887E+00;
    xtab[  4-1] =   -13.799722290211676634645246746673E+00;
    xtab[  5-1] =   -13.423518590070950062438258321855E+00;
    xtab[  6-1] =   -13.071208660474601901583995439649E+00;
    xtab[  7-1] =   -12.737235652415686338138003924072E+00;
    xtab[  8-1] =   -12.417939378869715805445879624069E+00;
    xtab[  9-1] =   -12.110749020947747600132123508132E+00;
    xtab[ 10-1] =   -11.813772198267727195134584136191E+00;
    xtab[ 11-1] =   -11.525565112572696599167888588564E+00;
    xtab[ 12-1] =   -11.244994583785543445194384194300E+00;
    xtab[ 13-1] =   -10.971150569840247423423040263881E+00;
    xtab[ 14-1] =   -10.703288201027481347670940744690E+00;
    xtab[ 15-1] =   -10.440787957772772867742591798027E+00;
    xtab[ 16-1] =   -10.183127473450343888624126450357E+00;
    xtab[ 17-1] =   -9.9298610495114250736847004273684E+00;
    xtab[ 18-1] =   -9.6806044412474728038150712732737E+00;
    xtab[ 19-1] =   -9.4350233389881650135019598506287E+00;
    xtab[ 20-1] =   -9.1928244988460305715774195052527E+00;
    xtab[ 21-1] =   -8.9537488108565404323807890169970E+00;
    xtab[ 22-1] =   -8.7175658087076307363833999548548E+00;
    xtab[ 23-1] =   -8.4840692689832473326097180339984E+00;
    xtab[ 24-1] =   -8.2530736454457156579694124243888E+00;
    xtab[ 25-1] =   -8.0244111514703375578594739796798E+00;
    xtab[ 26-1] =   -7.7979293513870105420829120455591E+00;
    xtab[ 27-1] =   -7.5734891556083454022834960763301E+00;
    xtab[ 28-1] =   -7.3509631392269052701961258043733E+00;
    xtab[ 29-1] =   -7.1302341220350710668064025713431E+00;
    xtab[ 30-1] =   -6.9111939615465713197465633109366E+00;
    xtab[ 31-1] =   -6.6937425208758294190074417381666E+00;
    xtab[ 32-1] =   -6.4777867811645365448144903821487E+00;
    xtab[ 33-1] =   -6.2632400742737354345609723857092E+00;
    xtab[ 34-1] =   -6.0500214161419845694465474482388E+00;
    xtab[ 35-1] =   -5.8380549248774187386601690807757E+00;
    xtab[ 36-1] =   -5.6272693105464816659423455794909E+00;
    xtab[ 37-1] =   -5.4175974259243240722848425872924E+00;
    xtab[ 38-1] =   -5.2089758693153983587570258372239E+00;
    xtab[ 39-1] =   -5.0013446320386360038520809107373E+00;
    xtab[ 40-1] =   -4.7946467843764925009748509930857E+00;
    xtab[ 41-1] =   -4.5888281947698372951606485031212E+00;
    xtab[ 42-1] =   -4.3838372778464736294253744407459E+00;
    xtab[ 43-1] =   -4.1796247675352031349421189892408E+00;
    xtab[ 44-1] =   -3.9761435120673355916035814195920E+00;
    xtab[ 45-1] =   -3.7733482881250526721004678400057E+00;
    xtab[ 46-1] =   -3.5711956317782180447199756485249E+00;
    xtab[ 47-1] =   -3.3696436841717397896643629240035E+00;
    xtab[ 48-1] =   -3.1686520501953630191857798261495E+00;
    xtab[ 49-1] =   -2.9681816685955910267761649521505E+00;
    xtab[ 50-1] =   -2.7681946921824058801226545958892E+00;
    xtab[ 51-1] =   -2.5686543769473501723144013022363E+00;
    xtab[ 52-1] =   -2.3695249790490401080012474645702E+00;
    xtab[ 53-1] =   -2.1707716587411506879498498083695E+00;
    xtab[ 54-1] =   -1.9723603904195020079324743227565E+00;
    xtab[ 55-1] =   -1.7742578780516791584676442103681E+00;
    xtab[ 56-1] =   -1.5764314753267801315519597621879E+00;
    xtab[ 57-1] =   -1.3788491099261778091441557053728E+00;
    xtab[ 58-1] =   -1.1814792113700685848678583598423E+00;
    xtab[ 59-1] =  -0.98429064194027277726568984213773E+00;
    xtab[ 60-1] =  -0.78725263021825034151596831878971E+00;
    xtab[ 61-1] =  -0.59033470680942102142230439346102E+00;
    xtab[ 62-1] =  -0.39350664185130136568037826200185E+00;
    xtab[ 63-1] =  -0.19673838392423251964272239737078E+00;
    xtab[ 64-1] =    0.0000000000000000000000000000000E+00;
    xtab[ 65-1] =   0.19673838392423251964272239737078E+00;
    xtab[ 66-1] =   0.39350664185130136568037826200185E+00;
    xtab[ 67-1] =   0.59033470680942102142230439346102E+00;
    xtab[ 68-1] =   0.78725263021825034151596831878971E+00;
    xtab[ 69-1] =   0.98429064194027277726568984213773E+00;
    xtab[ 70-1] =    1.1814792113700685848678583598423E+00;
    xtab[ 71-1] =    1.3788491099261778091441557053728E+00;
    xtab[ 72-1] =    1.5764314753267801315519597621879E+00;
    xtab[ 73-1] =    1.7742578780516791584676442103681E+00;
    xtab[ 74-1] =    1.9723603904195020079324743227565E+00;
    xtab[ 75-1] =    2.1707716587411506879498498083695E+00;
    xtab[ 76-1] =    2.3695249790490401080012474645702E+00;
    xtab[ 77-1] =    2.5686543769473501723144013022363E+00;
    xtab[ 78-1] =    2.7681946921824058801226545958892E+00;
    xtab[ 79-1] =    2.9681816685955910267761649521505E+00;
    xtab[ 80-1] =    3.1686520501953630191857798261495E+00;
    xtab[ 81-1] =    3.3696436841717397896643629240035E+00;
    xtab[ 82-1] =    3.5711956317782180447199756485249E+00;
    xtab[ 83-1] =    3.7733482881250526721004678400057E+00;
    xtab[ 84-1] =    3.9761435120673355916035814195920E+00;
    xtab[ 85-1] =    4.1796247675352031349421189892408E+00;
    xtab[ 86-1] =    4.3838372778464736294253744407459E+00;
    xtab[ 87-1] =    4.5888281947698372951606485031212E+00;
    xtab[ 88-1] =    4.7946467843764925009748509930857E+00;
    xtab[ 89-1] =    5.0013446320386360038520809107373E+00;
    xtab[ 90-1] =    5.2089758693153983587570258372239E+00;
    xtab[ 91-1] =    5.4175974259243240722848425872924E+00;
    xtab[ 92-1] =    5.6272693105464816659423455794909E+00;
    xtab[ 93-1] =    5.8380549248774187386601690807757E+00;
    xtab[ 94-1] =    6.0500214161419845694465474482388E+00;
    xtab[ 95-1] =    6.2632400742737354345609723857092E+00;
    xtab[ 96-1] =    6.4777867811645365448144903821487E+00;
    xtab[ 97-1] =    6.6937425208758294190074417381666E+00;
    xtab[ 98-1] =    6.9111939615465713197465633109366E+00;
    xtab[ 99-1] =    7.1302341220350710668064025713431E+00;
    xtab[100-1] =    7.3509631392269052701961258043733E+00;
    xtab[101-1] =    7.5734891556083454022834960763301E+00;
    xtab[102-1] =    7.7979293513870105420829120455591E+00;
    xtab[103-1] =    8.0244111514703375578594739796798E+00;
    xtab[104-1] =    8.2530736454457156579694124243888E+00;
    xtab[105-1] =    8.4840692689832473326097180339984E+00;
    xtab[106-1] =    8.7175658087076307363833999548548E+00;
    xtab[107-1] =    8.9537488108565404323807890169970E+00;
    xtab[108-1] =    9.1928244988460305715774195052527E+00;
    xtab[109-1] =    9.4350233389881650135019598506287E+00;
    xtab[110-1] =    9.6806044412474728038150712732737E+00;
    xtab[111-1] =    9.9298610495114250736847004273684E+00;
    xtab[112-1] =    10.183127473450343888624126450357E+00;
    xtab[113-1] =    10.440787957772772867742591798027E+00;
    xtab[114-1] =    10.703288201027481347670940744690E+00;
    xtab[115-1] =    10.971150569840247423423040263881E+00;
    xtab[116-1] =    11.244994583785543445194384194300E+00;
    xtab[117-1] =    11.525565112572696599167888588564E+00;
    xtab[118-1] =    11.813772198267727195134584136191E+00;
    xtab[119-1] =    12.110749020947747600132123508132E+00;
    xtab[120-1] =    12.417939378869715805445879624069E+00;
    xtab[121-1] =    12.737235652415686338138003924072E+00;
    xtab[122-1] =    13.071208660474601901583995439649E+00;
    xtab[123-1] =    13.423518590070950062438258321855E+00;
    xtab[124-1] =    13.799722290211676634645246746673E+00;
    xtab[125-1] =    14.209085995284870755168244250887E+00;
    xtab[126-1] =    14.669595158833972632746354112896E+00;
    xtab[127-1] =    15.228338148167350978246954433464E+00;

    weight[  1-1] =   0.12504497577050595552677230002883E-100;
    weight[  2-1] =   0.17272798059419131415318615789672E-93;
    weight[  3-1] =   0.89321681571986548608031150791499E-88;
    weight[  4-1] =   0.77306185240893578449625186483810E-83;
    weight[  5-1] =   0.20143957652648255497735460506196E-78;
    weight[  6-1] =   0.21503714733610239701351039429345E-74;
    weight[  7-1] =   0.11341924208594594813715533569504E-70;
    weight[  8-1] =   0.33489139011795051950683388483136E-67;
    weight[  9-1] =   0.60486548964016681064424451668405E-64;
    weight[ 10-1] =   0.71375092946352177824971347343892E-61;
    weight[ 11-1] =   0.57884563374885556636801095624030E-58;
    weight[ 12-1] =   0.33581166223858230300409326551248E-55;
    weight[ 13-1] =   0.14394641949253923568603163698953E-52;
    weight[ 14-1] =   0.46821808383216117724080263903889E-50;
    weight[ 15-1] =   0.11817054440684264071348471955361E-47;
    weight[ 16-1] =   0.23581659156008927203181682045005E-45;
    weight[ 17-1] =   0.37814427940797540210712758405540E-43;
    weight[ 18-1] =   0.49411031115771638145610738414006E-41;
    weight[ 19-1] =   0.53255303775425059266087298458297E-39;
    weight[ 20-1] =   0.47854390680131484999315199332765E-37;
    weight[ 21-1] =   0.36191883445952356128627543209554E-35;
    weight[ 22-1] =   0.23232083386343554805352497446119E-33;
    weight[ 23-1] =   0.12753331411008716683688974281454E-31;
    weight[ 24-1] =   0.60277753850758742112436095241270E-30;
    weight[ 25-1] =   0.24679773241777200207460855084439E-28;
    weight[ 26-1] =   0.88019567691698482573264198727415E-27;
    weight[ 27-1] =   0.27482489212040561315005725890593E-25;
    weight[ 28-1] =   0.75468218903085486125222816438456E-24;
    weight[ 29-1] =   0.18303134636280466270545996891835E-22;
    weight[ 30-1] =   0.39355990860860813085582448449811E-21;
    weight[ 31-1] =   0.75293161638581191068419292570042E-20;
    weight[ 32-1] =   0.12857997786722855037584105682618E-18;
    weight[ 33-1] =   0.19659326888445857792541925311450E-17;
    weight[ 34-1] =   0.26986511907214101894995783364250E-16;
    weight[ 35-1] =   0.33344414303198856330118301113874E-15;
    weight[ 36-1] =   0.37173303125150639885726463109574E-14;
    weight[ 37-1] =   0.37473954472839737091885387788983E-13;
    weight[ 38-1] =   0.34230094493397259538669512076007E-12;
    weight[ 39-1] =   0.28385303724993373166810860630552E-11;
    weight[ 40-1] =   0.21406920290454669208938772802828E-10;
    weight[ 41-1] =   0.14706331273431716244229273183839E-09;
    weight[ 42-1] =   0.92173940967434659264335883218167E-09;
    weight[ 43-1] =   0.52781663936972714041837056042506E-08;
    weight[ 44-1] =   0.27650497044951117835905283127679E-07;
    weight[ 45-1] =   0.13267855842539464770913063113371E-06;
    weight[ 46-1] =   0.58380944276113062188573331195042E-06;
    weight[ 47-1] =   0.23581561724775629112332165335800E-05;
    weight[ 48-1] =   0.87524468034280444703919485644809E-05;
    weight[ 49-1] =   0.29876790535909012274846532159647E-04;
    weight[ 50-1] =   0.93874435720072545206729594267039E-04;
    weight[ 51-1] =   0.27170762627931172053444716883938E-03;
    weight[ 52-1] =   0.72493929742498358979684249380921E-03;
    weight[ 53-1] =   0.17841208326763432884316727108264E-02;
    weight[ 54-1] =   0.40524855186046131499765636276283E-02;
    weight[ 55-1] =   0.85000263041544110385806705526917E-02;
    weight[ 56-1] =   0.16471142241609687824005585301760E-01;
    weight[ 57-1] =   0.29499296248213632269675010319119E-01;
    weight[ 58-1] =   0.48847387114300011006959603975676E-01;
    weight[ 59-1] =   0.74807989768583731416517226905270E-01;
    weight[ 60-1] =   0.10598520508090929403834368934301E+00;
    weight[ 61-1] =   0.13893945309051540832066283010510E+00;
    weight[ 62-1] =   0.16856236074207929740526975049765E+00;
    weight[ 63-1] =   0.18927849580120432177170145550076E+00;
    weight[ 64-1] =   0.19673340688823289786163676995151E+00;
    weight[ 65-1] =   0.18927849580120432177170145550076E+00;
    weight[ 66-1] =   0.16856236074207929740526975049765E+00;
    weight[ 67-1] =   0.13893945309051540832066283010510E+00;
    weight[ 68-1] =   0.10598520508090929403834368934301E+00;
    weight[ 69-1] =   0.74807989768583731416517226905270E-01;
    weight[ 70-1] =   0.48847387114300011006959603975676E-01;
    weight[ 71-1] =   0.29499296248213632269675010319119E-01;
    weight[ 72-1] =   0.16471142241609687824005585301760E-01;
    weight[ 73-1] =   0.85000263041544110385806705526917E-02;
    weight[ 74-1] =   0.40524855186046131499765636276283E-02;
    weight[ 75-1] =   0.17841208326763432884316727108264E-02;
    weight[ 76-1] =   0.72493929742498358979684249380921E-03;
    weight[ 77-1] =   0.27170762627931172053444716883938E-03;
    weight[ 78-1] =   0.93874435720072545206729594267039E-04;
    weight[ 79-1] =   0.29876790535909012274846532159647E-04;
    weight[ 80-1] =   0.87524468034280444703919485644809E-05;
    weight[ 81-1] =   0.23581561724775629112332165335800E-05;
    weight[ 82-1] =   0.58380944276113062188573331195042E-06;
    weight[ 83-1] =   0.13267855842539464770913063113371E-06;
    weight[ 84-1] =   0.27650497044951117835905283127679E-07;
    weight[ 85-1] =   0.52781663936972714041837056042506E-08;
    weight[ 86-1] =   0.92173940967434659264335883218167E-09;
    weight[ 87-1] =   0.14706331273431716244229273183839E-09;
    weight[ 88-1] =   0.21406920290454669208938772802828E-10;
    weight[ 89-1] =   0.28385303724993373166810860630552E-11;
    weight[ 90-1] =   0.34230094493397259538669512076007E-12;
    weight[ 91-1] =   0.37473954472839737091885387788983E-13;
    weight[ 92-1] =   0.37173303125150639885726463109574E-14;
    weight[ 93-1] =   0.33344414303198856330118301113874E-15;
    weight[ 94-1] =   0.26986511907214101894995783364250E-16;
    weight[ 95-1] =   0.19659326888445857792541925311450E-17;
    weight[ 96-1] =   0.12857997786722855037584105682618E-18;
    weight[ 97-1] =   0.75293161638581191068419292570042E-20;
    weight[ 98-1] =   0.39355990860860813085582448449811E-21;
    weight[ 99-1] =   0.18303134636280466270545996891835E-22;
    weight[100-1] =   0.75468218903085486125222816438456E-24;
    weight[101-1] =   0.27482489212040561315005725890593E-25;
    weight[102-1] =   0.88019567691698482573264198727415E-27;
    weight[103-1] =   0.24679773241777200207460855084439E-28;
    weight[104-1] =   0.60277753850758742112436095241270E-30;
    weight[105-1] =   0.12753331411008716683688974281454E-31;
    weight[106-1] =   0.23232083386343554805352497446119E-33;
    weight[107-1] =   0.36191883445952356128627543209554E-35;
    weight[108-1] =   0.47854390680131484999315199332765E-37;
    weight[109-1] =   0.53255303775425059266087298458297E-39;
    weight[110-1] =   0.49411031115771638145610738414006E-41;
    weight[111-1] =   0.37814427940797540210712758405540E-43;
    weight[112-1] =   0.23581659156008927203181682045005E-45;
    weight[113-1] =   0.11817054440684264071348471955361E-47;
    weight[114-1] =   0.46821808383216117724080263903889E-50;
    weight[115-1] =   0.14394641949253923568603163698953E-52;
    weight[116-1] =   0.33581166223858230300409326551248E-55;
    weight[117-1] =   0.57884563374885556636801095624030E-58;
    weight[118-1] =   0.71375092946352177824971347343892E-61;
    weight[119-1] =   0.60486548964016681064424451668405E-64;
    weight[120-1] =   0.33489139011795051950683388483136E-67;
    weight[121-1] =   0.11341924208594594813715533569504E-70;
    weight[122-1] =   0.21503714733610239701351039429345E-74;
    weight[123-1] =   0.20143957652648255497735460506196E-78;
    weight[124-1] =   0.77306185240893578449625186483810E-83;
    weight[125-1] =   0.89321681571986548608031150791499E-88;
    weight[126-1] =   0.17272798059419131415318615789672E-93;
    weight[127-1] =   0.12504497577050595552677230002883E-100;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "HERMITE_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 to 20,\n" );
    fprintf ( stderr, "  30, 31, 32, 40, 50, 60, 63, 64 or 127.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

int i4_factorial2 ( int n )

/******************************************************************************/
/*
  Purpose:
  
    I4_FACTORIAL2 computes the double factorial function N!!
  
  Formula:
  
    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
  
  Example:
  
     N    N!!
  
     0     1
     1     1
     2     2
     3     3
     4     8
     5    15
     6    48
     7   105
     8   384
     9   945
    10  3840
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    19 February 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the argument of the double factorial function.
    If N is less than 1, I4_FACTORIAL2 is returned as 1.
  
    Output, int I4_FACTORIAL2, the value of N!!.
*/
{
  int value;

  if ( n < 1 )
  {
    return 1;
  }

  value = 1;

  while ( 1 < n )
  {
    value = value * n;
    n = n - 2;
  }

  return value;
}
/******************************************************************************/

int i4_min ( int i1, int i2 )

/******************************************************************************/
/*
  Purpose:
  
    I4_MIN returns the smaller of two I4's.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    13 October 1998
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int I1, I2, two integers to be compared.
  
    Output, int I4_MIN, the smaller of I1 and I2.
  
*/
{
  if ( i1 < i2 )
  {
    return i1;
  }
  else
  {
    return i2;
  }

}
/******************************************************************************/

int i4_power ( int i, int j )

/******************************************************************************/
/*
  Purpose:
  
    I4_POWER returns the value of I^J.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    01 April 2004
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int I, J, the base and the power.  J should be nonnegative.
  
    Output, int I4_POWER, the value of I^J.
*/
{
  int k;
  int value;

  if ( j < 0 )
  {
    if ( i == 1 )
    {
      value = 1;
    }
    else if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J negative.\n" );
      exit ( 1 );
    }
    else
    {
      value = 0;
    }
  }
  else if ( j == 0 )
  {
    if ( i == 0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "I4_POWER - Fatal error!\n" );
      fprintf ( stderr, "  I^J requested, with I = 0 and J = 0.\n" );
      exit ( 1 );
    }
    else
    {
      value = 1;
    }
  }
  else if ( j == 1 )
  {
    value = i;
  }
  else
  {
    value = 1;
    for ( k = 1; k <= j; k++ )
    {
      value = value * i;
    }
  }
  return value;
}
/******************************************************************************/

void jacobi_compute ( int order, double alpha, double beta, double xtab[],
  double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    JACOBI_COMPUTE computes a Gauss-Jacobi quadrature rule.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x) = (1-X)^ALPHA * (1+X)^BETA.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) (1-X)**ALPHA * (1+X)**BETA * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
    Thanks to Xu Xiang of Fudan University for pointing out that
    an earlier implementation of this routine was incorrect!
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    14 May 2007
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input, int ORDER, the order of the quadrature rule to be computed.
  
    Input, double ALPHA, BETA, the exponents of (1-X) and
    (1+X) in the quadrature rule.  For simple Gauss-Legendre quadrature,
    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
  
    Output, double XTAB[ORDER], the abscissas.
  
    Output, double WEIGHT[ORDER], the weights.
*/
{
  double an;
  double bn;
  double cc;
  double delta;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double r3;
  double x;
/*
  Check ALPHA and BETA.
*/
  if ( alpha <= -1.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "JACOBI_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  -1.0 < ALPHA is required.\n" );
    exit ( 1 );
  }

  if ( beta <= -1.0 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "JACOBI_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  -1.0 < BETA is required.\n" );
    exit ( 1 );
  }

  double b[order];
  double c[order];
/*
  Set the recursion coefficients.
*/
  for ( i = 1; i <= order; i++ )
  {
    if ( alpha + beta == 0.0 || beta - alpha == 0.0 )
    {
      b[i-1] = 0.0;
    }
    else
    {
      b[i-1] = ( alpha + beta ) * ( beta - alpha ) /
             ( ( alpha + beta + ( double ) ( 2 * i ) )
             * ( alpha + beta + ( double ) ( 2 * i - 2 ) ) );
    }

    if ( i == 1 )
    {
      c[i-1] = 0.0;
    }
    else
    {
      c[i-1] = 4.0 * ( double ) ( i - 1 )
         * ( alpha + ( double ) ( i - 1 ) )
          * ( beta + ( double ) ( i - 1 ) )
            * ( alpha + beta + ( double ) ( i - 1 ) ) /
            ( ( alpha + beta + ( double ) ( 2 * i - 1 ) )
            * pow ( alpha + beta + ( double ) ( 2 * i - 2 ), 2 )
            * ( alpha + beta + ( double ) ( 2 * i - 3 ) ) );
    }
  }

  delta = exp ( log_gamma ( alpha        + 1.0 )
              + log_gamma (         beta + 1.0 )
              - log_gamma ( alpha + beta + 2.0 ) );
  prod = 1.0;
  for ( i = 2; i <= order; i++ )
  {
    prod = prod * c[i-1];
  }
  cc = delta * pow ( 2.0, alpha + beta + 1.0 ) * prod;

  for ( i = 1; i <= order; i++ )
  {
    if ( i == 1 )
    {
      an = alpha / ( double ) ( order );
      bn = beta / ( double ) ( order );

      r1 = ( 1.0 + alpha )
        * ( 2.78 / ( 4.0 + ( double ) ( order * order ) )
        + 0.768 * an / ( double ) ( order ) );

      r2 = 1.0 + 1.48 * an + 0.96 * bn
        + 0.452 * an * an + 0.83 * an * bn;

      x = ( r2 - r1 ) / r2;
    }
    else if ( i == 2 )
    {
      r1 = ( 4.1 + alpha ) /
        ( ( 1.0 + alpha ) * ( 1.0 + 0.156 * alpha ) );

      r2 = 1.0 + 0.06 * ( ( double ) ( order ) - 8.0 ) *
        ( 1.0 + 0.12 * alpha ) / ( double ) ( order );

      r3 = 1.0 + 0.012 * beta *
        ( 1.0 + 0.25 * r8_abs ( alpha ) ) / ( double ) ( order );

      x = x - r1 * r2 * r3 * ( 1.0 - x );
    }
    else if ( i == 3 )
    {
      r1 = ( 1.67 + 0.28 * alpha ) / ( 1.0 + 0.37 * alpha );

      r2 = 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 )
        / ( double ) ( order );

      r3 = 1.0 + 8.0 * beta /
        ( ( 6.28 + beta ) * ( double ) ( order * order ) );

      x = x - r1 * r2 * r3 * ( xtab[0] - x );
    }
    else if ( i < order - 1 )
    {
      x = 3.0 * xtab[i-2] - 3.0 * xtab[i-3] + xtab[i-4];
    }
    else if ( i == order - 1 )
    {
      r1 = ( 1.0 + 0.235 * beta ) / ( 0.766 + 0.119 * beta );

      r2 = 1.0 / ( 1.0 + 0.639
        * ( ( double ) ( order ) - 4.0 )
        / ( 1.0 + 0.71 * ( ( double ) ( order ) - 4.0 ) ) );

      r3 = 1.0 / ( 1.0 + 20.0 * alpha / ( ( 7.5 + alpha ) *
        ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }
    else if ( i == order )
    {
      r1 = ( 1.0 + 0.37 * beta ) / ( 1.67 + 0.28 * beta );

      r2 = 1.0 /
        ( 1.0 + 0.22 * ( ( double ) ( order ) - 8.0 )
        / ( double ) ( order ) );

      r3 = 1.0 / ( 1.0 + 8.0 * alpha /
        ( ( 6.28 + alpha ) * ( double ) ( order * order ) ) );

      x = x + r1 * r2 * r3 * ( x - xtab[i-3] );
    }
    jacobi_root ( &x, order, alpha, beta, &dp2, &p1, b, c );

    xtab[i-1] = x;
    weight[i-1] = cc / ( dp2 * p1 );
  }
/*
  Reverse the order of the values.
*/
  r8vec_reverse ( order, xtab );
  r8vec_reverse ( order, weight );

  return;
}
/******************************************************************************/

double jacobi_integral ( int expon, double alpha, double beta )

/******************************************************************************/
/*
  Purpose:
  
    JACOBI_INTEGRAL evaluates the integral of a monomial with Jacobi weight.
  
  Discussion:
  
    VALUE = Integral ( -1 <= X <= +1 ) x^EXPON (1-x)^ALPHA (1+x)^BETA dx
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    08 September 2007
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int EXPON, the exponent.
  
    Input, double ALPHA, the exponent of (1-X) in the weight factor.
  
    Input, double BETA, the exponent of (1+X) in the weight factor.
  
    Output, double JACOBI_INTEGRAL, the value of the integral.
*/
{
  double arg1;
  double arg2;
  double arg3;
  double arg4;
  double c;
  double s;
  double value;
  double value1;
  double value2;

  c = ( double ) ( expon );

  if ( ( expon % 2 ) == 0 )
  {
    s = +1.0;
  }
  else
  {
    s = -1.0;
  }

  arg1 = - alpha;
  arg2 =   1.0 + c;
  arg3 =   2.0 + beta + c;
  arg4 = - 1.0;

  value1 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  arg1 = - beta;
  arg2 =   1.0 + c;
  arg3 =   2.0 + alpha + c;
  arg4 = - 1.0;

  value2 = r8_hyper_2f1 ( arg1, arg2, arg3, arg4 );

  value = r8_gamma ( 1.0 + c ) * (
      s * r8_gamma ( 1.0 + beta  ) * value1
    / r8_gamma ( 2.0 + beta  + c )
    +     r8_gamma ( 1.0 + alpha ) * value2
    / r8_gamma ( 2.0 + alpha + c ) );

  return value;
}
/******************************************************************************/

void jacobi_recur ( double *p2, double *dp2, double *p1, double x, int order,
  double alpha, double beta, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:
  
    JACOBI_RECUR finds the value and derivative of a Jacobi polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    04 May 2006
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Output, double *P2, the value of J(ORDER)(X).
  
    Output, double *DP2, the value of J'(ORDER)(X).
  
    Output, double *P1, the value of J(ORDER-1)(X).
  
    Input, double X, the point at which polynomials are evaluated.
  
    Input, int ORDER, the order of the polynomial to be computed.
  
    Input, double ALPHA, BETA, the exponents of (1-X) and
    (1+X) in the quadrature rule.
  
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x + ( alpha - beta ) / ( alpha + beta + 2.0 );
  *dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2 = ( x - b[i-1] ) *  ( *p1 ) - c[i-1] * p0;
    *dp2 = ( x - b[i-1] ) * dp1 + ( *p1 ) - c[i-1] * dp0;
  }
  return;
}
/******************************************************************************/

void jacobi_root ( double *x, int order, double alpha, double beta,
  double *dp2, double *p1, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:
  
    JACOBI_ROOT improves an approximate root of a Jacobi polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    04 May 2006
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input/output, double *X, the approximate root, which
    should be improved on output.
  
    Input, int ORDER, the order of the polynomial to be computed.
  
    Input, double ALPHA, BETA, the exponents of (1-X) and
    (1+X) in the quadrature rule.
  
    Output, double *DP2, the value of J'(ORDER)(X).
  
    Output, double *P1, the value of J(ORDER-1)(X).
  
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    jacobi_recur ( &p2, dp2, p1, *x, order, alpha, beta, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      return;
    }
  }
  return;
}
/******************************************************************************/

void kronrod_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    KRONROD_SET sets abscissas and weights for Gauss-Kronrod quadrature.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x-1] = 1.0.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAb[I) )
  
    A Kronrod rule is used in conjunction with a lower order
    Gauss rule, and provides an efficient error estimation.
  
    The error may be estimated as the difference in the two integral
    approximations.
  
    The efficiency comes about because the Kronrod uses the abscissas
    of the Gauss rule, thus saving on the number of function evaluations
    necessary.  If the Kronrod rule were replaced by a Gauss rule of
    the same order, a higher precision integral estimate would be
    made, but the function would have to be evaluated at many more
    points.
  
    The Gauss Kronrod pair of rules involves an ( ORDER + 1 ) / 2
    point Gauss-Legendre rule and an ORDER point Kronrod rule.
    Thus, the 15 point Kronrod rule should be paired with the
    Gauss-Legendre 7 point rule.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    02 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Robert Piessens, Elise de Doncker-Kapenger,
    Christian Ueberhuber, David Kahaner,
    QUADPACK, A Subroutine Package for Automatic Integration,
    Springer Verlag, 1983.
  
  Parameters:
  
    Input, int ORDER, the order of the rule, which may be
    15, 21, 31 or 41, corresponding to Gauss-Legendre rules of
    order 7, 10, 15 or 20.
  
    Output, double XTAB[ORDER], the abscissas of the rule, which
    are symmetrically places in [-1,1].
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are positive, symmetric, and should sum to 2.
*/
{
  if ( order == 15 )
  {
    xtab[1-1] =  - 0.9914553711208126E+00;
    xtab[2-1] =  - 0.9491079123427585E+00;
    xtab[3-1] =  - 0.8648644233597691E+00;
    xtab[4-1] =  - 0.7415311855993944E+00;
    xtab[5-1] =  - 0.5860872354676911E+00;
    xtab[6-1] =  - 0.4058451513773972E+00;
    xtab[7-1] =  - 0.2077849550789850E+00;
    xtab[8-1] =    0.0E+00;
    xtab[9-1] =    0.2077849550789850E+00;
    xtab[10-1] =   0.4058451513773972E+00;
    xtab[11-1] =   0.5860872354676911E+00;
    xtab[12-1] =   0.7415311855993944E+00;
    xtab[13-1] =   0.8648644233597691E+00;
    xtab[14-1] =   0.9491079123427585E+00;
    xtab[15-1] =   0.9914553711208126E+00;

    weight[1-1] =  0.2293532201052922E-01;
    weight[2-1] =  0.6309209262997855E-01;
    weight[3-1] =  0.1047900103222502E+00;
    weight[4-1] =  0.1406532597155259E+00;
    weight[5-1] =  0.1690047266392679E+00;
    weight[6-1] =  0.1903505780647854E+00;
    weight[7-1] =  0.2044329400752989E+00;
    weight[8-1] =  0.2094821410847278E+00;
    weight[9-1] =  0.2044329400752989E+00;
    weight[10-1] = 0.1903505780647854E+00;
    weight[11-1] = 0.1690047266392679E+00;
    weight[12-1] = 0.1406532597155259E+00;
    weight[13-1] = 0.1047900103222502E+00;
    weight[14-1] = 0.6309209262997855E-01;
    weight[15-1] = 0.2293532201052922E-01;
  }
  else if ( order == 21 )
  {
    xtab[1-1] =  - 0.9956571630258081E+00;
    xtab[2-1] =  - 0.9739065285171717E+00;
    xtab[3-1] =  - 0.9301574913557082E+00;
    xtab[4-1] =  - 0.8650633666889845E+00;
    xtab[5-1] =  - 0.7808177265864169E+00;
    xtab[6-1] =  - 0.6794095682990244E+00;
    xtab[7-1] =  - 0.5627571346686047E+00;
    xtab[8-1] =  - 0.4333953941292472E+00;
    xtab[9-1] =  - 0.2943928627014602E+00;
    xtab[10-1] = - 0.1488743389816312E+00;
    xtab[11-1] =   0.0E+00;
    xtab[12-1] =   0.1488743389816312E+00;
    xtab[13-1] =   0.2943928627014602E+00;
    xtab[14-1] =   0.4333953941292472E+00;
    xtab[15-1] =   0.5627571346686047E+00;
    xtab[16-1] =   0.6794095682990244E+00;
    xtab[17-1] =   0.7808177265864169E+00;
    xtab[18-1] =   0.8650633666889845E+00;
    xtab[19-1] =   0.9301574913557082E+00;
    xtab[20-1] =   0.9739065285171717E+00;
    xtab[21-1] =   0.9956571630258081E+00;

    weight[1-1] =  0.1169463886737187E-01;
    weight[2-1] =  0.3255816230796473E-01;
    weight[3-1] =  0.5475589657435200E-01;
    weight[4-1] =  0.7503967481091995E-01;
    weight[5-1] =  0.9312545458369761E-01;
    weight[6-1] =  0.1093871588022976E+00;
    weight[7-1] =  0.1234919762620659E+00;
    weight[8-1] =  0.1347092173114733E+00;
    weight[9-1] =  0.1427759385770601E+00;
    weight[10-1] = 0.1477391049013385E+00;
    weight[11-1] = 0.1494455540029169E+00;
    weight[12-1] = 0.1477391049013385E+00;
    weight[13-1] = 0.1427759385770601E+00;
    weight[14-1] = 0.1347092173114733E+00;
    weight[15-1] = 0.1234919762620659E+00;
    weight[16-1] = 0.1093871588022976E+00;
    weight[17-1] = 0.9312545458369761E-01;
    weight[18-1] = 0.7503967481091995E-01;
    weight[19-1] = 0.5475589657435200E-01;
    weight[20-1] = 0.3255816230796473E-01;
    weight[21-1] = 0.1169463886737187E-01;
  }
  else if ( order == 31 )
  {
    xtab[1-1] =  - 0.9980022986933971E+00;
    xtab[2-1] =  - 0.9879925180204854E+00;
    xtab[3-1] =  - 0.9677390756791391E+00;
    xtab[4-1] =  - 0.9372733924007059E+00;
    xtab[5-1] =  - 0.8972645323440819E+00;
    xtab[6-1] =  - 0.8482065834104272E+00;
    xtab[7-1] =  - 0.7904185014424659E+00;
    xtab[8-1] =  - 0.7244177313601700E+00;
    xtab[9-1] =  - 0.6509967412974170E+00;
    xtab[10-1] = - 0.5709721726085388E+00;
    xtab[11-1] = - 0.4850818636402397E+00;
    xtab[12-1] = - 0.3941513470775634E+00;
    xtab[13-1] = - 0.2991800071531688E+00;
    xtab[14-1] = - 0.2011940939974345E+00;
    xtab[15-1] = - 0.1011420669187175E+00;
    xtab[16-1] =   0.0E+00;
    xtab[17-1] =   0.1011420669187175E+00;
    xtab[18-1] =   0.2011940939974345E+00;
    xtab[19-1] =   0.2991800071531688E+00;
    xtab[20-1] =   0.3941513470775634E+00;
    xtab[21-1] =   0.4850818636402397E+00;
    xtab[22-1] =   0.5709721726085388E+00;
    xtab[23-1] =   0.6509967412974170E+00;
    xtab[24-1] =   0.7244177313601700E+00;
    xtab[25-1] =   0.7904185014424659E+00;
    xtab[26-1] =   0.8482065834104272E+00;
    xtab[27-1] =   0.8972645323440819E+00;
    xtab[28-1] =   0.9372733924007059E+00;
    xtab[29-1] =   0.9677390756791391E+00;
    xtab[30-1] =   0.9879925180204854E+00;
    xtab[31-1] =   0.9980022986933971E+00;

    weight[1-1] =  0.5377479872923349E-02;
    weight[2-1] =  0.1500794732931612E-01;
    weight[3-1] =  0.2546084732671532E-01;
    weight[4-1] =  0.3534636079137585E-01;
    weight[5-1] =  0.4458975132476488E-01;
    weight[6-1] =  0.5348152469092809E-01;
    weight[7-1] =  0.6200956780067064E-01;
    weight[8-1] =  0.6985412131872826E-01;
    weight[9-1] =  0.7684968075772038E-01;
    weight[10-1] = 0.8308050282313302E-01;
    weight[11-1] = 0.8856444305621177E-01;
    weight[12-1] = 0.9312659817082532E-01;
    weight[13-1] = 0.9664272698362368E-01;
    weight[14-1] = 0.9917359872179196E-01;
    weight[15-1] = 0.1007698455238756E+00;
    weight[16-1] = 0.1013300070147915E+00;
    weight[17-1] = 0.1007698455238756E+00;
    weight[18-1] = 0.9917359872179196E-01;
    weight[19-1] = 0.9664272698362368E-01;
    weight[20-1] = 0.9312659817082532E-01;
    weight[21-1] = 0.8856444305621177E-01;
    weight[22-1] = 0.8308050282313302E-01;
    weight[23-1] = 0.7684968075772038E-01;
    weight[24-1] = 0.6985412131872826E-01;
    weight[25-1] = 0.6200956780067064E-01;
    weight[26-1] = 0.5348152469092809E-01;
    weight[27-1] = 0.4458975132476488E-01;
    weight[28-1] = 0.3534636079137585E-01;
    weight[29-1] = 0.2546084732671532E-01;
    weight[30-1] = 0.1500794732931612E-01;
    weight[31-1] = 0.5377479872923349E-02;
  }
  else if ( order == 41 )
  {
    xtab[1-1] =  - 0.9988590315882777E+00;
    xtab[2-1] =  - 0.9931285991850949E+00;
    xtab[3-1] =  - 0.9815078774502503E+00;
    xtab[4-1] =  - 0.9639719272779138E+00;
    xtab[5-1] =  - 0.9408226338317548E+00;
    xtab[6-1] =  - 0.9122344282513259E+00;
    xtab[7-1] =  - 0.8782768112522820E+00;
    xtab[8-1] =  - 0.8391169718222188E+00;
    xtab[9-1] =  - 0.7950414288375512E+00;
    xtab[10-1] = - 0.7463319064601508E+00;
    xtab[11-1] = - 0.6932376563347514E+00;
    xtab[12-1] = - 0.6360536807265150E+00;
    xtab[13-1] = - 0.5751404468197103E+00;
    xtab[14-1] = - 0.5108670019508271E+00;
    xtab[15-1] = - 0.4435931752387251E+00;
    xtab[16-1] = - 0.3737060887154196E+00;
    xtab[17-1] = - 0.3016278681149130E+00;
    xtab[18-1] = - 0.2277858511416451E+00;
    xtab[19-1] = - 0.1526054652409227E+00;
    xtab[20-1] = - 0.7652652113349733E-01;
    xtab[21-1] =   0.0E+00;
    xtab[22-1] =   0.7652652113349733E-01;
    xtab[23-1] =   0.1526054652409227E+00;
    xtab[24-1] =   0.2277858511416451E+00;
    xtab[25-1] =   0.3016278681149130E+00;
    xtab[26-1] =   0.3737060887154196E+00;
    xtab[27-1] =   0.4435931752387251E+00;
    xtab[28-1] =   0.5108670019508271E+00;
    xtab[29-1] =   0.5751404468197103E+00;
    xtab[30-1] =   0.6360536807265150E+00;
    xtab[31-1] =   0.6932376563347514E+00;
    xtab[32-1] =   0.7463319064601508E+00;
    xtab[33-1] =   0.7950414288375512E+00;
    xtab[34-1] =   0.8391169718222188E+00;
    xtab[35-1] =   0.8782768112522820E+00;
    xtab[36-1] =   0.9122344282513259E+00;
    xtab[37-1] =   0.9408226338317548E+00;
    xtab[38-1] =   0.9639719272779138E+00;
    xtab[39-1] =   0.9815078774502503E+00;
    xtab[40-1] =   0.9931285991850949E+00;
    xtab[41-1] =   0.9988590315882777E+00;

    weight[1-1] =  0.3073583718520532E-02;
    weight[2-1] =  0.8600269855642942E-02;
    weight[3-1] =  0.1462616925697125E-01;
    weight[4-1] =  0.2038837346126652E-01;
    weight[5-1] =  0.2588213360495116E-01;
    weight[6-1] =  0.3128730677703280E-01;
    weight[7-1] =  0.3660016975820080E-01;
    weight[8-1] =  0.4166887332797369E-01;
    weight[9-1] =  0.4643482186749767E-01;
    weight[10-1] = 0.5094457392372869E-01;
    weight[11-1] = 0.5519510534828599E-01;
    weight[12-1] = 0.5911140088063957E-01;
    weight[13-1] = 0.6265323755478117E-01;
    weight[14-1] = 0.6583459713361842E-01;
    weight[15-1] = 0.6864867292852162E-01;
    weight[16-1] = 0.7105442355344407E-01;
    weight[17-1] = 0.7303069033278667E-01;
    weight[18-1] = 0.7458287540049919E-01;
    weight[19-1] = 0.7570449768455667E-01;
    weight[20-1] = 0.7637786767208074E-01;
    weight[21-1] = 0.7660071191799966E-01;
    weight[22-1] = 0.7637786767208074E-01;
    weight[23-1] = 0.7570449768455667E-01;
    weight[24-1] = 0.7458287540049919E-01;
    weight[25-1] = 0.7303069033278667E-01;
    weight[26-1] = 0.7105442355344407E-01;
    weight[27-1] = 0.6864867292852162E-01;
    weight[28-1] = 0.6583459713361842E-01;
    weight[29-1] = 0.6265323755478117E-01;
    weight[30-1] = 0.5911140088063957E-01;
    weight[31-1] = 0.5519510534828599E-01;
    weight[32-1] = 0.5094457392372869E-01;
    weight[33-1] = 0.4643482186749767E-01;
    weight[34-1] = 0.4166887332797369E-01;
    weight[35-1] = 0.3660016975820080E-01;
    weight[36-1] = 0.3128730677703280E-01;
    weight[37-1] = 0.2588213360495116E-01;
    weight[38-1] = 0.2038837346126652E-01;
    weight[39-1] = 0.1462616925697125E-01;
    weight[40-1] = 0.8600269855642942E-02;
    weight[41-1] = 0.3073583718520532E-02;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "KRONROD_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 15, 21, 31 or 41.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void laguerre_compute ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LAGUERRE_COMPUTE computes a Gauss-Laguerre quadrature rule.
  
  Discussion:
  
    The integration interval is [ 0, +oo ).
  
    The weight function is w(x) = exp ( -x );.
  
    If the integral to approximate is:
  
        Integral ( 0 <= X < +oo ) EXP ( - X ) * F(X) dX
  
    then the quadrature rule is:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
    If the integral to approximate is:
  
        Integral ( A <= X < +oo ) F(X) dX
  
    then the quadrature rule is:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    18 February 2008
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input, int ORDER, the order of the quadrature rule to be computed.
    ORDER must be at least 1.
  
    Output, double XTAB[ORDER], the Gauss-Laguerre abscissas.
  
    Output, double WEIGHT[ORDER], the Gauss-Laguerre weights.
*/
{
  double cc;
  double dp2;
  int i;
  double p1;
  double prod;
  double r1;
  double r2;
  double ratio;
  double x;

  double b[order];
  double c[order];
/*
  Set the recursion coefficients.
*/
  for ( i = 0; i < order; i++ )
  {
    b[i] = ( double ) ( 2 * i + 1 );
  }

  for ( i = 0; i < order; i++ )
  {
    c[i] = ( double ) ( i * i );
  }
  prod = 1.0;
  for ( i = 1; i < order; i++ )
  {
    prod = prod * c[i];
  }
  cc = prod;

  for ( i = 0; i < order; i++ )
  {
/*
  Compute an estimate for the root.
*/
    if ( i == 0 )
    {
      x =  3.0 / ( 1.0 + 2.4 * ( double ) ( order ) );
    }
    else if ( i == 1 )
    {
      x = x + 15.0 / ( 1.0 + 2.5 * ( double ) ( order ) );
    }
    else
    {
      r1 = ( 1.0 + 2.55 * ( double ) ( i - 1 ) )
        / ( 1.9 * ( double ) ( i - 1 ) );

      x = x + r1 * ( x - xtab[i-2] );
    }
/*
  Use iteration to find the root.
*/
    laguerre_root ( &x, order, &dp2, &p1, b, c );
/*
  Set the abscissa and weight.
*/
    xtab[i] = x;
    weight[i] = ( cc / dp2 ) / p1;
  }

  return;
}
/******************************************************************************/

double laguerre_integral ( int expon )

/******************************************************************************/
/*
  Purpose:
  
    LAGUERRE_INTEGRAL evaluates a monomial Laguerre integral.
  
  Discussion:
  
    The integral being computed is
  
      integral ( 0 <= x < +oo ) x^n * exp ( -x ) dx
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    19 February 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int EXPON, the exponent.
    0 <= EXPON.
  
    Output, double EXACT, the value of the integral.
*/
{
  double exact;

  exact = r8_factorial ( expon );

  return exact;
}
/******************************************************************************/

void laguerre_recur ( double *p2, double *dp2, double *p1, double x,
  int order, double b[], double c[] )

/******************************************************************************/
/*
  Purpose:
  
    LAGUERRE_RECUR evaluates a Laguerre polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    18 February 2008
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Output, double *P2, the value of L(ORDER)(X).
  
    Output, double *DP2, the value of L'(ORDER)(X).
  
    Output, double *P1, the value of L(ORDER-1)(X).
  
    Input, double X, the point at which polynomials are evaluated.
  
    Input, int ORDER, the order of the polynomial to be computed.
  
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x - 1.0;
  *dp2 = 1.0;

  for ( i = 1; i < order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2  = ( x - b[i] ) * ( *p1 ) - c[i] * p0;
    *dp2 = ( x - b[i] ) * dp1 + ( *p1 ) - c[i] * dp0;
  }

  return;
}
/******************************************************************************/

void laguerre_root ( double *x, int order, double *dp2, double *p1,
  double b[], double c[] )

/******************************************************************************/
/*
  Purpose:
  
    LAGUERRE_ROOT improves a root of a Laguerre polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    18 February 2008
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input/output, double *X, the approximate root, which
    should be improved on output.
  
    Input, int ORDER, the order of the polynomial to be computed.
  
    Output, double *DP2, the value of L'(ORDER)(X).
  
    Output, double *P1, the value of L(ORDER-1)(X).
  
    Input, double B[ORDER], C[ORDER], the recursion coefficients.
*/
{
  double d;
  double eps;
  double p2;
  int step;
  int step_max = 10;

  eps = r8_epsilon ( );

  for ( step = 1; step <= step_max; step++ )
  {
    laguerre_recur ( &p2, dp2, p1, *x, order, b, c );

    d = p2 / ( *dp2 );
    *x = *x - d;

    if ( r8_abs ( d ) <= eps * ( r8_abs ( *x ) + 1.0 ) )
    {
      break;
    }
  }

  return;
}
/******************************************************************************/

void laguerre_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LAGUERRE_SET sets abscissas and weights for Laguerre quadrature.
  
  Discussion:
  
    The integration interval is [ 0, +oo ).
  
    The weight function is w(x-1] = exp ( -x ).
  
    The abscissas are the zeroes of the Laguerre polynomial L(ORDER)(X).
  
  
    If the integral to approximate is:
  
      Integral ( 0 <= X < +oo ) exp ( -X ) * F(X) dX
  
    then the quadrature rule is:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * f ( XTAB(I) )
  
    If the integral to approximate is:
  
      Integral ( 0 <= X < +oo ) F(X) dX
  
    then the quadrature rule is:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * exp ( XTAB(I) ) * f ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    05 October 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Vladimir Krylov,
    Approximate Calculation of Integrals,
    Dover, 2006,
    ISBN: 0486445798,
    LC: QA311.K713.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 20, 31, 32, 63, 64 or 127.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are positive, and should add to 1.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] =    1.0E+00;
    weight[1-1] =  1.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] =    0.585786437626904951198311275790E+00;
    xtab[2-1] =    0.341421356237309504880168872421E+01;

    weight[1-1] =  0.853553390593273762200422181052E+00;
    weight[2-1] =  0.146446609406726237799577818948E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] =    0.415774556783479083311533873128E+00;
    xtab[2-1] =    0.229428036027904171982205036136E+01;
    xtab[3-1] =    0.628994508293747919686641576551E+01;

    weight[1-1] =  0.711093009929173015449590191143E+00;
    weight[2-1] =  0.278517733569240848801444888457E+00;
    weight[3-1] =  0.103892565015861357489649204007E-01;
  }
  else if ( order == 4 )
  {
    xtab[1-1] =    0.322547689619392311800361943361E+00;
    xtab[2-1] =    0.174576110115834657568681671252E+01;
    xtab[3-1] =    0.453662029692112798327928538496E+01;
    xtab[4-1] =    0.939507091230113312923353644342E+01;

    weight[1-1] =  0.603154104341633601635966023818E+00;
    weight[2-1] =  0.357418692437799686641492017458E+00;
    weight[3-1] =  0.388879085150053842724381681562E-01;
    weight[4-1] =  0.539294705561327450103790567621E-03;
  }
  else if ( order == 5 )
  {
    xtab[1-1] =    0.263560319718140910203061943361E+00;
    xtab[2-1] =    0.141340305910651679221840798019E+01;
    xtab[3-1] =    0.359642577104072208122318658878E+01;
    xtab[4-1] =    0.708581000585883755692212418111E+01;
    xtab[5-1] =    0.126408008442757826594332193066E+02;

    weight[1-1] =  0.521755610582808652475860928792E+00;
    weight[2-1] =  0.398666811083175927454133348144E+00;
    weight[3-1] =  0.759424496817075953876533114055E-01;
    weight[4-1] =  0.361175867992204845446126257304E-02;
    weight[5-1] =  0.233699723857762278911490845516E-04;
  }
  else if ( order == 6 )
  {
    xtab[1-1] =    0.222846604179260689464354826787E+00;
    xtab[2-1] =    0.118893210167262303074315092194E+01;
    xtab[3-1] =    0.299273632605931407769132528451E+01;
    xtab[4-1] =    0.577514356910451050183983036943E+01;
    xtab[5-1] =    0.983746741838258991771554702994E+01;
    xtab[6-1] =    0.159828739806017017825457915674E+02;

    weight[1-1] =  0.458964673949963593568284877709E+00;
    weight[2-1] =  0.417000830772120994113377566193E+00;
    weight[3-1] =  0.113373382074044975738706185098E+00;
    weight[4-1] =  0.103991974531490748989133028469E-01;
    weight[5-1] =  0.261017202814932059479242860001E-03;
    weight[6-1] =  0.898547906429621238825292052825E-06;
  }
  else if ( order == 7 )
  {
    xtab[1-1] =    0.193043676560362413838247885004E+00;
    xtab[2-1] =    0.102666489533919195034519944317E+01;
    xtab[3-1] =    0.256787674495074620690778622666E+01;
    xtab[4-1] =    0.490035308452648456810171437810E+01;
    xtab[5-1] =    0.818215344456286079108182755123E+01;
    xtab[6-1] =    0.127341802917978137580126424582E+02;
    xtab[7-1] =    0.193957278622625403117125820576E+02;

    weight[1-1] =  0.409318951701273902130432880018E+00;
    weight[2-1] =  0.421831277861719779929281005417E+00;
    weight[3-1] =  0.147126348657505278395374184637E+00;
    weight[4-1] =  0.206335144687169398657056149642E-01;
    weight[5-1] =  0.107401014328074552213195962843E-02;
    weight[6-1] =  0.158654643485642012687326223234E-04;
    weight[7-1] =  0.317031547899558056227132215385E-07;
  }
  else if ( order == 8 )
  {
    xtab[1-1] =    0.170279632305100999788861856608E+00;
    xtab[2-1] =    0.903701776799379912186020223555E+00;
    xtab[3-1] =    0.225108662986613068930711836697E+01;
    xtab[4-1] =    0.426670017028765879364942182690E+01;
    xtab[5-1] =    0.704590540239346569727932548212E+01;
    xtab[6-1] =    0.107585160101809952240599567880E+02;
    xtab[7-1] =    0.157406786412780045780287611584E+02;
    xtab[8-1] =    0.228631317368892641057005342974E+02;

    weight[1-1] =  0.369188589341637529920582839376E+00;
    weight[2-1] =  0.418786780814342956076978581333E+00;
    weight[3-1] =  0.175794986637171805699659866777E+00;
    weight[4-1] =  0.333434922612156515221325349344E-01;
    weight[5-1] =  0.279453623522567252493892414793E-02;
    weight[6-1] =  0.907650877335821310423850149336E-04;
    weight[7-1] =  0.848574671627253154486801830893E-06;
    weight[8-1] =  0.104800117487151038161508853552E-08;
  }
  else if ( order == 9 )
  {
    xtab[1-1] =    0.152322227731808247428107073127E+00;
    xtab[2-1] =    0.807220022742255847741419210952E+00;
    xtab[3-1] =    0.200513515561934712298303324701E+01;
    xtab[4-1] =    0.378347397333123299167540609364E+01;
    xtab[5-1] =    0.620495677787661260697353521006E+01;
    xtab[6-1] =    0.937298525168757620180971073215E+01;
    xtab[7-1] =    0.134662369110920935710978818397E+02;
    xtab[8-1] =    0.188335977889916966141498992996E+02;
    xtab[9-1] =    0.263740718909273767961410072937E+02;

    weight[1-1] =  0.336126421797962519673467717606E+00;
    weight[2-1] =  0.411213980423984387309146942793E+00;
    weight[3-1] =  0.199287525370885580860575607212E+00;
    weight[4-1] =  0.474605627656515992621163600479E-01;
    weight[5-1] =  0.559962661079458317700419900556E-02;
    weight[6-1] =  0.305249767093210566305412824291E-03;
    weight[7-1] =  0.659212302607535239225572284875E-05;
    weight[8-1] =  0.411076933034954844290241040330E-07;
    weight[9-1] =  0.329087403035070757646681380323E-10;
  }
  else if ( order == 10 )
  {
    xtab[1-1] =    0.137793470540492430830772505653E+00;
    xtab[2-1] =    0.729454549503170498160373121676E+00;
    xtab[3-1] =    0.180834290174031604823292007575E+01;
    xtab[4-1] =    0.340143369785489951448253222141E+01;
    xtab[5-1] =    0.555249614006380363241755848687E+01;
    xtab[6-1] =    0.833015274676449670023876719727E+01;
    xtab[7-1] =    0.118437858379000655649185389191E+02;
    xtab[8-1] =    0.162792578313781020995326539358E+02;
    xtab[9-1] =    0.219965858119807619512770901956E+02;
    xtab[10-1] =   0.299206970122738915599087933408E+02;

    weight[1-1] =  0.308441115765020141547470834678E+00;
    weight[2-1] =  0.401119929155273551515780309913E+00;
    weight[3-1] =  0.218068287611809421588648523475E+00;
    weight[4-1] =  0.620874560986777473929021293135E-01;
    weight[5-1] =  0.950151697518110055383907219417E-02;
    weight[6-1] =  0.753008388587538775455964353676E-03;
    weight[7-1] =  0.282592334959956556742256382685E-04;
    weight[8-1] =  0.424931398496268637258657665975E-06;
    weight[9-1] =  0.183956482397963078092153522436E-08;
    weight[10-1] = 0.991182721960900855837754728324E-12;
  }
  else if ( order == 11 )
  {
    xtab[1-1] =    0.125796442187967522675794577516E+00;
    xtab[2-1] =    0.665418255839227841678127839420E+00;
    xtab[3-1] =    0.164715054587216930958700321365E+01;
    xtab[4-1] =    0.309113814303525495330195934259E+01;
    xtab[5-1] =    0.502928440157983321236999508366E+01;
    xtab[6-1] =    0.750988786380661681941099714450E+01;
    xtab[7-1] =    0.106059509995469677805559216457E+02;
    xtab[8-1] =    0.144316137580641855353200450349E+02;
    xtab[9-1] =    0.191788574032146786478174853989E+02;
    xtab[10-1] =   0.252177093396775611040909447797E+02;
    xtab[11-1] =   0.334971928471755372731917259395E+02;

    weight[1-1] =  0.284933212894200605056051024724E+00;
    weight[2-1] =  0.389720889527849377937553508048E+00;
    weight[3-1] =  0.232781831848991333940223795543E+00;
    weight[4-1] =  0.765644535461966864008541790132E-01;
    weight[5-1] =  0.143932827673506950918639187409E-01;
    weight[6-1] =  0.151888084648487306984777640042E-02;
    weight[7-1] =  0.851312243547192259720424170600E-04;
    weight[8-1] =  0.229240387957450407857683270709E-05;
    weight[9-1] =  0.248635370276779587373391491114E-07;
    weight[10-1] = 0.771262693369132047028152590222E-10;
    weight[11-1] = 0.288377586832362386159777761217E-13;
  }
  else if ( order == 12 )
  {
    xtab[1-1] =    0.115722117358020675267196428240E+00;
    xtab[2-1] =    0.611757484515130665391630053042E+00;
    xtab[3-1] =    0.151261026977641878678173792687E+01;
    xtab[4-1] =    0.283375133774350722862747177657E+01;
    xtab[5-1] =    0.459922763941834848460572922485E+01;
    xtab[6-1] =    0.684452545311517734775433041849E+01;
    xtab[7-1] =    0.962131684245686704391238234923E+01;
    xtab[8-1] =    0.130060549933063477203460524294E+02;
    xtab[9-1] =    0.171168551874622557281840528008E+02;
    xtab[10-1] =   0.221510903793970056699218950837E+02;
    xtab[11-1] =   0.284879672509840003125686072325E+02;
    xtab[12-1] =   0.370991210444669203366389142764E+02;

    weight[1-1] =  0.264731371055443190349738892056E+00;
    weight[2-1] =  0.377759275873137982024490556707E+00;
    weight[3-1] =  0.244082011319877564254870818274E+00;
    weight[4-1] =  0.904492222116809307275054934667E-01;
    weight[5-1] =  0.201023811546340965226612867827E-01;
    weight[6-1] =  0.266397354186531588105415760678E-02;
    weight[7-1] =  0.203231592662999392121432860438E-03;
    weight[8-1] =  0.836505585681979874533632766396E-05;
    weight[9-1] =  0.166849387654091026116989532619E-06;
    weight[10-1] = 0.134239103051500414552392025055E-08;
    weight[11-1] = 0.306160163503502078142407718971E-11;
    weight[12-1] = 0.814807746742624168247311868103E-15;
  }
  else if ( order == 13 )
  {
    xtab[1-1] =    0.107142388472252310648493376977E+00;
    xtab[2-1] =    0.566131899040401853406036347177E+00;
    xtab[3-1] =    0.139856433645101971792750259921E+01;
    xtab[4-1] =    0.261659710840641129808364008472E+01;
    xtab[5-1] =    0.423884592901703327937303389926E+01;
    xtab[6-1] =    0.629225627114007378039376523025E+01;
    xtab[7-1] =    0.881500194118697804733348868036E+01;
    xtab[8-1] =    0.118614035888112425762212021880E+02;
    xtab[9-1] =    0.155107620377037527818478532958E+02;
    xtab[10-1] =   0.198846356638802283332036594634E+02;
    xtab[11-1] =   0.251852638646777580842970297823E+02;
    xtab[12-1] =   0.318003863019472683713663283526E+02;
    xtab[13-1] =   0.407230086692655795658979667001E+02;

    weight[1-1] =  0.247188708429962621346249185964E+00;
    weight[2-1] =  0.365688822900521945306717530893E+00;
    weight[3-1] =  0.252562420057658502356824288815E+00;
    weight[4-1] =  0.103470758024183705114218631672E+00;
    weight[5-1] =  0.264327544155616157781587735702E-01;
    weight[6-1] =  0.422039604025475276555209292644E-02;
    weight[7-1] =  0.411881770472734774892472527082E-03;
    weight[8-1] =  0.235154739815532386882897300772E-04;
    weight[9-1] =  0.731731162024909910401047197761E-06;
    weight[10-1] = 0.110884162570398067979150974759E-07;
    weight[11-1] = 0.677082669220589884064621459082E-10;
    weight[12-1] = 0.115997995990507606094507145382E-12;
    weight[13-1] = 0.224509320389275841599187226865E-16;
  }
  else if ( order == 14 )
  {
    xtab[1-1] =    0.997475070325975745736829452514E-01;
    xtab[2-1] =    0.526857648851902896404583451502E+00;
    xtab[3-1] =    0.130062912125149648170842022116E+01;
    xtab[4-1] =    0.243080107873084463616999751038E+01;
    xtab[5-1] =    0.393210282229321888213134366778E+01;
    xtab[6-1] =    0.582553621830170841933899983898E+01;
    xtab[7-1] =    0.814024014156514503005978046052E+01;
    xtab[8-1] =    0.109164995073660188408130510904E+02;
    xtab[9-1] =    0.142108050111612886831059780825E+02;
    xtab[10-1] =   0.181048922202180984125546272083E+02;
    xtab[11-1] =   0.227233816282696248232280886985E+02;
    xtab[12-1] =   0.282729817232482056954158923218E+02;
    xtab[13-1] =   0.351494436605924265828643121364E+02;
    xtab[14-1] =   0.443660817111174230416312423666E+02;

    weight[1-1] =  0.231815577144864977840774861104E+00;
    weight[2-1] =  0.353784691597543151802331301273E+00;
    weight[3-1] =  0.258734610245428085987320561144E+00;
    weight[4-1] =  0.115482893556923210087304988673E+00;
    weight[5-1] =  0.331920921593373600387499587137E-01;
    weight[6-1] =  0.619286943700661021678785967675E-02;
    weight[7-1] =  0.739890377867385942425890907080E-03;
    weight[8-1] =  0.549071946684169837857331777667E-04;
    weight[9-1] =  0.240958576408537749675775256553E-05;
    weight[10-1] = 0.580154398167649518088619303904E-07;
    weight[11-1] = 0.681931469248497411961562387084E-09;
    weight[12-1] = 0.322120775189484793980885399656E-11;
    weight[13-1] = 0.422135244051658735159797335643E-14;
    weight[14-1] = 0.605237502228918880839870806281E-18;
  }
  else if ( order == 15 )
  {
    xtab[1-1] =    0.933078120172818047629030383672E-01;
    xtab[2-1] =    0.492691740301883908960101791412E+00;
    xtab[3-1] =    0.121559541207094946372992716488E+01;
    xtab[4-1] =    0.226994952620374320247421741375E+01;
    xtab[5-1] =    0.366762272175143727724905959436E+01;
    xtab[6-1] =    0.542533662741355316534358132596E+01;
    xtab[7-1] =    0.756591622661306786049739555812E+01;
    xtab[8-1] =    0.101202285680191127347927394568E+02;
    xtab[9-1] =    0.131302824821757235640991204176E+02;
    xtab[10-1] =   0.166544077083299578225202408430E+02;
    xtab[11-1] =   0.207764788994487667729157175676E+02;
    xtab[12-1] =   0.256238942267287801445868285977E+02;
    xtab[13-1] =   0.314075191697539385152432196202E+02;
    xtab[14-1] =   0.385306833064860094162515167595E+02;
    xtab[15-1] =   0.480260855726857943465734308508E+02;

    weight[1-1] =  0.218234885940086889856413236448E+00;
    weight[2-1] =  0.342210177922883329638948956807E+00;
    weight[3-1] =  0.263027577941680097414812275022E+00;
    weight[4-1] =  0.126425818105930535843030549378E+00;
    weight[5-1] =  0.402068649210009148415854789871E-01;
    weight[6-1] =  0.856387780361183836391575987649E-02;
    weight[7-1] =  0.121243614721425207621920522467E-02;
    weight[8-1] =  0.111674392344251941992578595518E-03;
    weight[9-1] =  0.645992676202290092465319025312E-05;
    weight[10-1] = 0.222631690709627263033182809179E-06;
    weight[11-1] = 0.422743038497936500735127949331E-08;
    weight[12-1] = 0.392189726704108929038460981949E-10;
    weight[13-1] = 0.145651526407312640633273963455E-12;
    weight[14-1] = 0.148302705111330133546164737187E-15;
    weight[15-1] = 0.160059490621113323104997812370E-19;
  }
  else if ( order == 16 )
  {
    xtab[1-1] =    0.876494104789278403601980973401E-01;
    xtab[2-1] =    0.462696328915080831880838260664E+00;
    xtab[3-1] =    0.114105777483122685687794501811E+01;
    xtab[4-1] =    0.212928364509838061632615907066E+01;
    xtab[5-1] =    0.343708663389320664523510701675E+01;
    xtab[6-1] =    0.507801861454976791292305830814E+01;
    xtab[7-1] =    0.707033853504823413039598947080E+01;
    xtab[8-1] =    0.943831433639193878394724672911E+01;
    xtab[9-1] =    0.122142233688661587369391246088E+02;
    xtab[10-1] =   0.154415273687816170767647741622E+02;
    xtab[11-1] =   0.191801568567531348546631409497E+02;
    xtab[12-1] =   0.235159056939919085318231872752E+02;
    xtab[13-1] =   0.285787297428821403675206137099E+02;
    xtab[14-1] =   0.345833987022866258145276871778E+02;
    xtab[15-1] =   0.419404526476883326354722330252E+02;
    xtab[16-1] =   0.517011603395433183643426971197E+02;

    weight[1-1] =  0.206151714957800994334273636741E+00;
    weight[2-1] =  0.331057854950884165992983098710E+00;
    weight[3-1] =  0.265795777644214152599502020650E+00;
    weight[4-1] =  0.136296934296377539975547513526E+00;
    weight[5-1] =  0.473289286941252189780623392781E-01;
    weight[6-1] =  0.112999000803394532312490459701E-01;
    weight[7-1] =  0.184907094352631086429176783252E-02;
    weight[8-1] =  0.204271915308278460126018338421E-03;
    weight[9-1] =  0.148445868739812987713515067551E-04;
    weight[10-1] = 0.682831933087119956439559590327E-06;
    weight[11-1] = 0.188102484107967321388159920418E-07;
    weight[12-1] = 0.286235024297388161963062629156E-09;
    weight[13-1] = 0.212707903322410296739033610978E-11;
    weight[14-1] = 0.629796700251786778717446214552E-14;
    weight[15-1] = 0.505047370003551282040213233303E-17;
    weight[16-1] = 0.416146237037285519042648356116E-21;
  }
  else if ( order == 17 )
  {
    xtab[1-1] =    0.826382147089476690543986151980E-01;
    xtab[2-1] =    0.436150323558710436375959029847E+00;
    xtab[3-1] =    0.107517657751142857732980316755E+01;
    xtab[4-1] =    0.200519353164923224070293371933E+01;
    xtab[5-1] =    0.323425612404744376157380120696E+01;
    xtab[6-1] =    0.477351351370019726480932076262E+01;
    xtab[7-1] =    0.663782920536495266541643929703E+01;
    xtab[8-1] =    0.884668551116980005369470571184E+01;
    xtab[9-1] =    0.114255293193733525869726151469E+02;
    xtab[10-1] =   0.144078230374813180021982874959E+02;
    xtab[11-1] =   0.178382847307011409290658752412E+02;
    xtab[12-1] =   0.217782682577222653261749080522E+02;
    xtab[13-1] =   0.263153178112487997766149598369E+02;
    xtab[14-1] =   0.315817716804567331343908517497E+02;
    xtab[15-1] =   0.377960938374771007286092846663E+02;
    xtab[16-1] =   0.453757165339889661829258363215E+02;
    xtab[17-1] =   0.553897517898396106640900199790E+02;

    weight[1-1] =  0.195332205251770832145927297697E+00;
    weight[2-1] =  0.320375357274540281336625631970E+00;
    weight[3-1] =  0.267329726357171097238809604160E+00;
    weight[4-1] =  0.145129854358758625407426447473E+00;
    weight[5-1] =  0.544369432453384577793805803066E-01;
    weight[6-1] =  0.143572977660618672917767247431E-01;
    weight[7-1] =  0.266282473557277256843236250006E-02;
    weight[8-1] =  0.343679727156299920611775097985E-03;
    weight[9-1] =  0.302755178378287010943703641131E-04;
    weight[10-1] = 0.176851505323167689538081156159E-05;
    weight[11-1] = 0.657627288681043332199222748162E-07;
    weight[12-1] = 0.146973093215954679034375821888E-08;
    weight[13-1] = 0.181691036255544979555476861323E-10;
    weight[14-1] = 0.109540138892868740297645078918E-12;
    weight[15-1] = 0.261737388222337042155132062413E-15;
    weight[16-1] = 0.167293569314615469085022374652E-18;
    weight[17-1] = 0.106562631627404278815253271162E-22;
  }
  else if ( order == 18 )
  {
    xtab[1-1] =    0.781691666697054712986747615334E-01;
    xtab[2-1] =    0.412490085259129291039101536536E+00;
    xtab[3-1] =    0.101652017962353968919093686187E+01;
    xtab[4-1] =    0.189488850996976091426727831954E+01;
    xtab[5-1] =    0.305435311320265975115241130719E+01;
    xtab[6-1] =    0.450420553888989282633795571455E+01;
    xtab[7-1] =    0.625672507394911145274209116326E+01;
    xtab[8-1] =    0.832782515660563002170470261564E+01;
    xtab[9-1] =    0.107379900477576093352179033397E+02;
    xtab[10-1] =   0.135136562075550898190863812108E+02;
    xtab[11-1] =   0.166893062819301059378183984163E+02;
    xtab[12-1] =   0.203107676262677428561313764553E+02;
    xtab[13-1] =   0.244406813592837027656442257980E+02;
    xtab[14-1] =   0.291682086625796161312980677805E+02;
    xtab[15-1] =   0.346279270656601721454012429438E+02;
    xtab[16-1] =   0.410418167728087581392948614284E+02;
    xtab[17-1] =   0.488339227160865227486586093290E+02;
    xtab[18-1] =   0.590905464359012507037157810181E+02;

    weight[1-1] =  0.185588603146918805623337752284E+00;
    weight[2-1] =  0.310181766370225293649597595713E+00;
    weight[3-1] =  0.267866567148536354820854394783E+00;
    weight[4-1] =  0.152979747468074906553843082053E+00;
    weight[5-1] =  0.614349178609616527076780103487E-01;
    weight[6-1] =  0.176872130807729312772600233761E-01;
    weight[7-1] =  0.366017976775991779802657207890E-02;
    weight[8-1] =  0.540622787007735323128416319257E-03;
    weight[9-1] =  0.561696505121423113817929049294E-04;
    weight[10-1] = 0.401530788370115755858883625279E-05;
    weight[11-1] = 0.191466985667567497969210011321E-06;
    weight[12-1] = 0.583609526863159412918086289717E-08;
    weight[13-1] = 0.107171126695539012772851317562E-09;
    weight[14-1] = 0.108909871388883385562011298291E-11;
    weight[15-1] = 0.538666474837830887608094323164E-14;
    weight[16-1] = 0.104986597803570340877859934846E-16;
    weight[17-1] = 0.540539845163105364356554467358E-20;
    weight[18-1] = 0.269165326920102862708377715980E-24;
  }
  else if ( order == 19 )
  {
    xtab[1-1] =    0.741587837572050877131369916024E-01;
    xtab[2-1] =    0.391268613319994607337648350299E+00;
    xtab[3-1] =    0.963957343997958058624879377130E+00;
    xtab[4-1] =    0.179617558206832812557725825252E+01;
    xtab[5-1] =    0.289365138187378399116494713237E+01;
    xtab[6-1] =    0.426421553962776647436040018167E+01;
    xtab[7-1] =    0.591814156164404855815360191408E+01;
    xtab[8-1] =    0.786861891533473373105668358176E+01;
    xtab[9-1] =    0.101324237168152659251627415800E+02;
    xtab[10-1] =   0.127308814638423980045092979656E+02;
    xtab[11-1] =   0.156912783398358885454136069861E+02;
    xtab[12-1] =   0.190489932098235501532136429732E+02;
    xtab[13-1] =   0.228508497608294829323930586693E+02;
    xtab[14-1] =   0.271606693274114488789963947149E+02;
    xtab[15-1] =   0.320691222518622423224362865906E+02;
    xtab[16-1] =   0.377129058012196494770647508283E+02;
    xtab[17-1] =   0.443173627958314961196067736013E+02;
    xtab[18-1] =   0.523129024574043831658644222420E+02;
    xtab[19-1] =   0.628024231535003758413504690673E+02;

    weight[1-1] =  0.176768474915912502251035479815E+00;
    weight[2-1] =  0.300478143607254379482156807712E+00;
    weight[3-1] =  0.267599547038175030772695440648E+00;
    weight[4-1] =  0.159913372135580216785512147895E+00;
    weight[5-1] =  0.682493799761491134552355368344E-01;
    weight[6-1] =  0.212393076065443249244062193091E-01;
    weight[7-1] =  0.484162735114839596725013121019E-02;
    weight[8-1] =  0.804912747381366766594647138204E-03;
    weight[9-1] =  0.965247209315350170843161738801E-04;
    weight[10-1] = 0.820730525805103054408982992869E-05;
    weight[11-1] = 0.483056672473077253944806671560E-06;
    weight[12-1] = 0.190499136112328569993615674552E-07;
    weight[13-1] = 0.481668463092806155766936380273E-09;
    weight[14-1] = 0.734825883955114437684376840171E-11;
    weight[15-1] = 0.620227538757261639893719012423E-13;
    weight[16-1] = 0.254143084301542272371866857954E-15;
    weight[17-1] = 0.407886129682571235007187465134E-18;
    weight[18-1] = 0.170775018759383706100412325084E-21;
    weight[19-1] = 0.671506464990818995998969111749E-26;
  }
  else if ( order == 20 )
  {
    xtab[1-1] =    0.705398896919887533666890045842E-01;
    xtab[2-1] =    0.372126818001611443794241388761E+00;
    xtab[3-1] =    0.916582102483273564667716277074E+00;
    xtab[4-1] =    0.170730653102834388068768966741E+01;
    xtab[5-1] =    0.274919925530943212964503046049E+01;
    xtab[6-1] =    0.404892531385088692237495336913E+01;
    xtab[7-1] =    0.561517497086161651410453988565E+01;
    xtab[8-1] =    0.745901745367106330976886021837E+01;
    xtab[9-1] =    0.959439286958109677247367273428E+01;
    xtab[10-1] =   0.120388025469643163096234092989E+02;
    xtab[11-1] =   0.148142934426307399785126797100E+02;
    xtab[12-1] =   0.179488955205193760173657909926E+02;
    xtab[13-1] =   0.214787882402850109757351703696E+02;
    xtab[14-1] =   0.254517027931869055035186774846E+02;
    xtab[15-1] =   0.299325546317006120067136561352E+02;
    xtab[16-1] =   0.350134342404790000062849359067E+02;
    xtab[17-1] =   0.408330570567285710620295677078E+02;
    xtab[18-1] =   0.476199940473465021399416271529E+02;
    xtab[19-1] =   0.558107957500638988907507734445E+02;
    xtab[20-1] =   0.665244165256157538186403187915E+02;

    weight[1-1] =  0.168746801851113862149223899689E+00;
    weight[2-1] =  0.291254362006068281716795323812E+00;
    weight[3-1] =  0.266686102867001288549520868998E+00;
    weight[4-1] =  0.166002453269506840031469127816E+00;
    weight[5-1] =  0.748260646687923705400624639615E-01;
    weight[6-1] =  0.249644173092832210728227383234E-01;
    weight[7-1] =  0.620255084457223684744754785395E-02;
    weight[8-1] =  0.114496238647690824203955356969E-02;
    weight[9-1] =  0.155741773027811974779809513214E-03;
    weight[10-1] = 0.154014408652249156893806714048E-04;
    weight[11-1] = 0.108648636651798235147970004439E-05;
    weight[12-1] = 0.533012090955671475092780244305E-07;
    weight[13-1] = 0.175798117905058200357787637840E-08;
    weight[14-1] = 0.372550240251232087262924585338E-10;
    weight[15-1] = 0.476752925157819052449488071613E-12;
    weight[16-1] = 0.337284424336243841236506064991E-14;
    weight[17-1] = 0.115501433950039883096396247181E-16;
    weight[18-1] = 0.153952214058234355346383319667E-19;
    weight[19-1] = 0.528644272556915782880273587683E-23;
    weight[20-1] = 0.165645661249902329590781908529E-27;
  }
  else if ( order == 31 )
  {
    xtab[  1-1] =   0.45901947621108290743496080275224E-01;
    xtab[  2-1] =   0.24198016382477204890408974151714E+00;
    xtab[  3-1] =   0.59525389422235073707330165005414E+00;
    xtab[  4-1] =    1.1066894995329987162111308789792E+00;
    xtab[  5-1] =    1.7775956928747727211593727482675E+00;
    xtab[  6-1] =    2.6097034152566806503893375925315E+00;
    xtab[  7-1] =    3.6051968023400442698805817554243E+00;
    xtab[  8-1] =    4.7667470844717611313629127271123E+00;
    xtab[  9-1] =    6.0975545671817409269925429328463E+00;
    xtab[ 10-1] =    7.6014009492331374229360106942867E+00;
    xtab[ 11-1] =    9.2827143134708894182536695297710E+00;
    xtab[ 12-1] =    11.146649755619291358993815629587E+00;
    xtab[ 13-1] =    13.199189576244998522464925028637E+00;
    xtab[ 14-1] =    15.447268315549310075809325891801E+00;
    xtab[ 15-1] =    17.898929826644757646725793817752E+00;
    xtab[ 16-1] =    20.563526336715822170743048968779E+00;
    xtab[ 17-1] =    23.451973482011858591050255575933E+00;
    xtab[ 18-1] =    26.577081352118260459975876986478E+00;
    xtab[ 19-1] =    29.953990872346445506951917840024E+00;
    xtab[ 20-1] =    33.600759532902202735410313885784E+00;
    xtab[ 21-1] =    37.539164407330440882887902558001E+00;
    xtab[ 22-1] =    41.795830870182219981347945853330E+00;
    xtab[ 23-1] =    46.403866806411123136029227604386E+00;
    xtab[ 24-1] =    51.405314476797755161861461088395E+00;
    xtab[ 25-1] =    56.854992868715843620511922055660E+00;
    xtab[ 26-1] =    62.826855908786321453677523304806E+00;
    xtab[ 27-1] =    69.425277191080345623322251656443E+00;
    xtab[ 28-1] =    76.807047763862732837609972285484E+00;
    xtab[ 29-1] =    85.230358607545669169387065607043E+00;
    xtab[ 30-1] =    95.188939891525629981308606853957E+00;
    xtab[ 31-1] =    107.95224382757871475002440117666E+00;

    weight[  1-1] =   0.11252789550372583820847728082801E+00;
    weight[  2-1] =   0.21552760818089123795222505285045E+00;
    weight[  3-1] =   0.23830825164569654731905788089234E+00;
    weight[  4-1] =   0.19538830929790229249915303390711E+00;
    weight[  5-1] =   0.12698283289306190143635272904602E+00;
    weight[  6-1] =   0.67186168923899300670929441993508E-01;
    weight[  7-1] =   0.29303224993879487404888669311974E-01;
    weight[  8-1] =   0.10597569915295736089529380314433E-01;
    weight[  9-1] =   0.31851272582386980320974842433019E-02;
    weight[ 10-1] =   0.79549548307940382922092149012477E-03;
    weight[ 11-1] =   0.16480052126636687317862967116412E-03;
    weight[ 12-1] =   0.28229237864310816393860971468993E-04;
    weight[ 13-1] =   0.39802902551008580387116174900106E-05;
    weight[ 14-1] =   0.45931839841801061673729694510289E-06;
    weight[ 15-1] =   0.43075545187731100930131457465897E-07;
    weight[ 16-1] =   0.32551249938271570855175749257884E-08;
    weight[ 17-1] =   0.19620246675410594996247151593142E-09;
    weight[ 18-1] =   0.93190499086617587129534716431331E-11;
    weight[ 19-1] =   0.34377541819411620520312597898311E-12;
    weight[ 20-1] =   0.96795247130446716997405035776206E-14;
    weight[ 21-1] =   0.20368066110115247398010624219291E-15;
    weight[ 22-1] =   0.31212687280713526831765358632585E-17;
    weight[ 23-1] =   0.33729581704161052453395678308350E-19;
    weight[ 24-1] =   0.24672796386616696011038363242541E-21;
    weight[ 25-1] =   0.11582201904525643634834564576593E-23;
    weight[ 26-1] =   0.32472922591425422434798022809020E-26;
    weight[ 27-1] =   0.49143017308057432740820076259666E-29;
    weight[ 28-1] =   0.34500071104808394132223135953806E-32;
    weight[ 29-1] =   0.87663710117162041472932760732881E-36;
    weight[ 30-1] =   0.50363643921161490411297172316582E-40;
    weight[ 31-1] =   0.19909984582531456482439549080330E-45;
  }
  else if ( order == 32 )
  {
    xtab[ 1-1] =   0.04448936583326720E+00;
    xtab[ 2-1] =   0.2345261095196173E+00;
    xtab[ 3-1] =   0.5768846293018861E+00;
    xtab[ 4-1] =    1.072448753817818E+00;
    xtab[ 5-1] =    1.722408776444646E+00;
    xtab[ 6-1] =    2.528336706425794E+00;
    xtab[ 7-1] =    3.492213273021993E+00;
    xtab[ 8-1] =    4.616456769749767E+00;
    xtab[ 9-1] =    5.903958504174245E+00;
    xtab[10-1] =    7.358126733186242E+00;
    xtab[11-1] =    8.982940924212595E+00;
    xtab[12-1] =    10.78301863253997E+00;
    xtab[13-1] =    12.76369798674272E+00;
    xtab[14-1] =    14.93113975552256E+00;
    xtab[15-1] =    17.29245433671532E+00;
    xtab[16-1] =    19.85586094033605E+00;
    xtab[17-1] =    22.63088901319678E+00;
    xtab[18-1] =    25.62863602245925E+00;
    xtab[19-1] =    28.86210181632347E+00;
    xtab[20-1] =    32.34662915396473E+00;
    xtab[21-1] =    36.10049480575197E+00;
    xtab[22-1] =    40.14571977153944E+00;
    xtab[23-1] =    44.50920799575494E+00;
    xtab[24-1] =    49.22439498730864E+00;
    xtab[25-1] =    54.33372133339691E+00;
    xtab[26-1] =    59.89250916213402E+00;
    xtab[27-1] =    65.97537728793505E+00;
    xtab[28-1] =    72.68762809066271E+00;
    xtab[29-1] =    80.18744697791352E+00;
    xtab[30-1] =    88.73534041789240E+00;
    xtab[31-1] =    98.82954286828397E+00;
    xtab[32-1] =    111.7513980979377E+00;

    weight[ 1-1] =   0.1092183419523677E+00;
    weight[ 2-1] =   0.2104431079388168E+00;
    weight[ 3-1] =   0.2352132296698471E+00;
    weight[ 4-1] =   0.1959033359728783E+00;
    weight[ 5-1] =   0.1299837862860714E+00;
    weight[ 6-1] =   0.7057862386571766E-01;
    weight[ 7-1] =   0.3176091250917535E-01;
    weight[ 8-1] =   0.1191821483483859E-01;
    weight[ 9-1] =   0.3738816294611523E-02;
    weight[10-1] =   0.9808033066149506E-03;
    weight[11-1] =   0.2148649188013644E-03;
    weight[12-1] =   0.3920341967987944E-04;
    weight[13-1] =   0.5934541612868650E-05;
    weight[14-1] =   0.7416404578667550E-06;
    weight[15-1] =   0.7604567879120788E-07;
    weight[16-1] =   0.6350602226625806E-08;
    weight[17-1] =   0.4281382971040924E-09;
    weight[18-1] =   0.2305899491891338E-10;
    weight[19-1] =   0.9799379288727102E-12;
    weight[20-1] =   0.3237801657729275E-13;
    weight[21-1] =   0.8171823443420754E-15;
    weight[22-1] =   0.1542133833393825E-16;
    weight[23-1] =   0.2119792290163629E-18;
    weight[24-1] =   0.2054429673788038E-20;
    weight[25-1] =   0.1346982586637393E-22;
    weight[26-1] =   0.5661294130397363E-25;
    weight[27-1] =   0.1418560545463052E-27;
    weight[28-1] =   0.1913375494454211E-30;
    weight[29-1] =   0.1192248760098218E-33;
    weight[30-1] =   0.2671511219240120E-37;
    weight[31-1] =   0.1338616942106269E-41;
    weight[32-1] =   0.4510536193898970E-47;
  }
  else if ( order == 63 )
  {
    xtab[  1-1] =   0.22768893732576153785994330248562E-01;
    xtab[  2-1] =   0.11998325242727824715771416426383E+00;
    xtab[  3-1] =   0.29494185444770149577427738517405E+00;
    xtab[  4-1] =   0.54779087896237725363865073775856E+00;
    xtab[  5-1] =   0.87869061179931901673895567052285E+00;
    xtab[  6-1] =    1.2878464335919706302309207788611E+00;
    xtab[  7-1] =    1.7755123815388553763979463268728E+00;
    xtab[  8-1] =    2.3419925567085989256055628337716E+00;
    xtab[  9-1] =    2.9876423223246473939976731053629E+00;
    xtab[ 10-1] =    3.7128695992018000346299637413422E+00;
    xtab[ 11-1] =    4.5181363349503584391105568561550E+00;
    xtab[ 12-1] =    5.4039601781825946286902599782736E+00;
    xtab[ 13-1] =    6.3709163787865330220392250891777E+00;
    xtab[ 14-1] =    7.4196399339311711154888493199004E+00;
    xtab[ 15-1] =    8.5508280008403328312589048722235E+00;
    xtab[ 16-1] =    9.7652425999245366807004592977996E+00;
    xtab[ 17-1] =    11.063713635140661736220550410604E+00;
    xtab[ 18-1] =    12.447142262356492749798687569289E+00;
    xtab[ 19-1] =    13.916504641057818562912967008183E+00;
    xtab[ 20-1] =    15.472856110036296424777143607779E+00;
    xtab[ 21-1] =    17.117335833863588753116900303886E+00;
    xtab[ 22-1] =    18.851171974154856850873483787506E+00;
    xtab[ 23-1] =    20.675687448056515660377265667433E+00;
    xtab[ 24-1] =    22.592306346311528381292277759986E+00;
    xtab[ 25-1] =    24.602561094972638883700642760037E+00;
    xtab[ 26-1] =    26.708100458737343969779087998829E+00;
    xtab[ 27-1] =    28.910698500451382640177718103234E+00;
    xtab[ 28-1] =    31.212264631175912885477773820802E+00;
    xtab[ 29-1] =    33.614854909101154836598842888345E+00;
    xtab[ 30-1] =    36.120684774484823056306328740825E+00;
    xtab[ 31-1] =    38.732143442933582145626041607663E+00;
    xtab[ 32-1] =    41.451810222318741191114726181363E+00;
    xtab[ 33-1] =    44.282473071479233839358857134636E+00;
    xtab[ 34-1] =    47.227149784295686898935095231536E+00;
    xtab[ 35-1] =    50.289112264240695761749021839419E+00;
    xtab[ 36-1] =    53.471914456788652808348280619542E+00;
    xtab[ 37-1] =    56.779424636342062213099781057119E+00;
    xtab[ 38-1] =    60.215862909019862886417550114424E+00;
    xtab[ 39-1] =    63.785845004235974631701139601836E+00;
    xtab[ 40-1] =    67.494433702293885830374325695045E+00;
    xtab[ 41-1] =    71.347199604295266286654803376075E+00;
    xtab[ 42-1] =    75.350293425653234254290504744279E+00;
    xtab[ 43-1] =    79.510532629986309149555391354778E+00;
    xtab[ 44-1] =    83.835506080872257843339817658508E+00;
    xtab[ 45-1] =    88.333701570354369086112766326498E+00;
    xtab[ 46-1] =    93.014662728558547405303399037100E+00;
    xtab[ 47-1] =    97.889184147578140043386727677112E+00;
    xtab[ 48-1] =    102.96955690741381650783952746778E+00;
    xtab[ 49-1] =    108.26988161961595392226350967206E+00;
    xtab[ 50-1] =    113.80647350287462738934485955901E+00;
    xtab[ 51-1] =    119.59839538830458666962452963285E+00;
    xtab[ 52-1] =    125.66817255856119431291196303280E+00;
    xtab[ 53-1] =    132.04277272091165746585590583045E+00;
    xtab[ 54-1] =    138.75498418103789078167590567526E+00;
    xtab[ 55-1] =    145.84541318313540358283994248439E+00;
    xtab[ 56-1] =    153.36548459497863623710815962660E+00;
    xtab[ 57-1] =    161.38215194813761243562172669592E+00;
    xtab[ 58-1] =    169.98570600665839438795175301156E+00;
    xtab[ 59-1] =    179.30366247401580910251827858515E+00;
    xtab[ 60-1] =    189.52789596532475473668721332981E+00;
    xtab[ 61-1] =    200.97521159924656741628671841018E+00;
    xtab[ 62-1] =    214.25368536638788642698056296400E+00;
    xtab[ 63-1] =    230.93465747089703971246562985079E+00;

    weight[  1-1] =   0.57118633213868979811587283390476E-01;
    weight[  2-1] =   0.12067476090640395283319932036351E+00;
    weight[  3-1] =   0.15925001096581873723870561096472E+00;
    weight[  4-1] =   0.16875178327560799234596192963585E+00;
    weight[  5-1] =   0.15366641977668956696193711310131E+00;
    weight[  6-1] =   0.12368770614716481641086652261948E+00;
    weight[  7-1] =   0.89275098854848671545279150057422E-01;
    weight[  8-1] =   0.58258485446105944957571825725160E-01;
    weight[  9-1] =   0.34546657545992580874717085812508E-01;
    weight[ 10-1] =   0.18675685985714656798286552591203E-01;
    weight[ 11-1] =   0.92233449044093536528490075241649E-02;
    weight[ 12-1] =   0.41671250684839592762582663470209E-02;
    weight[ 13-1] =   0.17238120299900582715386728541955E-02;
    weight[ 14-1] =   0.65320845029716311169340559359043E-03;
    weight[ 15-1] =   0.22677644670909586952405173207471E-03;
    weight[ 16-1] =   0.72127674154810668410750270234861E-04;
    weight[ 17-1] =   0.21011261180466484598811536851241E-04;
    weight[ 18-1] =   0.56035500893357212749181536071292E-05;
    weight[ 19-1] =   0.13673642785604888017836641282292E-05;
    weight[ 20-1] =   0.30507263930195817240736097189550E-06;
    weight[ 21-1] =   0.62180061839309763559981775409241E-07;
    weight[ 22-1] =   0.11566529551931711260022448996296E-07;
    weight[ 23-1] =   0.19614588267565478081534781863335E-08;
    weight[ 24-1] =   0.30286171195709411244334756404054E-09;
    weight[ 25-1] =   0.42521344539400686769012963452599E-10;
    weight[ 26-1] =   0.54202220578073819334698791381873E-11;
    weight[ 27-1] =   0.62627306838597672554166850420603E-12;
    weight[ 28-1] =   0.65474443156573322992307089591924E-13;
    weight[ 29-1] =   0.61815575808729181846302500000047E-14;
    weight[ 30-1] =   0.52592721363507381404263991342633E-15;
    weight[ 31-1] =   0.40230920092646484015391506025408E-16;
    weight[ 32-1] =   0.27600740511819536505013824207729E-17;
    weight[ 33-1] =   0.16936946756968296053322009855265E-18;
    weight[ 34-1] =   0.92689146872177087314963772462726E-20;
    weight[ 35-1] =   0.45093739060365632939780140603959E-21;
    weight[ 36-1] =   0.19435162876132376573629962695374E-22;
    weight[ 37-1] =   0.73926270895169207037999639194513E-24;
    weight[ 38-1] =   0.24714364154434632615980126000066E-25;
    weight[ 39-1] =   0.72288649446741597655145390616476E-27;
    weight[ 40-1] =   0.18407617292614039362985209905608E-28;
    weight[ 41-1] =   0.40583498566841960105759537058880E-30;
    weight[ 42-1] =   0.77000496416438368114463925286343E-32;
    weight[ 43-1] =   0.12488505764999334328843314866038E-33;
    weight[ 44-1] =   0.17185000226767010697663950619912E-35;
    weight[ 45-1] =   0.19896372636672396938013975755522E-37;
    weight[ 46-1] =   0.19199671378804058267713164416870E-39;
    weight[ 47-1] =   0.15278588285522166920459714708240E-41;
    weight[ 48-1] =   0.99054752688842142955854138884590E-44;
    weight[ 49-1] =   0.51597523673029211884228858692990E-46;
    weight[ 50-1] =   0.21249846664084111245693912887783E-48;
    weight[ 51-1] =   0.67903852766852910591172042494884E-51;
    weight[ 52-1] =   0.16466654148296177467908300517887E-53;
    weight[ 53-1] =   0.29509065402691055027053659375033E-56;
    weight[ 54-1] =   0.37838420647571051984882241014675E-59;
    weight[ 55-1] =   0.33358130068542431878174667995217E-62;
    weight[ 56-1] =   0.19223461022273880981363303073329E-65;
    weight[ 57-1] =   0.67812696961083016872779388922288E-69;
    weight[ 58-1] =   0.13404752802440604607620468935693E-72;
    weight[ 59-1] =   0.13109745101805029757648048223928E-76;
    weight[ 60-1] =   0.52624863881401787388694579143866E-81;
    weight[ 61-1] =   0.63780013856587414257760666006511E-86;
    weight[ 62-1] =   0.12997078942372924566347473916943E-91;
    weight[ 63-1] =   0.10008511496968754063443740168421E-98;
  }
  else if ( order == 64 )
  {
    xtab[ 1-1] =   0.02241587414670647E+00;
    xtab[ 2-1] =   0.1181225120967683E+00;
    xtab[ 3-1] =   0.2903657440180379E+00;
    xtab[ 4-1] =   0.5392862212279776E+00;
    xtab[ 5-1] =   0.8650370046481133E+00;
    xtab[ 6-1] =    1.267814040775243E+00;
    xtab[ 7-1] =    1.747859626059437E+00;
    xtab[ 8-1] =    2.305463739307508E+00;
    xtab[ 9-1] =    2.940965156725251E+00;
    xtab[10-1] =    3.654752650207290E+00;
    xtab[11-1] =    4.447266343313095E+00;
    xtab[12-1] =    5.318999254496392E+00;
    xtab[13-1] =    6.270499046923656E+00;
    xtab[14-1] =    7.302370002587396E+00;
    xtab[15-1] =    8.415275239483025E+00;
    xtab[16-1] =    9.609939192796110E+00;
    xtab[17-1] =    10.88715038388637E+00;
    xtab[18-1] =    12.24776450424430E+00;
    xtab[19-1] =    13.69270784554751E+00;
    xtab[20-1] =    15.22298111152473E+00;
    xtab[21-1] =    16.83966365264874E+00;
    xtab[22-1] =    18.54391817085919E+00;
    xtab[23-1] =    20.33699594873024E+00;
    xtab[24-1] =    22.22024266595088E+00;
    xtab[25-1] =    24.19510487593326E+00;
    xtab[26-1] =    26.26313722711849E+00;
    xtab[27-1] =    28.42601052750102E+00;
    xtab[28-1] =    30.68552076752597E+00;
    xtab[29-1] =    33.04359923643783E+00;
    xtab[30-1] =    35.50232389114121E+00;
    xtab[31-1] =    38.06393216564647E+00;
    xtab[32-1] =    40.73083544445863E+00;
    xtab[33-1] =    43.50563546642153E+00;
    xtab[34-1] =    46.39114297861619E+00;
    xtab[35-1] =    49.39039902562469E+00;
    xtab[36-1] =    52.50669934134631E+00;
    xtab[37-1] =    55.74362241327838E+00;
    xtab[38-1] =    59.10506191901710E+00;
    xtab[39-1] =    62.59526440015139E+00;
    xtab[40-1] =    66.21887325124756E+00;
    xtab[41-1] =    69.98098037714684E+00;
    xtab[42-1] =    73.88718723248296E+00;
    xtab[43-1] =    77.94367743446313E+00;
    xtab[44-1] =    82.15730377831930E+00;
    xtab[45-1] =    86.53569334945652E+00;
    xtab[46-1] =    91.08737561313309E+00;
    xtab[47-1] =    95.82194001552074E+00;
    xtab[48-1] =    100.7502319695140E+00;
    xtab[49-1] =    105.8845994687999E+00;
    xtab[50-1] =    111.2392075244396E+00;
    xtab[51-1] =    116.8304450513065E+00;
    xtab[52-1] =    122.6774602685386E+00;
    xtab[53-1] =    128.8028787692377E+00;
    xtab[54-1] =    135.2337879495258E+00;
    xtab[55-1] =    142.0031214899315E+00;
    xtab[56-1] =    149.1516659000494E+00;
    xtab[57-1] =    156.7310751326712E+00;
    xtab[58-1] =    164.8086026551505E+00;
    xtab[59-1] =    173.4749468364243E+00;
    xtab[60-1] =    182.8582046914315E+00;
    xtab[61-1] =    193.1511360370729E+00;
    xtab[62-1] =    204.6720284850595E+00;
    xtab[63-1] =    218.0318519353285E+00;
    xtab[64-1] =    234.8095791713262E+00;

    weight[ 1-1] =   0.05625284233891167E+00;
    weight[ 2-1] =   0.1190239873124324E+00;
    weight[ 3-1] =   0.1574964038621455E+00;
    weight[ 4-1] =   0.1675470504157663E+00;
    weight[ 5-1] =   0.1533528557792350E+00;
    weight[ 6-1] =   0.1242210536093002E+00;
    weight[ 7-1] =   0.9034230098648614E-01;
    weight[ 8-1] =   0.5947775576835436E-01;
    weight[ 9-1] =   0.3562751890403622E-01;
    weight[10-1] =   0.1948041043116633E-01;
    weight[11-1] =   0.9743594899381976E-02;
    weight[12-1] =   0.4464310364166273E-02;
    weight[13-1] =   0.1875359581323072E-02;
    weight[14-1] =   0.7226469815750038E-03;
    weight[15-1] =   0.2554875328334950E-03;
    weight[16-1] =   0.8287143534396861E-04;
    weight[17-1] =   0.2465686396788580E-04;
    weight[18-1] =   0.6726713878829671E-05;
    weight[19-1] =   0.1681785369964099E-05;
    weight[20-1] =   0.3850812981546702E-06;
    weight[21-1] =   0.8068728040990521E-07;
    weight[22-1] =   0.1545723706757701E-07;
    weight[23-1] =   0.2704480147617480E-08;
    weight[24-1] =   0.4316775475427202E-09;
    weight[25-1] =   0.6277752541761450E-10;
    weight[26-1] =   0.8306317376288941E-11;
    weight[27-1] =   0.9984031787220194E-12;
    weight[28-1] =   0.1088353887116663E-12;
    weight[29-1] =   0.1074017403441591E-13;
    weight[30-1] =   0.9575737231574444E-15;
    weight[31-1] =   0.7697028023648578E-16;
    weight[32-1] =   0.5564881137454025E-17;
    weight[33-1] =   0.3609756409010444E-18;
    weight[34-1] =   0.2095095369548946E-19;
    weight[35-1] =   0.1084793301097549E-20;
    weight[36-1] =   0.4994699486363812E-22;
    weight[37-1] =   0.2037836974598821E-23;
    weight[38-1] =   0.7339537564278870E-25;
    weight[39-1] =   0.2323783082198704E-26;
    weight[40-1] =   0.6438234706908816E-28;
    weight[41-1] =   0.1553121095788267E-29;
    weight[42-1] =   0.3244250092019555E-31;
    weight[43-1] =   0.5832386267836164E-33;
    weight[44-1] =   0.8963254833102871E-35;
    weight[45-1] =   0.1168703989550731E-36;
    weight[46-1] =   0.1282055984359980E-38;
    weight[47-1] =   0.1172094937404993E-40;
    weight[48-1] =   0.8835339672328639E-43;
    weight[49-1] =   0.5424955590306203E-45;
    weight[50-1] =   0.2675542666678911E-47;
    weight[51-1] =   0.1042917031411369E-49;
    weight[52-1] =   0.3152902351957753E-52;
    weight[53-1] =   0.7229541910647446E-55;
    weight[54-1] =   0.1224235301230072E-57;
    weight[55-1] =   0.1482168504901918E-60;
    weight[56-1] =   0.1232519348814539E-63;
    weight[57-1] =   0.6691499004571088E-67;
    weight[58-1] =   0.2220465941850503E-70;
    weight[59-1] =   0.4120946094738935E-74;
    weight[60-1] =   0.3774399061896465E-78;
    weight[61-1] =   0.1414115052917608E-82;
    weight[62-1] =   0.1591833064041367E-87;
    weight[63-1] =   0.2989484348860626E-93;
    weight[64-1] =   0.02089063508436960E-99;
  }
  else if ( order == 127 )
  {
    xtab[  1-1] =   0.11339635298518611691893169631306E-01;
    xtab[  2-1] =   0.59749753435726620281348237057387E-01;
    xtab[  3-1] =   0.14685098690746167612388223687431E+00;
    xtab[  4-1] =   0.27267590735859553131378008278900E+00;
    xtab[  5-1] =   0.43724600644192665554577035869932E+00;
    xtab[  6-1] =   0.64058688222566929533576416399983E+00;
    xtab[  7-1] =   0.88272968639058364481487653650042E+00;
    xtab[  8-1] =    1.1637114160166537661560584700951E+00;
    xtab[  9-1] =    1.4835750152834613891313584861012E+00;
    xtab[ 10-1] =    1.8423694351613565380686320809853E+00;
    xtab[ 11-1] =    2.2401496839579024244513315656522E+00;
    xtab[ 12-1] =    2.6769768780141303692167869961238E+00;
    xtab[ 13-1] =    3.1529182957082825565771508308846E+00;
    xtab[ 14-1] =    3.6680474360304752540226339926515E+00;
    xtab[ 15-1] =    4.2224440823301888455977876667425E+00;
    xtab[ 16-1] =    4.8161943715870502475665535087286E+00;
    xtab[ 17-1] =    5.4493908694559416755862178908416E+00;
    xtab[ 18-1] =    6.1221326512997254193944584763155E+00;
    xtab[ 19-1] =    6.8345253894122668112237994973336E+00;
    xtab[ 20-1] =    7.5866814466367472174205986836847E+00;
    xtab[ 21-1] =    8.3787199765932725254842120659452E+00;
    xtab[ 22-1] =    9.2107670307426558777922506102445E+00;
    xtab[ 23-1] =    10.082955672528643809166439353647E+00;
    xtab[ 24-1] =    10.995426098858125429803147358780E+00;
    xtab[ 25-1] =    11.948325769197725997610605127857E+00;
    xtab[ 26-1] =    12.941809542585531053723381098192E+00;
    xtab[ 27-1] =    13.976039822878506520014405668679E+00;
    xtab[ 28-1] =    15.051186712579523631574796365435E+00;
    xtab[ 29-1] =    16.167428175612852922977395051768E+00;
    xtab[ 30-1] =    17.324950209443673446561163712616E+00;
    xtab[ 31-1] =    18.523947026965688560811711309349E+00;
    xtab[ 32-1] =    19.764621248611504104071669386884E+00;
    xtab[ 33-1] =    21.047184105173183606877044020054E+00;
    xtab[ 34-1] =    22.371855651855542817648123918101E+00;
    xtab[ 35-1] =    23.738864994122497183652313788712E+00;
    xtab[ 36-1] =    25.148450525937368234077278385644E+00;
    xtab[ 37-1] =    26.600860181041749607253384279755E+00;
    xtab[ 38-1] =    28.096351697964619201753961292129E+00;
    xtab[ 39-1] =    29.635192899504178910610227138642E+00;
    xtab[ 40-1] =    31.217661987479759144214467152615E+00;
    xtab[ 41-1] =    32.844047853610430460522951341338E+00;
    xtab[ 42-1] =    34.514650407441149149105635947422E+00;
    xtab[ 43-1] =    36.229780922306804019615388508885E+00;
    xtab[ 44-1] =    37.989762400399956435968780140278E+00;
    xtab[ 45-1] =    39.794929958089961778396437141707E+00;
    xtab[ 46-1] =    41.645631232730180705153990897484E+00;
    xtab[ 47-1] =    43.542226812286859549950892993822E+00;
    xtab[ 48-1] =    45.485090689228791137996151336673E+00;
    xtab[ 49-1] =    47.474610740231964719468766599146E+00;
    xtab[ 50-1] =    49.511189233379087716728884584381E+00;
    xtab[ 51-1] =    51.595243364671244443182771266934E+00;
    xtab[ 52-1] =    53.727205825819316758288140069145E+00;
    xtab[ 53-1] =    55.907525405447553305830605991732E+00;
    xtab[ 54-1] =    58.136667626022439197077526025660E+00;
    xtab[ 55-1] =    60.415115419018590295707192053805E+00;
    xtab[ 56-1] =    62.743369841051809700207126742685E+00;
    xtab[ 57-1] =    65.121950833949996311956025417139E+00;
    xtab[ 58-1] =    67.551398031997886314411872443149E+00;
    xtab[ 59-1] =    70.032271619884584511229871192030E+00;
    xtab[ 60-1] =    72.565153245206849090888669416801E+00;
    xtab[ 61-1] =    75.150646989739935299354362325096E+00;
    xtab[ 62-1] =    77.789380404085816000647405462136E+00;
    xtab[ 63-1] =    80.482005610750729205803962926758E+00;
    xtab[ 64-1] =    83.229200481195914886796120019048E+00;
    xtab[ 65-1] =    86.031669892953582966798238732643E+00;
    xtab[ 66-1] =    88.890147073512051099652518544282E+00;
    xtab[ 67-1] =    91.805395038358177994971250170499E+00;
    xtab[ 68-1] =    94.778208131331583205387031034825E+00;
    xtab[ 69-1] =    97.809413676305116411054110115424E+00;
    xtab[ 70-1] =    100.89987375017285940371939762172E+00;
    xtab[ 71-1] =    104.05048708821598934704076845022E+00;
    xtab[ 72-1] =    107.26219113414600428423116401414E+00;
    xtab[ 73-1] =    110.53596424851500530602771351277E+00;
    xtab[ 74-1] =    113.87282809075839485348376187652E+00;
    xtab[ 75-1] =    117.27385019192517774095477886379E+00;
    xtab[ 76-1] =    120.74014673718880106173978002719E+00;
    xtab[ 77-1] =    124.27288557955698354259506446928E+00;
    xtab[ 78-1] =    127.87328950885942645093841745425E+00;
    xtab[ 79-1] =    131.54263980314366921809377742137E+00;
    xtab[ 80-1] =    135.28228009311836970132738106369E+00;
    xtab[ 81-1] =    139.09362057432970013964422086977E+00;
    xtab[ 82-1] =    142.97814260643601776808227753574E+00;
    xtab[ 83-1] =    146.93740374437366549441080969072E+00;
    xtab[ 84-1] =    150.97304325252187127492511437460E+00;
    xtab[ 85-1] =    155.08678816034612572229641420609E+00;
    xtab[ 86-1] =    159.28045992663288235401956989889E+00;
    xtab[ 87-1] =    163.55598178957571104015967182053E+00;
    xtab[ 88-1] =    167.91538689194360134245547184721E+00;
    xtab[ 89-1] =    172.36082728473812536838156191681E+00;
    xtab[ 90-1] =    176.89458392960192176311674993508E+00;
    xtab[ 91-1] =    181.51907784036813069227528834025E+00;
    xtab[ 92-1] =    186.23688252828112373861202530357E+00;
    xtab[ 93-1] =    191.05073794450929196790836610789E+00;
    xtab[ 94-1] =    195.96356614879879837839002542988E+00;
    xtab[ 95-1] =    200.97848897600025153696475526130E+00;
    xtab[ 96-1] =    206.09884802468871112127283042753E+00;
    xtab[ 97-1] =    211.32822735671655260572377256981E+00;
    xtab[ 98-1] =    216.67047937658230323477089465777E+00;
    xtab[ 99-1] =    222.12975445929687246267304963754E+00;
    xtab[100-1] =    227.71053502072232419089132431317E+00;
    xtab[101-1] =    233.41767488282602453367775322563E+00;
    xtab[102-1] =    239.25644498830308620018749667089E+00;
    xtab[103-1] =    245.23258677871567172531254018984E+00;
    xtab[104-1] =    251.35237488718128030005500991754E+00;
    xtab[105-1] =    257.62269123792061413076191882313E+00;
    xtab[106-1] =    264.05111322908240551754377241831E+00;
    xtab[107-1] =    270.64601945722796749299111718606E+00;
    xtab[108-1] =    277.41671750163651071798388218104E+00;
    xtab[109-1] =    284.37359974220870326674402873120E+00;
    xtab[110-1] =    291.52833521346495719581282021650E+00;
    xtab[111-1] =    298.89410837028248600878895615414E+00;
    xtab[112-1] =    306.48591978262611320418112423947E+00;
    xtab[113-1] =    314.32096986471177487400007507615E+00;
    xtab[114-1] =    322.41915589128679683349440361344E+00;
    xtab[115-1] =    330.80372663802405651933847334878E+00;
    xtab[116-1] =    339.50216127832433747735367595958E+00;
    xtab[117-1] =    348.54737559472697355480761787441E+00;
    xtab[118-1] =    357.97942028029845454049007443090E+00;
    xtab[119-1] =    367.84794520076004578858341422871E+00;
    xtab[120-1] =    378.21590623135532818332979188889E+00;
    xtab[121-1] =    389.16539141251004101579475325153E+00;
    xtab[122-1] =    400.80729331451702589996361286427E+00;
    xtab[123-1] =    413.29853681779384418008260081859E+00;
    xtab[124-1] =    426.87579153663675538288509017051E+00;
    xtab[125-1] =    441.93085485310841412460309271842E+00;
    xtab[126-1] =    459.21804639888429981971267313224E+00;
    xtab[127-1] =    480.69378263388373859704269229304E+00;

    weight[  1-1] =   0.28773246692000124355770010301506E-01;
    weight[  2-1] =   0.63817468175134649363480949265236E-01;
    weight[  3-1] =   0.91919669721570571389864194652717E-01;
    weight[  4-1] =   0.11054167914413766381245463002967E+00;
    weight[  5-1] =   0.11879771633375850188328329422643E+00;
    weight[  6-1] =   0.11737818530052695148804451630074E+00;
    weight[  7-1] =   0.10819305984180551488335145581193E+00;
    weight[  8-1] =   0.93827075290489628080377261401107E-01;
    weight[  9-1] =   0.76966450960588843995822485928431E-01;
    weight[ 10-1] =   0.59934903912939714332570730063476E-01;
    weight[ 11-1] =   0.44417742073889001371708316272923E-01;
    weight[ 12-1] =   0.31385080966252320983009372215062E-01;
    weight[ 13-1] =   0.21172316041924506411370709025015E-01;
    weight[ 14-1] =   0.13650145364230541652171185564626E-01;
    weight[ 15-1] =   0.84172852710599172279366657385445E-02;
    weight[ 16-1] =   0.49674990059882760515912858620175E-02;
    weight[ 17-1] =   0.28069903895001884631961957446400E-02;
    weight[ 18-1] =   0.15192951003941952460445341057817E-02;
    weight[ 19-1] =   0.78789028751796084086217287140548E-03;
    weight[ 20-1] =   0.39156751064868450584507324648999E-03;
    weight[ 21-1] =   0.18652434268825860550093566260060E-03;
    weight[ 22-1] =   0.85173160415576621908809828160247E-04;
    weight[ 23-1] =   0.37285639197853037712145321577724E-04;
    weight[ 24-1] =   0.15648416791712993947447805296768E-04;
    weight[ 25-1] =   0.62964340695224829035692735524979E-05;
    weight[ 26-1] =   0.24288929711328724574541379938222E-05;
    weight[ 27-1] =   0.89824607890051007201922871545035E-06;
    weight[ 28-1] =   0.31844174740760353710742966328091E-06;
    weight[ 29-1] =   0.10821272905566839211861807542741E-06;
    weight[ 30-1] =   0.35245076750635536015902779085340E-07;
    weight[ 31-1] =   0.11001224365719347407063839761738E-07;
    weight[ 32-1] =   0.32904079616717932125329343003261E-08;
    weight[ 33-1] =   0.94289145237889976419772700772988E-09;
    weight[ 34-1] =   0.25882578904668318184050195309296E-09;
    weight[ 35-1] =   0.68047437103370762630942259017560E-10;
    weight[ 36-1] =   0.17131398805120837835399564475632E-10;
    weight[ 37-1] =   0.41291744524052865469443922304935E-11;
    weight[ 38-1] =   0.95264189718807273220707664873469E-12;
    weight[ 39-1] =   0.21032604432442425932962942047474E-12;
    weight[ 40-1] =   0.44427151938729352860940434285789E-13;
    weight[ 41-1] =   0.89760500362833703323319846405449E-14;
    weight[ 42-1] =   0.17341511407769287074627948346848E-14;
    weight[ 43-1] =   0.32028099548988356631494379835210E-15;
    weight[ 44-1] =   0.56531388950793682022660742095189E-16;
    weight[ 45-1] =   0.95329672799026591234588044025896E-17;
    weight[ 46-1] =   0.15353453477310142565288509437552E-17;
    weight[ 47-1] =   0.23608962179467365686057842132176E-18;
    weight[ 48-1] =   0.34648742794456611332193876653230E-19;
    weight[ 49-1] =   0.48515241897086461320126957663545E-20;
    weight[ 50-1] =   0.64786228633519813428137373790678E-21;
    weight[ 51-1] =   0.82476020965403242936448553126316E-22;
    weight[ 52-1] =   0.10005361880214719793491658282977E-22;
    weight[ 53-1] =   0.11561395116207304954233181263632E-23;
    weight[ 54-1] =   0.12719342731167922655612134264961E-24;
    weight[ 55-1] =   0.13316584714165372967340004160814E-25;
    weight[ 56-1] =   0.13261218454678944033646108509198E-26;
    weight[ 57-1] =   0.12554995447643949807286074138324E-27;
    weight[ 58-1] =   0.11294412178579462703240913107219E-28;
    weight[ 59-1] =   0.96491020279562119228500608131696E-30;
    weight[ 60-1] =   0.78241846768302099396733076955632E-31;
    weight[ 61-1] =   0.60181503542219626658249939076636E-32;
    weight[ 62-1] =   0.43882482704961741551510518054138E-33;
    weight[ 63-1] =   0.30314137647517256304035802501863E-34;
    weight[ 64-1] =   0.19826016543944539545224676057020E-35;
    weight[ 65-1] =   0.12267623373665926559013654872402E-36;
    weight[ 66-1] =   0.71763931692508888943812834967620E-38;
    weight[ 67-1] =   0.39659378833836963584113716149270E-39;
    weight[ 68-1] =   0.20688970553868040099581951696677E-40;
    weight[ 69-1] =   0.10179587017979517245268418427523E-41;
    weight[ 70-1] =   0.47200827745986374625714293679649E-43;
    weight[ 71-1] =   0.20606828985553374825744353490744E-44;
    weight[ 72-1] =   0.84627575907305987245899032156188E-46;
    weight[ 73-1] =   0.32661123687088798658026998931647E-47;
    weight[ 74-1] =   0.11833939207883162380564134612682E-48;
    weight[ 75-1] =   0.40211209123895013807243250164050E-50;
    weight[ 76-1] =   0.12799824394111125389430292847476E-51;
    weight[ 77-1] =   0.38123877747548846504399051365162E-53;
    weight[ 78-1] =   0.10612057542701156767898551949650E-54;
    weight[ 79-1] =   0.27571446947200403594113572720812E-56;
    weight[ 80-1] =   0.66772544240928492881306904862856E-58;
    weight[ 81-1] =   0.15052438383868234954068178600268E-59;
    weight[ 82-1] =   0.31538986800113758526689068500772E-61;
    weight[ 83-1] =   0.61326614299483180785237418887960E-63;
    weight[ 84-1] =   0.11048510030324810567549119229368E-64;
    weight[ 85-1] =   0.18410563538091348076979665543900E-66;
    weight[ 86-1] =   0.28323926570052832195543883237652E-68;
    weight[ 87-1] =   0.40154409843763655508670978777418E-70;
    weight[ 88-1] =   0.52351530215683708779772201956106E-72;
    weight[ 89-1] =   0.62634476665005100555787696642851E-74;
    weight[ 90-1] =   0.68612210535666530365348093803922E-76;
    weight[ 91-1] =   0.68651298840956019297134099761855E-78;
    weight[ 92-1] =   0.62581388433728084867318704240915E-80;
    weight[ 93-1] =   0.51833271237514904046803469968027E-82;
    weight[ 94-1] =   0.38893621571918443533108973497673E-84;
    weight[ 95-1] =   0.26357711379476932781525533730623E-86;
    weight[ 96-1] =   0.16078851293917979699005509638883E-88;
    weight[ 97-1] =   0.87978042070968939637972577886624E-91;
    weight[ 98-1] =   0.43013405077495109903408697802188E-93;
    weight[ 99-1] =   0.18713435881342838527144321803729E-95;
    weight[100-1] =   0.72125744708060471675805761366523E-98;
    weight[101-1] =   0.24508746062177874383231742333023E-100;
    weight[102-1] =   0.73042094619470875777647865078327E-103;
    weight[103-1] =   0.18983290818383463537886818579820E-105;
    weight[104-1] =   0.42757400244246684123093264825902E-108;
    weight[105-1] =   0.82894681420515755691423485228897E-111;
    weight[106-1] =   0.13729432219324400013067050156048E-113;
    weight[107-1] =   0.19265464126404973222043166489406E-116;
    weight[108-1] =   0.22693344503301354826140809941334E-119;
    weight[109-1] =   0.22209290603717355061909071271535E-122;
    weight[110-1] =   0.17851087685544512662856555121755E-125;
    weight[111-1] =   0.11630931990387164467431190485525E-128;
    weight[112-1] =   0.60524443584652392290952805077893E-132;
    weight[113-1] =   0.24729569115063528647628375096400E-135;
    weight[114-1] =   0.77789065006489410364997205809045E-139;
    weight[115-1] =   0.18409738662712607039570678274636E-142;
    weight[116-1] =   0.31900921131079114970179071968597E-146;
    weight[117-1] =   0.39179487139174199737617666077555E-150;
    weight[118-1] =   0.32782158394188697053774429820559E-154;
    weight[119-1] =   0.17793590713138888062819640128739E-158;
    weight[120-1] =   0.58882353408932623157467835381214E-163;
    weight[121-1] =   0.10957236509071169877747203273886E-167;
    weight[122-1] =   0.10281621114867000898285076975760E-172;
    weight[123-1] =   0.41704725557697758145816510853967E-178;
    weight[124-1] =   0.58002877720316101774638319601971E-184;
    weight[125-1] =   0.18873507745825517106171619101120E-190;
    weight[126-1] =   0.69106601826730911682786705950895E-198;
    weight[127-1] =   0.43506813201105855628383313334402E-207;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LAGUERRE_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 to 20, 31, 32, 63, 64 or 127.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

double laguerre_sum ( double func ( double x ), double a, int order,
  double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LAGUERRE_SUM carries out Laguerre quadrature over [ A, +oo ).
  
  Discussion:
  
    The simplest Laguerre integral to approximate is the
    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
    it is easy to modify the rule to approximate the integral from
    A to +oo as well.
  
    Another common Laguerre integral to approximate is the
    integral from 0 to +oo of EXP(-X) * X**ALPHA * F(X).
    This routine may be used to sum up the terms of the Laguerre
    rule for such an integral as well.  However, if ALPHA is nonzero,
    then there is no simple way to extend the rule to approximate the
    integral from A to +oo.  The simplest procedures would be
    to approximate the integral from 0 to A.
  
    The integration interval is [ A, +oo ) or [ 0, +oo ).
  
    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x**alpha.
  
    The integral to approximate:
  
      Integral ( A <= X <= +oo ) EXP ( -X ) * F(X) dX
    or
      Integral ( 0 <= X <= +oo ) EXP ( -X ) * X**ALPHA * F(X) dX
  
    The quadrature rule:
  
      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) + A )
    or
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    02 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, double FUNC ( double X ), the name of the function which
    evaluates the integrand.
  
    Input, double A, the beginning of the integration interval.
  
    Input, int ORDER, the order of the rule.
  
    Input, double XTAB[ORDER], the abscissas of the rule.
  
    Input, double WEIGHT[ORDER], the weights of the rule.
  
    Output, double LAGUERRE_SUM, the approximate value of the integral.
*/
{
  int i;
  double result;

  if ( order < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LAGUERRE_SUM - Fatal error!\n" );
    fprintf ( stderr, "  Nonpositive ORDER = %d\n", order );
    exit ( 1 );
  }

  result = 0.0;
  for ( i = 0; i < order; i++ )
  {
    result = result + weight[i] * func ( xtab[i] + a );
  }
  result = exp ( - a ) * result;

  return result;
}
/******************************************************************************/

void legendre_compute ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_COMPUTE computes a Gauss-Legendre quadrature rule.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x) = 1.0.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 August 2007
  
  Author:
  
    Original FORTRAN77 version by Philip Davis, Philip Rabinowitz.
    C version by John Burkardt.
  
  Reference:
  
    Philip Davis, Philip Rabinowitz,
    Methods of Numerical Integration,
    Second Edition,
    Dover, 2007,
    ISBN: 0486453391,
    LC: QA299.3.D28.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be greater than 0.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are positive, symmetric, and should sum to 2.
*/
{
  double d1;
  double d2pn;
  double d3pn;
  double d4pn;
  double dp;
  double dpn;
  double e1;
  double fx;
  double h;
  int i;
  int iback;
  int k;
  int m;
  int mp1mi;
  int ncopy;
  int nmove;
  double p;
  double pi = 3.141592653589793;
  double pk;
  double pkm1;
  double pkp1;
  double t;
  double u;
  double v;
  double x0;
  double xtemp;

  if ( order < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  e1 = ( double ) ( order * ( order + 1 ) );

  m = ( order + 1 ) / 2;

  for ( i = 1; i <= m; i++ )
  {
    mp1mi = m + 1 - i;

    t = ( double ) ( 4 * i - 1 ) * pi / ( double ) ( 4 * order + 2 );

    x0 = cos ( t ) * ( 1.0 - ( 1.0 - 1.0 / ( double ) ( order ) )
      / ( double ) ( 8 * order * order ) );

    pkm1 = 1.0;
    pk = x0;

    for ( k = 2; k <= order; k++ )
    {
      pkp1 = 2.0 * x0 * pk - pkm1 - ( x0 * pk - pkm1 ) / ( double ) ( k );
      pkm1 = pk;
      pk = pkp1;
    }

    d1 = ( double ) ( order ) * ( pkm1 - x0 * pk );

    dpn = d1 / ( 1.0 - x0 * x0 );

    d2pn = ( 2.0 * x0 * dpn - e1 * pk ) / ( 1.0 - x0 * x0 );

    d3pn = ( 4.0 * x0 * d2pn + ( 2.0 - e1 ) * dpn ) / ( 1.0 - x0 * x0 );

    d4pn = ( 6.0 * x0 * d3pn + ( 6.0 - e1 ) * d2pn ) / ( 1.0 - x0 * x0 );

    u = pk / dpn;
    v = d2pn / dpn;
/*
  Initial approximation H:
*/
    h = -u * ( 1.0 + 0.5 * u * ( v + u * ( v * v - d3pn / ( 3.0 * dpn ) ) ) );
/*
  Refine H using one step of Newton's method:
*/
    p = pk + h * ( dpn + 0.5 * h * ( d2pn + h / 3.0
      * ( d3pn + 0.25 * h * d4pn ) ) );

    dp = dpn + h * ( d2pn + 0.5 * h * ( d3pn + h * d4pn / 3.0 ) );

    h = h - p / dp;

    xtemp = x0 + h;

    xtab[mp1mi-1] = xtemp;

    fx = d1 - h * e1 * ( pk + 0.5 * h * ( dpn + h / 3.0
      * ( d2pn + 0.25 * h * ( d3pn + 0.2 * h * d4pn ) ) ) );

    weight[mp1mi-1] = 2.0 * ( 1.0 - xtemp * xtemp ) / ( fx * fx );
  }

  if ( ( order % 2 ) == 1 )
  {
    xtab[0] = 0.0;
  }
/*
  Shift the data up.
*/
  nmove = ( order + 1 ) / 2;
  ncopy = order - nmove;

  for ( i = 1; i <= nmove; i++ )
  {
    iback = order + 1 - i;
    xtab[iback-1] = xtab[iback-ncopy-1];
    weight[iback-1] = weight[iback-ncopy-1];
  }
/*
  Reflect values for the negative abscissas.
*/
  for ( i = 1; i <= order - nmove; i++ )
  {
    xtab[i-1] = - xtab[order-i];
    weight[i-1] = weight[order-i];
  }

  return;
}
/******************************************************************************/

double legendre_integral ( int expon )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_INTEGRAL evaluates a monomial Legendre integral.
  
  Discussion:
  
    To test a Legendre quadrature rule, we use it to approximate the
    integral of a monomial:
  
      integral ( -1 <= x <= +1 ) x^n dx
  
    This routine is given the value of the exponent, and returns the
    exact value of the integral.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    19 February 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int EXPON, the exponent.
  
    Output, double LEGENDRE_INTEGRAL, the value of the exact integral.
*/
{
  double exact;
/*
  Get the exact value of the integral.
*/
  if ( ( expon % 2 ) == 0 )
  {
    exact = 2.0 / ( double ) ( expon + 1 );
  }
  else
  {
    exact = 0.0;
  }

  return exact;
}
/******************************************************************************/

void legendre_recur ( double *p2, double *dp2, double *p1, double x, int order )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_RECUR finds the value and derivative of a Legendre polynomial.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    03 May 2006
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Output, double *P2, the value of P(ORDER)(X).
  
    Output, double *DP2, the value of P'(ORDER)(X).
  
    Output, double *P1, the value of P(ORDER-1)(X).
  
    Input, double X, the point at which polynomials are evaluated.
  
    Input, int ORDER, the order of the polynomial to be computed.
*/
{
  double dp0;
  double dp1;
  int i;
  double p0;

  *p1 = 1.0;
  dp1 = 0.0;

  *p2 = x;
  *dp2 = 1.0;

  for ( i = 2; i <= order; i++ )
  {
    p0 = *p1;
    dp0 = dp1;

    *p1 = *p2;
    dp1 = *dp2;

    *p2 = ( ( double ) ( 2 * i - 1 ) * x * ( *p1 )
          + ( double ) (   - i + 1 )     * p0 )
          / ( double ) (     i     );

    *dp2 = ( ( double ) ( 2 * i - 1 ) * ( ( *p1 ) + x * dp1 )
           - ( double ) (     i - 1 ) * dp0 )
           / ( double ) (     i     );
  }

  return;
}
/******************************************************************************/

void legendre_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET sets abscissas and weights for Gauss-Legendre quadrature.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function w(x-1] = 1.0;
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    Quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
    The quadrature rule will integrate exactly all polynomials up to
    X**(2*ORDER-1).
  
    The abscissas of the rule are the zeroes of the Legendre polynomial
    P(ORDER)(X).
  
    The integral produced by a Gauss-Legendre rule is equal to the
    integral of the unique polynomial of degree ORDER-1 which
    agrees with the function at the ORDER abscissas of the rule.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    02 October 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Vladimir Krylov,
    Approximate Calculation of Integrals,
    Dover, 2006,
    ISBN: 0486445798.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 33, 63, 64, 65 or 127.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are positive, symmetric and should sum to 2.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] =   0.0;

    weight[1-1] = 2.0;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = - 0.577350269189625764509148780502;
    xtab[2-1] =   0.577350269189625764509148780502;

    weight[1-1] = 1.0;
    weight[2-1] = 1.0;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = - 0.774596669241483377035853079956;
    xtab[2-1] =   0.0;
    xtab[3-1] =   0.774596669241483377035853079956;

    weight[1-1] = 5.0 / 9.0;
    weight[2-1] = 8.0 / 9.0;
    weight[3-1] = 5.0 / 9.0;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = - 0.861136311594052575223946488893;
    xtab[2-1] = - 0.339981043584856264802665759103;
    xtab[3-1] =   0.339981043584856264802665759103;
    xtab[4-1] =   0.861136311594052575223946488893;

    weight[1-1] = 0.347854845137453857373063949222;
    weight[2-1] = 0.652145154862546142626936050778;
    weight[3-1] = 0.652145154862546142626936050778;
    weight[4-1] = 0.347854845137453857373063949222;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = - 0.906179845938663992797626878299;
    xtab[2-1] = - 0.538469310105683091036314420700;
    xtab[3-1] =   0.0;
    xtab[4-1] =   0.538469310105683091036314420700;
    xtab[5-1] =   0.906179845938663992797626878299;

    weight[1-1] = 0.236926885056189087514264040720;
    weight[2-1] = 0.478628670499366468041291514836;
    weight[3-1] = 0.568888888888888888888888888889;
    weight[4-1] = 0.478628670499366468041291514836;
    weight[5-1] = 0.236926885056189087514264040720;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = - 0.932469514203152027812301554494;
    xtab[2-1] = - 0.661209386466264513661399595020;
    xtab[3-1] = - 0.238619186083196908630501721681;
    xtab[4-1] =   0.238619186083196908630501721681;
    xtab[5-1] =   0.661209386466264513661399595020;
    xtab[6-1] =   0.932469514203152027812301554494;

    weight[1-1] = 0.171324492379170345040296142173;
    weight[2-1] = 0.360761573048138607569833513838;
    weight[3-1] = 0.467913934572691047389870343990;
    weight[4-1] = 0.467913934572691047389870343990;
    weight[5-1] = 0.360761573048138607569833513838;
    weight[6-1] = 0.171324492379170345040296142173;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = - 0.949107912342758524526189684048;
    xtab[2-1] = - 0.741531185599394439863864773281;
    xtab[3-1] = - 0.405845151377397166906606412077;
    xtab[4-1] =   0.0;
    xtab[5-1] =   0.405845151377397166906606412077;
    xtab[6-1] =   0.741531185599394439863864773281;
    xtab[7-1] =   0.949107912342758524526189684048;

    weight[1-1] = 0.129484966168869693270611432679;
    weight[2-1] = 0.279705391489276667901467771424;
    weight[3-1] = 0.381830050505118944950369775489;
    weight[4-1] = 0.417959183673469387755102040816;
    weight[5-1] = 0.381830050505118944950369775489;
    weight[6-1] = 0.279705391489276667901467771424;
    weight[7-1] = 0.129484966168869693270611432679;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = - 0.960289856497536231683560868569;
    xtab[2-1] = - 0.796666477413626739591553936476;
    xtab[3-1] = - 0.525532409916328985817739049189;
    xtab[4-1] = - 0.183434642495649804939476142360;
    xtab[5-1] =   0.183434642495649804939476142360;
    xtab[6-1] =   0.525532409916328985817739049189;
    xtab[7-1] =   0.796666477413626739591553936476;
    xtab[8-1] =   0.960289856497536231683560868569;

    weight[1-1] = 0.101228536290376259152531354310;
    weight[2-1] = 0.222381034453374470544355994426;
    weight[3-1] = 0.313706645877887287337962201987;
    weight[4-1] = 0.362683783378361982965150449277;
    weight[5-1] = 0.362683783378361982965150449277;
    weight[6-1] = 0.313706645877887287337962201987;
    weight[7-1] = 0.222381034453374470544355994426;
    weight[8-1] = 0.101228536290376259152531354310;
  }
  else if ( order == 9 )
  {
    xtab[1-1] = - 0.968160239507626089835576202904;
    xtab[2-1] = - 0.836031107326635794299429788070;
    xtab[3-1] = - 0.613371432700590397308702039341;
    xtab[4-1] = - 0.324253423403808929038538014643;
    xtab[5-1] =   0.0;
    xtab[6-1] =   0.324253423403808929038538014643;
    xtab[7-1] =   0.613371432700590397308702039341;
    xtab[8-1] =   0.836031107326635794299429788070;
    xtab[9-1] =   0.968160239507626089835576202904;

    weight[1-1] = 0.812743883615744119718921581105E-01;
    weight[2-1] = 0.180648160694857404058472031243;
    weight[3-1] = 0.260610696402935462318742869419;
    weight[4-1] = 0.312347077040002840068630406584;
    weight[5-1] = 0.330239355001259763164525069287;
    weight[6-1] = 0.312347077040002840068630406584;
    weight[7-1] = 0.260610696402935462318742869419;
    weight[8-1] = 0.180648160694857404058472031243;
    weight[9-1] = 0.812743883615744119718921581105E-01;
  }
  else if ( order == 10 )
  {
    xtab[1-1] =  - 0.973906528517171720077964012084;
    xtab[2-1] =  - 0.865063366688984510732096688423;
    xtab[3-1] =  - 0.679409568299024406234327365115;
    xtab[4-1] =  - 0.433395394129247190799265943166;
    xtab[5-1] =  - 0.148874338981631210884826001130;
    xtab[6-1] =    0.148874338981631210884826001130;
    xtab[7-1] =    0.433395394129247190799265943166;
    xtab[8-1] =    0.679409568299024406234327365115;
    xtab[9-1] =    0.865063366688984510732096688423;
    xtab[10-1] =   0.973906528517171720077964012084;

    weight[1-1] =  0.666713443086881375935688098933E-01;
    weight[2-1] =  0.149451349150580593145776339658;
    weight[3-1] =  0.219086362515982043995534934228;
    weight[4-1] =  0.269266719309996355091226921569;
    weight[5-1] =  0.295524224714752870173892994651;
    weight[6-1] =  0.295524224714752870173892994651;
    weight[7-1] =  0.269266719309996355091226921569;
    weight[8-1] =  0.219086362515982043995534934228;
    weight[9-1] =  0.149451349150580593145776339658;
    weight[10-1] = 0.666713443086881375935688098933E-01;
  }
  else if ( order == 11 )
  {
    xtab[1-1] =  - 0.978228658146056992803938001123;
    xtab[2-1] =  - 0.887062599768095299075157769304;
    xtab[3-1] =  - 0.730152005574049324093416252031;
    xtab[4-1] =  - 0.519096129206811815925725669459;
    xtab[5-1] =  - 0.269543155952344972331531985401;
    xtab[6-1] =    0.0;
    xtab[7-1] =    0.269543155952344972331531985401;
    xtab[8-1] =    0.519096129206811815925725669459;
    xtab[9-1] =    0.730152005574049324093416252031;
    xtab[10-1] =   0.887062599768095299075157769304;
    xtab[11-1] =   0.978228658146056992803938001123;

    weight[1-1] =  0.556685671161736664827537204425E-01;
    weight[2-1] =  0.125580369464904624634694299224;
    weight[3-1] =  0.186290210927734251426097641432;
    weight[4-1] =  0.233193764591990479918523704843;
    weight[5-1] =  0.262804544510246662180688869891;
    weight[6-1] =  0.272925086777900630714483528336;
    weight[7-1] =  0.262804544510246662180688869891;
    weight[8-1] =  0.233193764591990479918523704843;
    weight[9-1] =  0.186290210927734251426097641432;
    weight[10-1] = 0.125580369464904624634694299224;
    weight[11-1] = 0.556685671161736664827537204425E-01;
  }
  else if ( order == 12 )
  {
    xtab[1-1] =  - 0.981560634246719250690549090149;
    xtab[2-1] =  - 0.904117256370474856678465866119;
    xtab[3-1] =  - 0.769902674194304687036893833213;
    xtab[4-1] =  - 0.587317954286617447296702418941;
    xtab[5-1] =  - 0.367831498998180193752691536644;
    xtab[6-1] =  - 0.125233408511468915472441369464;
    xtab[7-1] =    0.125233408511468915472441369464;
    xtab[8-1] =    0.367831498998180193752691536644;
    xtab[9-1] =    0.587317954286617447296702418941;
    xtab[10-1] =   0.769902674194304687036893833213;
    xtab[11-1] =   0.904117256370474856678465866119;
    xtab[12-1] =   0.981560634246719250690549090149;

    weight[1-1] =  0.471753363865118271946159614850E-01;
    weight[2-1] =  0.106939325995318430960254718194;
    weight[3-1] =  0.160078328543346226334652529543;
    weight[4-1] =  0.203167426723065921749064455810;
    weight[5-1] =  0.233492536538354808760849898925;
    weight[6-1] =  0.249147045813402785000562436043;
    weight[7-1] =  0.249147045813402785000562436043;
    weight[8-1] =  0.233492536538354808760849898925;
    weight[9-1] =  0.203167426723065921749064455810;
    weight[10-1] = 0.160078328543346226334652529543;
    weight[11-1] = 0.106939325995318430960254718194;
    weight[12-1] = 0.471753363865118271946159614850E-01;
  }
  else if ( order == 13 )
  {
    xtab[1-1] =  - 0.984183054718588149472829448807;
    xtab[2-1] =  - 0.917598399222977965206547836501;
    xtab[3-1] =  - 0.801578090733309912794206489583;
    xtab[4-1] =  - 0.642349339440340220643984606996;
    xtab[5-1] =  - 0.448492751036446852877912852128;
    xtab[6-1] =  - 0.230458315955134794065528121098;
    xtab[7-1] =    0.0;
    xtab[8-1] =    0.230458315955134794065528121098;
    xtab[9-1] =    0.448492751036446852877912852128;
    xtab[10-1] =   0.642349339440340220643984606996;
    xtab[11-1] =   0.801578090733309912794206489583;
    xtab[12-1] =   0.917598399222977965206547836501;
    xtab[13-1] =   0.984183054718588149472829448807;

    weight[1-1] =  0.404840047653158795200215922010E-01;
    weight[2-1] =  0.921214998377284479144217759538E-01;
    weight[3-1] =  0.138873510219787238463601776869;
    weight[4-1] =  0.178145980761945738280046691996;
    weight[5-1] =  0.207816047536888502312523219306;
    weight[6-1] =  0.226283180262897238412090186040;
    weight[7-1] =  0.232551553230873910194589515269;
    weight[8-1] =  0.226283180262897238412090186040;
    weight[9-1] =  0.207816047536888502312523219306;
    weight[10-1] = 0.178145980761945738280046691996;
    weight[11-1] = 0.138873510219787238463601776869;
    weight[12-1] = 0.921214998377284479144217759538E-01;
    weight[13-1] = 0.404840047653158795200215922010E-01;
  }
  else if ( order == 14 )
  {
    xtab[1-1] =  - 0.986283808696812338841597266704;
    xtab[2-1] =  - 0.928434883663573517336391139378;
    xtab[3-1] =  - 0.827201315069764993189794742650;
    xtab[4-1] =  - 0.687292904811685470148019803019;
    xtab[5-1] =  - 0.515248636358154091965290718551;
    xtab[6-1] =  - 0.319112368927889760435671824168;
    xtab[7-1] =  - 0.108054948707343662066244650220;
    xtab[8-1] =    0.108054948707343662066244650220;
    xtab[9-1] =    0.319112368927889760435671824168;
    xtab[10-1] =   0.515248636358154091965290718551;
    xtab[11-1] =   0.687292904811685470148019803019;
    xtab[12-1] =   0.827201315069764993189794742650;
    xtab[13-1] =   0.928434883663573517336391139378;
    xtab[14-1] =   0.986283808696812338841597266704;

    weight[1-1] =  0.351194603317518630318328761382E-01;
    weight[2-1] =  0.801580871597602098056332770629E-01;
    weight[3-1] =  0.121518570687903184689414809072;
    weight[4-1] =  0.157203167158193534569601938624;
    weight[5-1] =  0.185538397477937813741716590125;
    weight[6-1] =  0.205198463721295603965924065661;
    weight[7-1] =  0.215263853463157790195876443316;
    weight[8-1] =  0.215263853463157790195876443316;
    weight[9-1] =  0.205198463721295603965924065661;
    weight[10-1] = 0.185538397477937813741716590125;
    weight[11-1] = 0.157203167158193534569601938624;
    weight[12-1] = 0.121518570687903184689414809072;
    weight[13-1] = 0.801580871597602098056332770629E-01;
    weight[14-1] = 0.351194603317518630318328761382E-01;
  }
  else if ( order == 15 )
  {
    xtab[1-1] =  - 0.987992518020485428489565718587;
    xtab[2-1] =  - 0.937273392400705904307758947710;
    xtab[3-1] =  - 0.848206583410427216200648320774;
    xtab[4-1] =  - 0.724417731360170047416186054614;
    xtab[5-1] =  - 0.570972172608538847537226737254;
    xtab[6-1] =  - 0.394151347077563369897207370981;
    xtab[7-1] =  - 0.201194093997434522300628303395;
    xtab[8-1] =    0.0;
    xtab[9-1] =    0.201194093997434522300628303395;
    xtab[10-1] =   0.394151347077563369897207370981;
    xtab[11-1] =   0.570972172608538847537226737254;
    xtab[12-1] =   0.724417731360170047416186054614;
    xtab[13-1] =   0.848206583410427216200648320774;
    xtab[14-1] =   0.937273392400705904307758947710;
    xtab[15-1] =   0.987992518020485428489565718587;

    weight[1-1] =  0.307532419961172683546283935772E-01;
    weight[2-1] =  0.703660474881081247092674164507E-01;
    weight[3-1] =  0.107159220467171935011869546686;
    weight[4-1] =  0.139570677926154314447804794511;
    weight[5-1] =  0.166269205816993933553200860481;
    weight[6-1] =  0.186161000015562211026800561866;
    weight[7-1] =  0.198431485327111576456118326444;
    weight[8-1] =  0.202578241925561272880620199968;
    weight[9-1] =  0.198431485327111576456118326444;
    weight[10-1] = 0.186161000015562211026800561866;
    weight[11-1] = 0.166269205816993933553200860481;
    weight[12-1] = 0.139570677926154314447804794511;
    weight[13-1] = 0.107159220467171935011869546686;
    weight[14-1] = 0.703660474881081247092674164507E-01;
    weight[15-1] = 0.307532419961172683546283935772E-01;
  }
  else if ( order == 16 )
  {
    xtab[1-1] =  - 0.989400934991649932596154173450;
    xtab[2-1] =  - 0.944575023073232576077988415535;
    xtab[3-1] =  - 0.865631202387831743880467897712;
    xtab[4-1] =  - 0.755404408355003033895101194847;
    xtab[5-1] =  - 0.617876244402643748446671764049;
    xtab[6-1] =  - 0.458016777657227386342419442984;
    xtab[7-1] =  - 0.281603550779258913230460501460;
    xtab[8-1] =  - 0.950125098376374401853193354250E-01;
    xtab[9-1] =    0.950125098376374401853193354250E-01;
    xtab[10-1] =   0.281603550779258913230460501460;
    xtab[11-1] =   0.458016777657227386342419442984;
    xtab[12-1] =   0.617876244402643748446671764049;
    xtab[13-1] =   0.755404408355003033895101194847;
    xtab[14-1] =   0.865631202387831743880467897712;
    xtab[15-1] =   0.944575023073232576077988415535;
    xtab[16-1] =   0.989400934991649932596154173450;

    weight[1-1] =  0.271524594117540948517805724560E-01;
    weight[2-1] =  0.622535239386478928628438369944E-01;
    weight[3-1] =  0.951585116824927848099251076022E-01;
    weight[4-1] =  0.124628971255533872052476282192;
    weight[5-1] =  0.149595988816576732081501730547;
    weight[6-1] =  0.169156519395002538189312079030;
    weight[7-1] =  0.182603415044923588866763667969;
    weight[8-1] =  0.189450610455068496285396723208;
    weight[9-1] =  0.189450610455068496285396723208;
    weight[10-1] = 0.182603415044923588866763667969;
    weight[11-1] = 0.169156519395002538189312079030;
    weight[12-1] = 0.149595988816576732081501730547;
    weight[13-1] = 0.124628971255533872052476282192;
    weight[14-1] = 0.951585116824927848099251076022E-01;
    weight[15-1] = 0.622535239386478928628438369944E-01;
    weight[16-1] = 0.271524594117540948517805724560E-01;
  }
  else if ( order == 17 )
  {
    xtab[1-1] =  - 0.990575475314417335675434019941;
    xtab[2-1] =  - 0.950675521768767761222716957896;
    xtab[3-1] =  - 0.880239153726985902122955694488;
    xtab[4-1] =  - 0.781514003896801406925230055520;
    xtab[5-1] =  - 0.657671159216690765850302216643;
    xtab[6-1] =  - 0.512690537086476967886246568630;
    xtab[7-1] =  - 0.351231763453876315297185517095;
    xtab[8-1] =  - 0.178484181495847855850677493654;
    xtab[9-1] =    0.0;
    xtab[10-1] =   0.178484181495847855850677493654;
    xtab[11-1] =   0.351231763453876315297185517095;
    xtab[12-1] =   0.512690537086476967886246568630;
    xtab[13-1] =   0.657671159216690765850302216643;
    xtab[14-1] =   0.781514003896801406925230055520;
    xtab[15-1] =   0.880239153726985902122955694488;
    xtab[16-1] =   0.950675521768767761222716957896;
    xtab[17-1] =   0.990575475314417335675434019941;

    weight[1-1] =  0.241483028685479319601100262876E-01;
    weight[2-1] =  0.554595293739872011294401653582E-01;
    weight[3-1] =  0.850361483171791808835353701911E-01;
    weight[4-1] =  0.111883847193403971094788385626;
    weight[5-1] =  0.135136368468525473286319981702;
    weight[6-1] =  0.154045761076810288081431594802;
    weight[7-1] =  0.168004102156450044509970663788;
    weight[8-1] =  0.176562705366992646325270990113;
    weight[9-1] =  0.179446470356206525458265644262;
    weight[10-1] = 0.176562705366992646325270990113;
    weight[11-1] = 0.168004102156450044509970663788;
    weight[12-1] = 0.154045761076810288081431594802;
    weight[13-1] = 0.135136368468525473286319981702;
    weight[14-1] = 0.111883847193403971094788385626;
    weight[15-1] = 0.850361483171791808835353701911E-01;
    weight[16-1] = 0.554595293739872011294401653582E-01;
    weight[17-1] = 0.241483028685479319601100262876E-01;
  }
  else if ( order == 18 )
  {
    xtab[1-1] =  - 0.991565168420930946730016004706;
    xtab[2-1] =  - 0.955823949571397755181195892930;
    xtab[3-1] =  - 0.892602466497555739206060591127;
    xtab[4-1] =  - 0.803704958972523115682417455015;
    xtab[5-1] =  - 0.691687043060353207874891081289;
    xtab[6-1] =  - 0.559770831073947534607871548525;
    xtab[7-1] =  - 0.411751161462842646035931793833;
    xtab[8-1] =  - 0.251886225691505509588972854878;
    xtab[9-1] =  - 0.847750130417353012422618529358E-01;
    xtab[10-1] =   0.847750130417353012422618529358E-01;
    xtab[11-1] =   0.251886225691505509588972854878;
    xtab[12-1] =   0.411751161462842646035931793833;
    xtab[13-1] =   0.559770831073947534607871548525;
    xtab[14-1] =   0.691687043060353207874891081289;
    xtab[15-1] =   0.803704958972523115682417455015;
    xtab[16-1] =   0.892602466497555739206060591127;
    xtab[17-1] =   0.955823949571397755181195892930;
    xtab[18-1] =   0.991565168420930946730016004706;

    weight[1-1] =  0.216160135264833103133427102665E-01;
    weight[2-1] =  0.497145488949697964533349462026E-01;
    weight[3-1] =  0.764257302548890565291296776166E-01;
    weight[4-1] =  0.100942044106287165562813984925;
    weight[5-1] =  0.122555206711478460184519126800;
    weight[6-1] =  0.140642914670650651204731303752;
    weight[7-1] =  0.154684675126265244925418003836;
    weight[8-1] =  0.164276483745832722986053776466;
    weight[9-1] =  0.169142382963143591840656470135;
    weight[10-1] = 0.169142382963143591840656470135;
    weight[11-1] = 0.164276483745832722986053776466;
    weight[12-1] = 0.154684675126265244925418003836;
    weight[13-1] = 0.140642914670650651204731303752;
    weight[14-1] = 0.122555206711478460184519126800;
    weight[15-1] = 0.100942044106287165562813984925;
    weight[16-1] = 0.764257302548890565291296776166E-01;
    weight[17-1] = 0.497145488949697964533349462026E-01;
    weight[18-1] = 0.216160135264833103133427102665E-01;
  }
  else if ( order == 19 )
  {
    xtab[1-1] =  - 0.992406843843584403189017670253;
    xtab[2-1] =  - 0.960208152134830030852778840688;
    xtab[3-1] =  - 0.903155903614817901642660928532;
    xtab[4-1] =  - 0.822714656537142824978922486713;
    xtab[5-1] =  - 0.720966177335229378617095860824;
    xtab[6-1] =  - 0.600545304661681023469638164946;
    xtab[7-1] =  - 0.464570741375960945717267148104;
    xtab[8-1] =  - 0.316564099963629831990117328850;
    xtab[9-1] =  - 0.160358645640225375868096115741;
    xtab[10-1] =   0.0;
    xtab[11-1] =   0.160358645640225375868096115741;
    xtab[12-1] =   0.316564099963629831990117328850;
    xtab[13-1] =   0.464570741375960945717267148104;
    xtab[14-1] =   0.600545304661681023469638164946;
    xtab[15-1] =   0.720966177335229378617095860824;
    xtab[16-1] =   0.822714656537142824978922486713;
    xtab[17-1] =   0.903155903614817901642660928532;
    xtab[18-1] =   0.960208152134830030852778840688;
    xtab[19-1] =   0.992406843843584403189017670253;

    weight[1-1] =  0.194617882297264770363120414644E-01;
    weight[2-1] =  0.448142267656996003328381574020E-01;
    weight[3-1] =  0.690445427376412265807082580060E-01;
    weight[4-1] =  0.914900216224499994644620941238E-01;
    weight[5-1] =  0.111566645547333994716023901682;
    weight[6-1] =  0.128753962539336227675515784857;
    weight[7-1] =  0.142606702173606611775746109442;
    weight[8-1] =  0.152766042065859666778855400898;
    weight[9-1] =  0.158968843393954347649956439465;
    weight[10-1] = 0.161054449848783695979163625321;
    weight[11-1] = 0.158968843393954347649956439465;
    weight[12-1] = 0.152766042065859666778855400898;
    weight[13-1] = 0.142606702173606611775746109442;
    weight[14-1] = 0.128753962539336227675515784857;
    weight[15-1] = 0.111566645547333994716023901682;
    weight[16-1] = 0.914900216224499994644620941238E-01;
    weight[17-1] = 0.690445427376412265807082580060E-01;
    weight[18-1] = 0.448142267656996003328381574020E-01;
    weight[19-1] = 0.194617882297264770363120414644E-01;
  }
  else if ( order == 20 )
  {
    xtab[1-1] =  - 0.993128599185094924786122388471;
    xtab[2-1] =  - 0.963971927277913791267666131197;
    xtab[3-1] =  - 0.912234428251325905867752441203;
    xtab[4-1] =  - 0.839116971822218823394529061702;
    xtab[5-1] =  - 0.746331906460150792614305070356;
    xtab[6-1] =  - 0.636053680726515025452836696226;
    xtab[7-1] =  - 0.510867001950827098004364050955;
    xtab[8-1] =  - 0.373706088715419560672548177025;
    xtab[9-1] =  - 0.227785851141645078080496195369;
    xtab[10-1] = - 0.765265211334973337546404093988E-01;
    xtab[11-1] =   0.765265211334973337546404093988E-01;
    xtab[12-1] =   0.227785851141645078080496195369;
    xtab[13-1] =   0.373706088715419560672548177025;
    xtab[14-1] =   0.510867001950827098004364050955;
    xtab[15-1] =   0.636053680726515025452836696226;
    xtab[16-1] =   0.746331906460150792614305070356;
    xtab[17-1] =   0.839116971822218823394529061702;
    xtab[18-1] =   0.912234428251325905867752441203;
    xtab[19-1] =   0.963971927277913791267666131197;
    xtab[20-1] =   0.993128599185094924786122388471;

    weight[1-1] =  0.176140071391521183118619623519E-01;
    weight[2-1] =  0.406014298003869413310399522749E-01;
    weight[3-1] =  0.626720483341090635695065351870E-01;
    weight[4-1] =  0.832767415767047487247581432220E-01;
    weight[5-1] =  0.101930119817240435036750135480;
    weight[6-1] =  0.118194531961518417312377377711;
    weight[7-1] =  0.131688638449176626898494499748;
    weight[8-1] =  0.142096109318382051329298325067;
    weight[9-1] =  0.149172986472603746787828737002;
    weight[10-1] = 0.152753387130725850698084331955;
    weight[11-1] = 0.152753387130725850698084331955;
    weight[12-1] = 0.149172986472603746787828737002;
    weight[13-1] = 0.142096109318382051329298325067;
    weight[14-1] = 0.131688638449176626898494499748;
    weight[15-1] = 0.118194531961518417312377377711;
    weight[16-1] = 0.101930119817240435036750135480;
    weight[17-1] = 0.832767415767047487247581432220E-01;
    weight[18-1] = 0.626720483341090635695065351870E-01;
    weight[19-1] = 0.406014298003869413310399522749E-01;
    weight[20-1] = 0.176140071391521183118619623519E-01;
  }
  else if ( order == 21 )
  {
    xtab[ 1-1] =  -0.9937521706203896E+00;
    xtab[ 2-1] =  -0.9672268385663063E+00;
    xtab[ 3-1] =  -0.9200993341504008E+00;
    xtab[ 4-1] =  -0.8533633645833173E+00;
    xtab[ 5-1] =  -0.7684399634756779E+00;
    xtab[ 6-1] =  -0.6671388041974123E+00;
    xtab[ 7-1] =  -0.5516188358872198E+00;
    xtab[ 8-1] =  -0.4243421202074388E+00;
    xtab[ 9-1] =  -0.2880213168024011E+00;
    xtab[10-1] =  -0.1455618541608951E+00;
    xtab[11-1] =   0.0000000000000000E+00;
    xtab[12-1] =   0.1455618541608951E+00;
    xtab[13-1] =   0.2880213168024011E+00;
    xtab[14-1] =   0.4243421202074388E+00;
    xtab[15-1] =   0.5516188358872198E+00;
    xtab[16-1] =   0.6671388041974123E+00;
    xtab[17-1] =   0.7684399634756779E+00;
    xtab[18-1] =   0.8533633645833173E+00;
    xtab[19-1] =   0.9200993341504008E+00;
    xtab[20-1] =   0.9672268385663063E+00;
    xtab[21-1] =   0.9937521706203896E+00;

    weight[ 1-1] =   0.1601722825777420E-01;
    weight[ 2-1] =   0.3695378977085242E-01;
    weight[ 3-1] =   0.5713442542685715E-01;
    weight[ 4-1] =   0.7610011362837928E-01;
    weight[ 5-1] =   0.9344442345603393E-01;
    weight[ 6-1] =   0.1087972991671484E+00;
    weight[ 7-1] =   0.1218314160537285E+00;
    weight[ 8-1] =   0.1322689386333373E+00;
    weight[ 9-1] =   0.1398873947910731E+00;
    weight[10-1] =   0.1445244039899700E+00;
    weight[11-1] =   0.1460811336496904E+00;
    weight[12-1] =   0.1445244039899700E+00;
    weight[13-1] =   0.1398873947910731E+00;
    weight[14-1] =   0.1322689386333373E+00;
    weight[15-1] =   0.1218314160537285E+00;
    weight[16-1] =   0.1087972991671484E+00;
    weight[17-1] =   0.9344442345603393E-01;
    weight[18-1] =   0.7610011362837928E-01;
    weight[19-1] =   0.5713442542685715E-01;
    weight[20-1] =   0.3695378977085242E-01;
    weight[21-1] =   0.1601722825777420E-01;
  }
  else if ( order == 22 )
  {
    xtab[ 1-1] =  -0.9942945854823994E+00;
    xtab[ 2-1] =  -0.9700604978354287E+00;
    xtab[ 3-1] =  -0.9269567721871740E+00;
    xtab[ 4-1] =  -0.8658125777203002E+00;
    xtab[ 5-1] =  -0.7878168059792081E+00;
    xtab[ 6-1] =  -0.6944872631866827E+00;
    xtab[ 7-1] =  -0.5876404035069116E+00;
    xtab[ 8-1] =  -0.4693558379867570E+00;
    xtab[ 9-1] =  -0.3419358208920842E+00;
    xtab[10-1] =  -0.2078604266882213E+00;
    xtab[11-1] =  -0.6973927331972223E-01;
    xtab[12-1] =   0.6973927331972223E-01;
    xtab[13-1] =   0.2078604266882213E+00;
    xtab[14-1] =   0.3419358208920842E+00;
    xtab[15-1] =   0.4693558379867570E+00;
    xtab[16-1] =   0.5876404035069116E+00;
    xtab[17-1] =   0.6944872631866827E+00;
    xtab[18-1] =   0.7878168059792081E+00;
    xtab[19-1] =   0.8658125777203002E+00;
    xtab[20-1] =   0.9269567721871740E+00;
    xtab[21-1] =   0.9700604978354287E+00;
    xtab[22-1] =   0.9942945854823994E+00;

    weight[ 1-1] =   0.1462799529827203E-01;
    weight[ 2-1] =   0.3377490158481413E-01;
    weight[ 3-1] =   0.5229333515268327E-01;
    weight[ 4-1] =   0.6979646842452038E-01;
    weight[ 5-1] =   0.8594160621706777E-01;
    weight[ 6-1] =   0.1004141444428809E+00;
    weight[ 7-1] =   0.1129322960805392E+00;
    weight[ 8-1] =   0.1232523768105124E+00;
    weight[ 9-1] =   0.1311735047870623E+00;
    weight[10-1] =   0.1365414983460152E+00;
    weight[11-1] =   0.1392518728556321E+00;
    weight[12-1] =   0.1392518728556321E+00;
    weight[13-1] =   0.1365414983460152E+00;
    weight[14-1] =   0.1311735047870623E+00;
    weight[15-1] =   0.1232523768105124E+00;
    weight[16-1] =   0.1129322960805392E+00;
    weight[17-1] =   0.1004141444428809E+00;
    weight[18-1] =   0.8594160621706777E-01;
    weight[19-1] =   0.6979646842452038E-01;
    weight[20-1] =   0.5229333515268327E-01;
    weight[21-1] =   0.3377490158481413E-01;
    weight[22-1] =   0.1462799529827203E-01;
  }
  else if ( order == 23 )
  {
    xtab[ 1-1] =  -0.9947693349975522E+00;
    xtab[ 2-1] =  -0.9725424712181152E+00;
    xtab[ 3-1] =  -0.9329710868260161E+00;
    xtab[ 4-1] =  -0.8767523582704416E+00;
    xtab[ 5-1] =  -0.8048884016188399E+00;
    xtab[ 6-1] =  -0.7186613631319502E+00;
    xtab[ 7-1] =  -0.6196098757636461E+00;
    xtab[ 8-1] =  -0.5095014778460075E+00;
    xtab[ 9-1] =  -0.3903010380302908E+00;
    xtab[10-1] =  -0.2641356809703449E+00;
    xtab[11-1] =  -0.1332568242984661E+00;
    xtab[12-1] =   0.0000000000000000E+00;
    xtab[13-1] =   0.1332568242984661E+00;
    xtab[14-1] =   0.2641356809703449E+00;
    xtab[15-1] =   0.3903010380302908E+00;
    xtab[16-1] =   0.5095014778460075E+00;
    xtab[17-1] =   0.6196098757636461E+00;
    xtab[18-1] =   0.7186613631319502E+00;
    xtab[19-1] =   0.8048884016188399E+00;
    xtab[20-1] =   0.8767523582704416E+00;
    xtab[21-1] =   0.9329710868260161E+00;
    xtab[22-1] =   0.9725424712181152E+00;
    xtab[23-1] =   0.9947693349975522E+00;

    weight[ 1-1] =   0.1341185948714167E-01;
    weight[ 2-1] =   0.3098800585697944E-01;
    weight[ 3-1] =   0.4803767173108464E-01;
    weight[ 4-1] =   0.6423242140852586E-01;
    weight[ 5-1] =   0.7928141177671895E-01;
    weight[ 6-1] =   0.9291576606003514E-01;
    weight[ 7-1] =   0.1048920914645414E+00;
    weight[ 8-1] =   0.1149966402224114E+00;
    weight[ 9-1] =   0.1230490843067295E+00;
    weight[10-1] =   0.1289057221880822E+00;
    weight[11-1] =   0.1324620394046967E+00;
    weight[12-1] =   0.1336545721861062E+00;
    weight[13-1] =   0.1324620394046967E+00;
    weight[14-1] =   0.1289057221880822E+00;
    weight[15-1] =   0.1230490843067295E+00;
    weight[16-1] =   0.1149966402224114E+00;
    weight[17-1] =   0.1048920914645414E+00;
    weight[18-1] =   0.9291576606003514E-01;
    weight[19-1] =   0.7928141177671895E-01;
    weight[20-1] =   0.6423242140852586E-01;
    weight[21-1] =   0.4803767173108464E-01;
    weight[22-1] =   0.3098800585697944E-01;
    weight[23-1] =   0.1341185948714167E-01;
  }
  else if ( order == 24 )
  {
    xtab[ 1-1] =  -0.9951872199970213E+00;
    xtab[ 2-1] =  -0.9747285559713095E+00;
    xtab[ 3-1] =  -0.9382745520027327E+00;
    xtab[ 4-1] =  -0.8864155270044011E+00;
    xtab[ 5-1] =  -0.8200019859739029E+00;
    xtab[ 6-1] =  -0.7401241915785544E+00;
    xtab[ 7-1] =  -0.6480936519369755E+00;
    xtab[ 8-1] =  -0.5454214713888396E+00;
    xtab[ 9-1] =  -0.4337935076260451E+00;
    xtab[10-1] =  -0.3150426796961634E+00;
    xtab[11-1] =  -0.1911188674736163E+00;
    xtab[12-1] =  -0.6405689286260562E-01;
    xtab[13-1] =   0.6405689286260562E-01;
    xtab[14-1] =   0.1911188674736163E+00;
    xtab[15-1] =   0.3150426796961634E+00;
    xtab[16-1] =   0.4337935076260451E+00;
    xtab[17-1] =   0.5454214713888396E+00;
    xtab[18-1] =   0.6480936519369755E+00;
    xtab[19-1] =   0.7401241915785544E+00;
    xtab[20-1] =   0.8200019859739029E+00;
    xtab[21-1] =   0.8864155270044011E+00;
    xtab[22-1] =   0.9382745520027327E+00;
    xtab[23-1] =   0.9747285559713095E+00;
    xtab[24-1] =   0.9951872199970213E+00;

    weight[ 1-1] =   0.1234122979998730E-01;
    weight[ 2-1] =   0.2853138862893375E-01;
    weight[ 3-1] =   0.4427743881741982E-01;
    weight[ 4-1] =   0.5929858491543672E-01;
    weight[ 5-1] =   0.7334648141108031E-01;
    weight[ 6-1] =   0.8619016153195320E-01;
    weight[ 7-1] =   0.9761865210411380E-01;
    weight[ 8-1] =   0.1074442701159656E+00;
    weight[ 9-1] =   0.1155056680537256E+00;
    weight[10-1] =   0.1216704729278035E+00;
    weight[11-1] =   0.1258374563468283E+00;
    weight[12-1] =   0.1279381953467521E+00;
    weight[13-1] =   0.1279381953467521E+00;
    weight[14-1] =   0.1258374563468283E+00;
    weight[15-1] =   0.1216704729278035E+00;
    weight[16-1] =   0.1155056680537256E+00;
    weight[17-1] =   0.1074442701159656E+00;
    weight[18-1] =   0.9761865210411380E-01;
    weight[19-1] =   0.8619016153195320E-01;
    weight[20-1] =   0.7334648141108031E-01;
    weight[21-1] =   0.5929858491543672E-01;
    weight[22-1] =   0.4427743881741982E-01;
    weight[23-1] =   0.2853138862893375E-01;
    weight[24-1] =   0.1234122979998730E-01;
  }
  else if ( order == 25 )
  {
    xtab[ 1-1] =  -0.9955569697904981E+00;
    xtab[ 2-1] =  -0.9766639214595175E+00;
    xtab[ 3-1] =  -0.9429745712289743E+00;
    xtab[ 4-1] =  -0.8949919978782754E+00;
    xtab[ 5-1] =  -0.8334426287608340E+00;
    xtab[ 6-1] =  -0.7592592630373577E+00;
    xtab[ 7-1] =  -0.6735663684734684E+00;
    xtab[ 8-1] =  -0.5776629302412229E+00;
    xtab[ 9-1] =  -0.4730027314457150E+00;
    xtab[10-1] =  -0.3611723058093879E+00;
    xtab[11-1] =  -0.2438668837209884E+00;
    xtab[12-1] =  -0.1228646926107104E+00;
    xtab[13-1] =   0.0000000000000000E+00;
    xtab[14-1] =   0.1228646926107104E+00;
    xtab[15-1] =   0.2438668837209884E+00;
    xtab[16-1] =   0.3611723058093879E+00;
    xtab[17-1] =   0.4730027314457150E+00;
    xtab[18-1] =   0.5776629302412229E+00;
    xtab[19-1] =   0.6735663684734684E+00;
    xtab[20-1] =   0.7592592630373577E+00;
    xtab[21-1] =   0.8334426287608340E+00;
    xtab[22-1] =   0.8949919978782754E+00;
    xtab[23-1] =   0.9429745712289743E+00;
    xtab[24-1] =   0.9766639214595175E+00;
    xtab[25-1] =   0.9955569697904981E+00;

    weight[ 1-1] =   0.1139379850102617E-01;
    weight[ 2-1] =   0.2635498661503214E-01;
    weight[ 3-1] =   0.4093915670130639E-01;
    weight[ 4-1] =   0.5490469597583517E-01;
    weight[ 5-1] =   0.6803833381235694E-01;
    weight[ 6-1] =   0.8014070033500101E-01;
    weight[ 7-1] =   0.9102826198296370E-01;
    weight[ 8-1] =   0.1005359490670506E+00;
    weight[ 9-1] =   0.1085196244742637E+00;
    weight[10-1] =   0.1148582591457116E+00;
    weight[11-1] =   0.1194557635357847E+00;
    weight[12-1] =   0.1222424429903101E+00;
    weight[13-1] =   0.1231760537267154E+00;
    weight[14-1] =   0.1222424429903101E+00;
    weight[15-1] =   0.1194557635357847E+00;
    weight[16-1] =   0.1148582591457116E+00;
    weight[17-1] =   0.1085196244742637E+00;
    weight[18-1] =   0.1005359490670506E+00;
    weight[19-1] =   0.9102826198296370E-01;
    weight[20-1] =   0.8014070033500101E-01;
    weight[21-1] =   0.6803833381235694E-01;
    weight[22-1] =   0.5490469597583517E-01;
    weight[23-1] =   0.4093915670130639E-01;
    weight[24-1] =   0.2635498661503214E-01;
    weight[25-1] =   0.1139379850102617E-01;
  }
  else if ( order == 26 )
  {
    xtab[ 1-1] =  -0.9958857011456169E+00;
    xtab[ 2-1] =  -0.9783854459564710E+00;
    xtab[ 3-1] =  -0.9471590666617142E+00;
    xtab[ 4-1] =  -0.9026378619843071E+00;
    xtab[ 5-1] =  -0.8454459427884981E+00;
    xtab[ 6-1] =  -0.7763859488206789E+00;
    xtab[ 7-1] =  -0.6964272604199573E+00;
    xtab[ 8-1] =  -0.6066922930176181E+00;
    xtab[ 9-1] =  -0.5084407148245057E+00;
    xtab[10-1] =  -0.4030517551234863E+00;
    xtab[11-1] =  -0.2920048394859569E+00;
    xtab[12-1] =  -0.1768588203568902E+00;
    xtab[13-1] =  -0.5923009342931320E-01;
    xtab[14-1] =   0.5923009342931320E-01;
    xtab[15-1] =   0.1768588203568902E+00;
    xtab[16-1] =   0.2920048394859569E+00;
    xtab[17-1] =   0.4030517551234863E+00;
    xtab[18-1] =   0.5084407148245057E+00;
    xtab[19-1] =   0.6066922930176181E+00;
    xtab[20-1] =   0.6964272604199573E+00;
    xtab[21-1] =   0.7763859488206789E+00;
    xtab[22-1] =   0.8454459427884981E+00;
    xtab[23-1] =   0.9026378619843071E+00;
    xtab[24-1] =   0.9471590666617142E+00;
    xtab[25-1] =   0.9783854459564710E+00;
    xtab[26-1] =   0.9958857011456169E+00;

    weight[ 1-1] =   0.1055137261734304E-01;
    weight[ 2-1] =   0.2441785109263173E-01;
    weight[ 3-1] =   0.3796238329436282E-01;
    weight[ 4-1] =   0.5097582529714782E-01;
    weight[ 5-1] =   0.6327404632957484E-01;
    weight[ 6-1] =   0.7468414976565967E-01;
    weight[ 7-1] =   0.8504589431348521E-01;
    weight[ 8-1] =   0.9421380035591416E-01;
    weight[ 9-1] =   0.1020591610944255E+00;
    weight[10-1] =   0.1084718405285765E+00;
    weight[11-1] =   0.1133618165463197E+00;
    weight[12-1] =   0.1166604434852967E+00;
    weight[13-1] =   0.1183214152792622E+00;
    weight[14-1] =   0.1183214152792622E+00;
    weight[15-1] =   0.1166604434852967E+00;
    weight[16-1] =   0.1133618165463197E+00;
    weight[17-1] =   0.1084718405285765E+00;
    weight[18-1] =   0.1020591610944255E+00;
    weight[19-1] =   0.9421380035591416E-01;
    weight[20-1] =   0.8504589431348521E-01;
    weight[21-1] =   0.7468414976565967E-01;
    weight[22-1] =   0.6327404632957484E-01;
    weight[23-1] =   0.5097582529714782E-01;
    weight[24-1] =   0.3796238329436282E-01;
    weight[25-1] =   0.2441785109263173E-01;
    weight[26-1] =   0.1055137261734304E-01;
  }
  else if ( order == 27 )
  {
    xtab[ 1-1] =  -0.9961792628889886E+00;
    xtab[ 2-1] =  -0.9799234759615012E+00;
    xtab[ 3-1] =  -0.9509005578147051E+00;
    xtab[ 4-1] =  -0.9094823206774911E+00;
    xtab[ 5-1] =  -0.8562079080182945E+00;
    xtab[ 6-1] =  -0.7917716390705082E+00;
    xtab[ 7-1] =  -0.7170134737394237E+00;
    xtab[ 8-1] =  -0.6329079719464952E+00;
    xtab[ 9-1] =  -0.5405515645794569E+00;
    xtab[10-1] =  -0.4411482517500269E+00;
    xtab[11-1] =  -0.3359939036385089E+00;
    xtab[12-1] =  -0.2264593654395369E+00;
    xtab[13-1] =  -0.1139725856095300E+00;
    xtab[14-1] =   0.0000000000000000E+00;
    xtab[15-1] =   0.1139725856095300E+00;
    xtab[16-1] =   0.2264593654395369E+00;
    xtab[17-1] =   0.3359939036385089E+00;
    xtab[18-1] =   0.4411482517500269E+00;
    xtab[19-1] =   0.5405515645794569E+00;
    xtab[20-1] =   0.6329079719464952E+00;
    xtab[21-1] =   0.7170134737394237E+00;
    xtab[22-1] =   0.7917716390705082E+00;
    xtab[23-1] =   0.8562079080182945E+00;
    xtab[24-1] =   0.9094823206774911E+00;
    xtab[25-1] =   0.9509005578147051E+00;
    xtab[26-1] =   0.9799234759615012E+00;
    xtab[27-1] =   0.9961792628889886E+00;

    weight[ 1-1] =   0.9798996051294232E-02;
    weight[ 2-1] =   0.2268623159618062E-01;
    weight[ 3-1] =   0.3529705375741969E-01;
    weight[ 4-1] =   0.4744941252061504E-01;
    weight[ 5-1] =   0.5898353685983366E-01;
    weight[ 6-1] =   0.6974882376624561E-01;
    weight[ 7-1] =   0.7960486777305781E-01;
    weight[ 8-1] =   0.8842315854375689E-01;
    weight[ 9-1] =   0.9608872737002842E-01;
    weight[10-1] =   0.1025016378177459E+00;
    weight[11-1] =   0.1075782857885332E+00;
    weight[12-1] =   0.1112524883568452E+00;
    weight[13-1] =   0.1134763461089651E+00;
    weight[14-1] =   0.1142208673789570E+00;
    weight[15-1] =   0.1134763461089651E+00;
    weight[16-1] =   0.1112524883568452E+00;
    weight[17-1] =   0.1075782857885332E+00;
    weight[18-1] =   0.1025016378177459E+00;
    weight[19-1] =   0.9608872737002842E-01;
    weight[20-1] =   0.8842315854375689E-01;
    weight[21-1] =   0.7960486777305781E-01;
    weight[22-1] =   0.6974882376624561E-01;
    weight[23-1] =   0.5898353685983366E-01;
    weight[24-1] =   0.4744941252061504E-01;
    weight[25-1] =   0.3529705375741969E-01;
    weight[26-1] =   0.2268623159618062E-01;
    weight[27-1] =   0.9798996051294232E-02;
  }
  else if ( order == 28 )
  {
    xtab[ 1-1] =  -0.9964424975739544E+00;
    xtab[ 2-1] =  -0.9813031653708728E+00;
    xtab[ 3-1] =  -0.9542592806289382E+00;
    xtab[ 4-1] =  -0.9156330263921321E+00;
    xtab[ 5-1] =  -0.8658925225743951E+00;
    xtab[ 6-1] =  -0.8056413709171791E+00;
    xtab[ 7-1] =  -0.7356108780136318E+00;
    xtab[ 8-1] =  -0.6566510940388650E+00;
    xtab[ 9-1] =  -0.5697204718114017E+00;
    xtab[10-1] =  -0.4758742249551183E+00;
    xtab[11-1] =  -0.3762515160890787E+00;
    xtab[12-1] =  -0.2720616276351780E+00;
    xtab[13-1] =  -0.1645692821333808E+00;
    xtab[14-1] =  -0.5507928988403427E-01;
    xtab[15-1] =   0.5507928988403427E-01;
    xtab[16-1] =   0.1645692821333808E+00;
    xtab[17-1] =   0.2720616276351780E+00;
    xtab[18-1] =   0.3762515160890787E+00;
    xtab[19-1] =   0.4758742249551183E+00;
    xtab[20-1] =   0.5697204718114017E+00;
    xtab[21-1] =   0.6566510940388650E+00;
    xtab[22-1] =   0.7356108780136318E+00;
    xtab[23-1] =   0.8056413709171791E+00;
    xtab[24-1] =   0.8658925225743951E+00;
    xtab[25-1] =   0.9156330263921321E+00;
    xtab[26-1] =   0.9542592806289382E+00;
    xtab[27-1] =   0.9813031653708728E+00;
    xtab[28-1] =   0.9964424975739544E+00;

    weight[ 1-1] =   0.9124282593094672E-02;
    weight[ 2-1] =   0.2113211259277118E-01;
    weight[ 3-1] =   0.3290142778230441E-01;
    weight[ 4-1] =   0.4427293475900429E-01;
    weight[ 5-1] =   0.5510734567571667E-01;
    weight[ 6-1] =   0.6527292396699959E-01;
    weight[ 7-1] =   0.7464621423456877E-01;
    weight[ 8-1] =   0.8311341722890127E-01;
    weight[ 9-1] =   0.9057174439303289E-01;
    weight[10-1] =   0.9693065799792999E-01;
    weight[11-1] =   0.1021129675780608E+00;
    weight[12-1] =   0.1060557659228464E+00;
    weight[13-1] =   0.1087111922582942E+00;
    weight[14-1] =   0.1100470130164752E+00;
    weight[15-1] =   0.1100470130164752E+00;
    weight[16-1] =   0.1087111922582942E+00;
    weight[17-1] =   0.1060557659228464E+00;
    weight[18-1] =   0.1021129675780608E+00;
    weight[19-1] =   0.9693065799792999E-01;
    weight[20-1] =   0.9057174439303289E-01;
    weight[21-1] =   0.8311341722890127E-01;
    weight[22-1] =   0.7464621423456877E-01;
    weight[23-1] =   0.6527292396699959E-01;
    weight[24-1] =   0.5510734567571667E-01;
    weight[25-1] =   0.4427293475900429E-01;
    weight[26-1] =   0.3290142778230441E-01;
    weight[27-1] =   0.2113211259277118E-01;
    weight[28-1] =   0.9124282593094672E-02;
  }
  else if ( order == 29 )
  {
    xtab[ 1-1] =  -0.9966794422605966E+00;
    xtab[ 2-1] =  -0.9825455052614132E+00;
    xtab[ 3-1] =  -0.9572855957780877E+00;
    xtab[ 4-1] =  -0.9211802329530588E+00;
    xtab[ 5-1] =  -0.8746378049201028E+00;
    xtab[ 6-1] =  -0.8181854876152524E+00;
    xtab[ 7-1] =  -0.7524628517344771E+00;
    xtab[ 8-1] =  -0.6782145376026865E+00;
    xtab[ 9-1] =  -0.5962817971382278E+00;
    xtab[10-1] =  -0.5075929551242276E+00;
    xtab[11-1] =  -0.4131528881740087E+00;
    xtab[12-1] =  -0.3140316378676399E+00;
    xtab[13-1] =  -0.2113522861660011E+00;
    xtab[14-1] =  -0.1062782301326792E+00;
    xtab[15-1] =   0.0000000000000000E+00;
    xtab[16-1] =   0.1062782301326792E+00;
    xtab[17-1] =   0.2113522861660011E+00;
    xtab[18-1] =   0.3140316378676399E+00;
    xtab[19-1] =   0.4131528881740087E+00;
    xtab[20-1] =   0.5075929551242276E+00;
    xtab[21-1] =   0.5962817971382278E+00;
    xtab[22-1] =   0.6782145376026865E+00;
    xtab[23-1] =   0.7524628517344771E+00;
    xtab[24-1] =   0.8181854876152524E+00;
    xtab[25-1] =   0.8746378049201028E+00;
    xtab[26-1] =   0.9211802329530588E+00;
    xtab[27-1] =   0.9572855957780877E+00;
    xtab[28-1] =   0.9825455052614132E+00;
    xtab[29-1] =   0.9966794422605966E+00;

    weight[ 1-1] =   0.8516903878746365E-02;
    weight[ 2-1] =   0.1973208505612276E-01;
    weight[ 3-1] =   0.3074049220209360E-01;
    weight[ 4-1] =   0.4140206251868281E-01;
    weight[ 5-1] =   0.5159482690249799E-01;
    weight[ 6-1] =   0.6120309065707916E-01;
    weight[ 7-1] =   0.7011793325505125E-01;
    weight[ 8-1] =   0.7823832713576385E-01;
    weight[ 9-1] =   0.8547225736617248E-01;
    weight[10-1] =   0.9173775713925882E-01;
    weight[11-1] =   0.9696383409440862E-01;
    weight[12-1] =   0.1010912737599150E+00;
    weight[13-1] =   0.1040733100777293E+00;
    weight[14-1] =   0.1058761550973210E+00;
    weight[15-1] =   0.1064793817183143E+00;
    weight[16-1] =   0.1058761550973210E+00;
    weight[17-1] =   0.1040733100777293E+00;
    weight[18-1] =   0.1010912737599150E+00;
    weight[19-1] =   0.9696383409440862E-01;
    weight[20-1] =   0.9173775713925882E-01;
    weight[21-1] =   0.8547225736617248E-01;
    weight[22-1] =   0.7823832713576385E-01;
    weight[23-1] =   0.7011793325505125E-01;
    weight[24-1] =   0.6120309065707916E-01;
    weight[25-1] =   0.5159482690249799E-01;
    weight[26-1] =   0.4140206251868281E-01;
    weight[27-1] =   0.3074049220209360E-01;
    weight[28-1] =   0.1973208505612276E-01;
    weight[29-1] =   0.8516903878746365E-02;
  }
  else if ( order == 30 )
  {
    xtab[ 1-1] =  -0.9968934840746495E+00;
    xtab[ 2-1] =  -0.9836681232797472E+00;
    xtab[ 3-1] =  -0.9600218649683075E+00;
    xtab[ 4-1] =  -0.9262000474292743E+00;
    xtab[ 5-1] =  -0.8825605357920526E+00;
    xtab[ 6-1] =  -0.8295657623827684E+00;
    xtab[ 7-1] =  -0.7677774321048262E+00;
    xtab[ 8-1] =  -0.6978504947933158E+00;
    xtab[ 9-1] =  -0.6205261829892429E+00;
    xtab[10-1] =  -0.5366241481420199E+00;
    xtab[11-1] =  -0.4470337695380892E+00;
    xtab[12-1] =  -0.3527047255308781E+00;
    xtab[13-1] =  -0.2546369261678899E+00;
    xtab[14-1] =  -0.1538699136085835E+00;
    xtab[15-1] =  -0.5147184255531770E-01;
    xtab[16-1] =   0.5147184255531770E-01;
    xtab[17-1] =   0.1538699136085835E+00;
    xtab[18-1] =   0.2546369261678899E+00;
    xtab[19-1] =   0.3527047255308781E+00;
    xtab[20-1] =   0.4470337695380892E+00;
    xtab[21-1] =   0.5366241481420199E+00;
    xtab[22-1] =   0.6205261829892429E+00;
    xtab[23-1] =   0.6978504947933158E+00;
    xtab[24-1] =   0.7677774321048262E+00;
    xtab[25-1] =   0.8295657623827684E+00;
    xtab[26-1] =   0.8825605357920526E+00;
    xtab[27-1] =   0.9262000474292743E+00;
    xtab[28-1] =   0.9600218649683075E+00;
    xtab[29-1] =   0.9836681232797472E+00;
    xtab[30-1] =   0.9968934840746495E+00;

    weight[ 1-1] =   0.7968192496166648E-02;
    weight[ 2-1] =   0.1846646831109099E-01;
    weight[ 3-1] =   0.2878470788332330E-01;
    weight[ 4-1] =   0.3879919256962704E-01;
    weight[ 5-1] =   0.4840267283059405E-01;
    weight[ 6-1] =   0.5749315621761905E-01;
    weight[ 7-1] =   0.6597422988218052E-01;
    weight[ 8-1] =   0.7375597473770516E-01;
    weight[ 9-1] =   0.8075589522942023E-01;
    weight[10-1] =   0.8689978720108314E-01;
    weight[11-1] =   0.9212252223778619E-01;
    weight[12-1] =   0.9636873717464424E-01;
    weight[13-1] =   0.9959342058679524E-01;
    weight[14-1] =   0.1017623897484056E+00;
    weight[15-1] =   0.1028526528935587E+00;
    weight[16-1] =   0.1028526528935587E+00;
    weight[17-1] =   0.1017623897484056E+00;
    weight[18-1] =   0.9959342058679524E-01;
    weight[19-1] =   0.9636873717464424E-01;
    weight[20-1] =   0.9212252223778619E-01;
    weight[21-1] =   0.8689978720108314E-01;
    weight[22-1] =   0.8075589522942023E-01;
    weight[23-1] =   0.7375597473770516E-01;
    weight[24-1] =   0.6597422988218052E-01;
    weight[25-1] =   0.5749315621761905E-01;
    weight[26-1] =   0.4840267283059405E-01;
    weight[27-1] =   0.3879919256962704E-01;
    weight[28-1] =   0.2878470788332330E-01;
    weight[29-1] =   0.1846646831109099E-01;
    weight[30-1] =   0.7968192496166648E-02;
  }
  else if ( order == 31 )
  {
    xtab[ 1-1] =  -0.99708748181947707454263838179654;
    xtab[ 2-1] =  -0.98468590966515248400211329970113;
    xtab[ 3-1] =  -0.96250392509294966178905249675943;
    xtab[ 4-1] =  -0.93075699789664816495694576311725;
    xtab[ 5-1] =  -0.88976002994827104337419200908023;
    xtab[ 6-1] =  -0.83992032014626734008690453594388;
    xtab[ 7-1] =  -0.78173314841662494040636002019484;
    xtab[ 8-1] =  -0.71577678458685328390597086536649;
    xtab[ 9-1] =  -0.64270672292426034618441820323250;
    xtab[10-1] =  -0.56324916140714926272094492359516;
    xtab[11-1] =  -0.47819378204490248044059403935649;
    xtab[12-1] =  -0.38838590160823294306135146128752;
    xtab[13-1] =  -0.29471806998170161661790389767170;
    xtab[14-1] =  -0.19812119933557062877241299603283;
    xtab[15-1] =  -0.99555312152341520325174790118941E-01;
    xtab[16-1] =   0.00000000000000000000000000000000;
    xtab[17-1] =   0.99555312152341520325174790118941E-01;
    xtab[18-1] =   0.19812119933557062877241299603283;
    xtab[19-1] =   0.29471806998170161661790389767170;
    xtab[20-1] =   0.38838590160823294306135146128752;
    xtab[21-1] =   0.47819378204490248044059403935649;
    xtab[22-1] =   0.56324916140714926272094492359516;
    xtab[23-1] =   0.64270672292426034618441820323250;
    xtab[24-1] =   0.71577678458685328390597086536649;
    xtab[25-1] =   0.78173314841662494040636002019484;
    xtab[26-1] =   0.83992032014626734008690453594388;
    xtab[27-1] =   0.88976002994827104337419200908023;
    xtab[28-1] =   0.93075699789664816495694576311725;
    xtab[29-1] =   0.96250392509294966178905249675943;
    xtab[30-1] =   0.98468590966515248400211329970113;
    xtab[31-1] =   0.99708748181947707454263838179654;

    weight[ 1-1] =   0.74708315792487746093913218970494E-02;
    weight[ 2-1] =   0.17318620790310582463552990782414E-01;
    weight[ 3-1] =   0.27009019184979421800608642617676E-01;
    weight[ 4-1] =   0.36432273912385464024392008749009E-01;
    weight[ 5-1] =   0.45493707527201102902315857856518E-01;
    weight[ 6-1] =   0.54103082424916853711666259085477E-01;
    weight[ 7-1] =   0.62174786561028426910343543686657E-01;
    weight[ 8-1] =   0.69628583235410366167756126255124E-01;
    weight[ 9-1] =   0.76390386598776616426357674901331E-01;
    weight[10-1] =   0.82392991761589263903823367431962E-01;
    weight[11-1] =   0.87576740608477876126198069695333E-01;
    weight[12-1] =   0.91890113893641478215362871607150E-01;
    weight[13-1] =   0.95290242912319512807204197487597E-01;
    weight[14-1] =   0.97743335386328725093474010978997E-01;
    weight[15-1] =   0.99225011226672307874875514428615E-01;
    weight[16-1] =   0.99720544793426451427533833734349E-01;
    weight[17-1] =   0.99225011226672307874875514428615E-01;
    weight[18-1] =   0.97743335386328725093474010978997E-01;
    weight[19-1] =   0.95290242912319512807204197487597E-01;
    weight[20-1] =   0.91890113893641478215362871607150E-01;
    weight[21-1] =   0.87576740608477876126198069695333E-01;
    weight[22-1] =   0.82392991761589263903823367431962E-01;
    weight[23-1] =   0.76390386598776616426357674901331E-01;
    weight[24-1] =   0.69628583235410366167756126255124E-01;
    weight[25-1] =   0.62174786561028426910343543686657E-01;
    weight[26-1] =   0.54103082424916853711666259085477E-01;
    weight[27-1] =   0.45493707527201102902315857856518E-01;
    weight[28-1] =   0.36432273912385464024392008749009E-01;
    weight[29-1] =   0.27009019184979421800608642617676E-01;
    weight[30-1] =   0.17318620790310582463552990782414E-01;
    weight[31-1] =   0.74708315792487746093913218970494E-02;
  }
  else if ( order == 32 )
  {
    xtab[1-1] =  - 0.997263861849481563544981128665;
    xtab[2-1] =  - 0.985611511545268335400175044631;
    xtab[3-1] =  - 0.964762255587506430773811928118;
    xtab[4-1] =  - 0.934906075937739689170919134835;
    xtab[5-1] =  - 0.896321155766052123965307243719;
    xtab[6-1] =  - 0.849367613732569970133693004968;
    xtab[7-1] =  - 0.794483795967942406963097298970;
    xtab[8-1] =  - 0.732182118740289680387426665091;
    xtab[9-1] =  - 0.663044266930215200975115168663;
    xtab[10-1] = - 0.587715757240762329040745476402;
    xtab[11-1] = - 0.506899908932229390023747474378;
    xtab[12-1] = - 0.421351276130635345364119436172;
    xtab[13-1] = - 0.331868602282127649779916805730;
    xtab[14-1] = - 0.239287362252137074544603209166;
    xtab[15-1] = - 0.144471961582796493485186373599;
    xtab[16-1] = - 0.483076656877383162348125704405E-01;
    xtab[17-1] =   0.483076656877383162348125704405E-01;
    xtab[18-1] =   0.144471961582796493485186373599;
    xtab[19-1] =   0.239287362252137074544603209166;
    xtab[20-1] =   0.331868602282127649779916805730;
    xtab[21-1] =   0.421351276130635345364119436172;
    xtab[22-1] =   0.506899908932229390023747474378;
    xtab[23-1] =   0.587715757240762329040745476402;
    xtab[24-1] =   0.663044266930215200975115168663;
    xtab[25-1] =   0.732182118740289680387426665091;
    xtab[26-1] =   0.794483795967942406963097298970;
    xtab[27-1] =   0.849367613732569970133693004968;
    xtab[28-1] =   0.896321155766052123965307243719;
    xtab[29-1] =   0.934906075937739689170919134835;
    xtab[30-1] =   0.964762255587506430773811928118;
    xtab[31-1] =   0.985611511545268335400175044631;
    xtab[32-1] =   0.997263861849481563544981128665;

    weight[1-1] =  0.701861000947009660040706373885E-02;
    weight[2-1] =  0.162743947309056706051705622064E-01;
    weight[3-1] =  0.253920653092620594557525897892E-01;
    weight[4-1] =  0.342738629130214331026877322524E-01;
    weight[5-1] =  0.428358980222266806568786466061E-01;
    weight[6-1] =  0.509980592623761761961632446895E-01;
    weight[7-1] =  0.586840934785355471452836373002E-01;
    weight[8-1] =  0.658222227763618468376500637069E-01;
    weight[9-1] =  0.723457941088485062253993564785E-01;
    weight[10-1] = 0.781938957870703064717409188283E-01;
    weight[11-1] = 0.833119242269467552221990746043E-01;
    weight[12-1] = 0.876520930044038111427714627518E-01;
    weight[13-1] = 0.911738786957638847128685771116E-01;
    weight[14-1] = 0.938443990808045656391802376681E-01;
    weight[15-1] = 0.956387200792748594190820022041E-01;
    weight[16-1] = 0.965400885147278005667648300636E-01;
    weight[17-1] = 0.965400885147278005667648300636E-01;
    weight[18-1] = 0.956387200792748594190820022041E-01;
    weight[19-1] = 0.938443990808045656391802376681E-01;
    weight[20-1] = 0.911738786957638847128685771116E-01;
    weight[21-1] = 0.876520930044038111427714627518E-01;
    weight[22-1] = 0.833119242269467552221990746043E-01;
    weight[23-1] = 0.781938957870703064717409188283E-01;
    weight[24-1] = 0.723457941088485062253993564785E-01;
    weight[25-1] = 0.658222227763618468376500637069E-01;
    weight[26-1] = 0.586840934785355471452836373002E-01;
    weight[27-1] = 0.509980592623761761961632446895E-01;
    weight[28-1] = 0.428358980222266806568786466061E-01;
    weight[29-1] = 0.342738629130214331026877322524E-01;
    weight[30-1] = 0.253920653092620594557525897892E-01;
    weight[31-1] = 0.162743947309056706051705622064E-01;
    weight[32-1] = 0.701861000947009660040706373885E-02;
  }
  else if ( order == 33 )
  {
    xtab[ 1-1] =  -0.9974246942464552;
    xtab[ 2-1] =  -0.9864557262306425;
    xtab[ 3-1] =  -0.9668229096899927;
    xtab[ 4-1] =  -0.9386943726111684;
    xtab[ 5-1] =  -0.9023167677434336;
    xtab[ 6-1] =  -0.8580096526765041;
    xtab[ 7-1] =  -0.8061623562741665;
    xtab[ 8-1] =  -0.7472304964495622;
    xtab[ 9-1] =  -0.6817319599697428;
    xtab[10-1] =  -0.6102423458363790;
    xtab[11-1] =  -0.5333899047863476;
    xtab[12-1] =  -0.4518500172724507;
    xtab[13-1] =  -0.3663392577480734;
    xtab[14-1] =  -0.2776090971524970;
    xtab[15-1] =  -0.1864392988279916;
    xtab[16-1] =  -0.09363106585473338;
    xtab[17-1] =   0.000000000000000;
    xtab[18-1] =   0.09363106585473338;
    xtab[19-1] =   0.1864392988279916;
    xtab[20-1] =   0.2776090971524970;
    xtab[21-1] =   0.3663392577480734;
    xtab[22-1] =   0.4518500172724507;
    xtab[23-1] =   0.5333899047863476;
    xtab[24-1] =   0.6102423458363790;
    xtab[25-1] =   0.6817319599697428;
    xtab[26-1] =   0.7472304964495622;
    xtab[27-1] =   0.8061623562741665;
    xtab[28-1] =   0.8580096526765041;
    xtab[29-1] =   0.9023167677434336;
    xtab[30-1] =   0.9386943726111684;
    xtab[31-1] =   0.9668229096899927;
    xtab[32-1] =   0.9864557262306425;
    xtab[33-1] =   0.9974246942464552;

    weight[ 1-1] =   0.6606227847587558E-02;
    weight[ 2-1] =   0.1532170151293465E-01;
    weight[ 3-1] =   0.2391554810174960E-01;
    weight[ 4-1] =   0.3230035863232891E-01;
    weight[ 5-1] =   0.4040154133166965E-01;
    weight[ 6-1] =   0.4814774281871162E-01;
    weight[ 7-1] =   0.5547084663166357E-01;
    weight[ 8-1] =   0.6230648253031755E-01;
    weight[ 9-1] =   0.6859457281865676E-01;
    weight[10-1] =   0.7427985484395420E-01;
    weight[11-1] =   0.7931236479488685E-01;
    weight[12-1] =   0.8364787606703869E-01;
    weight[13-1] =   0.8724828761884425E-01;
    weight[14-1] =   0.9008195866063859E-01;
    weight[15-1] =   0.9212398664331678E-01;
    weight[16-1] =   0.9335642606559612E-01;
    weight[17-1] =   0.9376844616020999E-01;
    weight[18-1] =   0.9335642606559612E-01;
    weight[19-1] =   0.9212398664331678E-01;
    weight[20-1] =   0.9008195866063859E-01;
    weight[21-1] =   0.8724828761884425E-01;
    weight[22-1] =   0.8364787606703869E-01;
    weight[23-1] =   0.7931236479488685E-01;
    weight[24-1] =   0.7427985484395420E-01;
    weight[25-1] =   0.6859457281865676E-01;
    weight[26-1] =   0.6230648253031755E-01;
    weight[27-1] =   0.5547084663166357E-01;
    weight[28-1] =   0.4814774281871162E-01;
    weight[29-1] =   0.4040154133166965E-01;
    weight[30-1] =   0.3230035863232891E-01;
    weight[31-1] =   0.2391554810174960E-01;
    weight[32-1] =   0.1532170151293465E-01;
    weight[33-1] =   0.6606227847587558E-02;
  }
  else if ( order == 64 )
  {
    xtab[1-1] =  - 0.999305041735772139456905624346;
    xtab[2-1] =  - 0.996340116771955279346924500676;
    xtab[3-1] =  - 0.991013371476744320739382383443;
    xtab[4-1] =  - 0.983336253884625956931299302157;
    xtab[5-1] =  - 0.973326827789910963741853507352;
    xtab[6-1] =  - 0.961008799652053718918614121897;
    xtab[7-1] =  - 0.946411374858402816062481491347;
    xtab[8-1] =  - 0.929569172131939575821490154559;
    xtab[9-1] =  - 0.910522137078502805756380668008;
    xtab[10-1] = - 0.889315445995114105853404038273;
    xtab[11-1] = - 0.865999398154092819760783385070;
    xtab[12-1] = - 0.840629296252580362751691544696;
    xtab[13-1] = - 0.813265315122797559741923338086;
    xtab[14-1] = - 0.783972358943341407610220525214;
    xtab[15-1] = - 0.752819907260531896611863774886;
    xtab[16-1] = - 0.719881850171610826848940217832;
    xtab[17-1] = - 0.685236313054233242563558371031;
    xtab[18-1] = - 0.648965471254657339857761231993;
    xtab[19-1] = - 0.611155355172393250248852971019;
    xtab[20-1] = - 0.571895646202634034283878116659;
    xtab[21-1] = - 0.531279464019894545658013903544;
    xtab[22-1] = - 0.489403145707052957478526307022;
    xtab[23-1] = - 0.446366017253464087984947714759;
    xtab[24-1] = - 0.402270157963991603695766771260;
    xtab[25-1] = - 0.357220158337668115950442615046;
    xtab[26-1] = - 0.311322871990210956157512698560;
    xtab[27-1] = - 0.264687162208767416373964172510;
    xtab[28-1] = - 0.217423643740007084149648748989;
    xtab[29-1] = - 0.169644420423992818037313629748;
    xtab[30-1] = - 0.121462819296120554470376463492;
    xtab[31-1] = - 0.729931217877990394495429419403E-01;
    xtab[32-1] = - 0.243502926634244325089558428537E-01;
    xtab[33-1] =   0.243502926634244325089558428537E-01;
    xtab[34-1] =   0.729931217877990394495429419403E-01;
    xtab[35-1] =   0.121462819296120554470376463492;
    xtab[36-1] =   0.169644420423992818037313629748;
    xtab[37-1] =   0.217423643740007084149648748989;
    xtab[38-1] =   0.264687162208767416373964172510;
    xtab[39-1] =   0.311322871990210956157512698560;
    xtab[40-1] =   0.357220158337668115950442615046;
    xtab[41-1] =   0.402270157963991603695766771260;
    xtab[42-1] =   0.446366017253464087984947714759;
    xtab[43-1] =   0.489403145707052957478526307022;
    xtab[44-1] =   0.531279464019894545658013903544;
    xtab[45-1] =   0.571895646202634034283878116659;
    xtab[46-1] =   0.611155355172393250248852971019;
    xtab[47-1] =   0.648965471254657339857761231993;
    xtab[48-1] =   0.685236313054233242563558371031;
    xtab[49-1] =   0.719881850171610826848940217832;
    xtab[50-1] =   0.752819907260531896611863774886;
    xtab[51-1] =   0.783972358943341407610220525214;
    xtab[52-1] =   0.813265315122797559741923338086;
    xtab[53-1] =   0.840629296252580362751691544696;
    xtab[54-1] =   0.865999398154092819760783385070;
    xtab[55-1] =   0.889315445995114105853404038273;
    xtab[56-1] =   0.910522137078502805756380668008;
    xtab[57-1] =   0.929569172131939575821490154559;
    xtab[58-1] =   0.946411374858402816062481491347;
    xtab[59-1] =   0.961008799652053718918614121897;
    xtab[60-1] =   0.973326827789910963741853507352;
    xtab[61-1] =   0.983336253884625956931299302157;
    xtab[62-1] =   0.991013371476744320739382383443;
    xtab[63-1] =   0.996340116771955279346924500676;
    xtab[64-1] =   0.999305041735772139456905624346;

    weight[1-1] =  0.178328072169643294729607914497E-02;
    weight[2-1] =  0.414703326056246763528753572855E-02;
    weight[3-1] =  0.650445796897836285611736039998E-02;
    weight[4-1] =  0.884675982636394772303091465973E-02;
    weight[5-1] =  0.111681394601311288185904930192E-01;
    weight[6-1] =  0.134630478967186425980607666860E-01;
    weight[7-1] =  0.157260304760247193219659952975E-01;
    weight[8-1] =  0.179517157756973430850453020011E-01;
    weight[9-1] =  0.201348231535302093723403167285E-01;
    weight[10-1] = 0.222701738083832541592983303842E-01;
    weight[11-1] = 0.243527025687108733381775504091E-01;
    weight[12-1] = 0.263774697150546586716917926252E-01;
    weight[13-1] = 0.283396726142594832275113052002E-01;
    weight[14-1] = 0.302346570724024788679740598195E-01;
    weight[15-1] = 0.320579283548515535854675043479E-01;
    weight[16-1] = 0.338051618371416093915654821107E-01;
    weight[17-1] = 0.354722132568823838106931467152E-01;
    weight[18-1] = 0.370551285402400460404151018096E-01;
    weight[19-1] = 0.385501531786156291289624969468E-01;
    weight[20-1] = 0.399537411327203413866569261283E-01;
    weight[21-1] = 0.412625632426235286101562974736E-01;
    weight[22-1] = 0.424735151236535890073397679088E-01;
    weight[23-1] = 0.435837245293234533768278609737E-01;
    weight[24-1] = 0.445905581637565630601347100309E-01;
    weight[25-1] = 0.454916279274181444797709969713E-01;
    weight[26-1] = 0.462847965813144172959532492323E-01;
    weight[27-1] = 0.469681828162100173253262857546E-01;
    weight[28-1] = 0.475401657148303086622822069442E-01;
    weight[29-1] = 0.479993885964583077281261798713E-01;
    weight[30-1] = 0.483447622348029571697695271580E-01;
    weight[31-1] = 0.485754674415034269347990667840E-01;
    weight[32-1] = 0.486909570091397203833653907347E-01;
    weight[33-1] = 0.486909570091397203833653907347E-01;
    weight[34-1] = 0.485754674415034269347990667840E-01;
    weight[35-1] = 0.483447622348029571697695271580E-01;
    weight[36-1] = 0.479993885964583077281261798713E-01;
    weight[37-1] = 0.475401657148303086622822069442E-01;
    weight[38-1] = 0.469681828162100173253262857546E-01;
    weight[39-1] = 0.462847965813144172959532492323E-01;
    weight[40-1] = 0.454916279274181444797709969713E-01;
    weight[41-1] = 0.445905581637565630601347100309E-01;
    weight[42-1] = 0.435837245293234533768278609737E-01;
    weight[43-1] = 0.424735151236535890073397679088E-01;
    weight[44-1] = 0.412625632426235286101562974736E-01;
    weight[45-1] = 0.399537411327203413866569261283E-01;
    weight[46-1] = 0.385501531786156291289624969468E-01;
    weight[47-1] = 0.370551285402400460404151018096E-01;
    weight[48-1] = 0.354722132568823838106931467152E-01;
    weight[49-1] = 0.338051618371416093915654821107E-01;
    weight[50-1] = 0.320579283548515535854675043479E-01;
    weight[51-1] = 0.302346570724024788679740598195E-01;
    weight[52-1] = 0.283396726142594832275113052002E-01;
    weight[53-1] = 0.263774697150546586716917926252E-01;
    weight[54-1] = 0.243527025687108733381775504091E-01;
    weight[55-1] = 0.222701738083832541592983303842E-01;
    weight[56-1] = 0.201348231535302093723403167285E-01;
    weight[57-1] = 0.179517157756973430850453020011E-01;
    weight[58-1] = 0.157260304760247193219659952975E-01;
    weight[59-1] = 0.134630478967186425980607666860E-01;
    weight[60-1] = 0.111681394601311288185904930192E-01;
    weight[61-1] = 0.884675982636394772303091465973E-02;
    weight[62-1] = 0.650445796897836285611736039998E-02;
    weight[63-1] = 0.414703326056246763528753572855E-02;
    weight[64-1] = 0.178328072169643294729607914497E-02;
  }
  else if ( order == 65 )
  {
    xtab[ 1-1] =  -0.9993260970754129;
    xtab[ 2-1] =  -0.9964509480618492;
    xtab[ 3-1] =  -0.9912852761768016;
    xtab[ 4-1] =  -0.9838398121870350;
    xtab[ 5-1] =  -0.9741315398335512;
    xtab[ 6-1] =  -0.9621827547180553;
    xtab[ 7-1] =  -0.9480209281684076;
    xtab[ 8-1] =  -0.9316786282287494;
    xtab[ 9-1] =  -0.9131934405428462;
    xtab[10-1] =  -0.8926078805047389;
    xtab[11-1] =  -0.8699692949264071;
    xtab[12-1] =  -0.8453297528999303;
    xtab[13-1] =  -0.8187459259226514;
    xtab[14-1] =  -0.7902789574921218;
    xtab[15-1] =  -0.7599943224419998;
    xtab[16-1] =  -0.7279616763294247;
    xtab[17-1] =  -0.6942546952139916;
    xtab[18-1] =  -0.6589509061936252;
    xtab[19-1] =  -0.6221315090854003;
    xtab[20-1] =  -0.5838811896604873;
    xtab[21-1] =  -0.5442879248622271;
    xtab[22-1] =  -0.5034427804550069;
    xtab[23-1] =  -0.4614397015691450;
    xtab[24-1] =  -0.4183752966234090;
    xtab[25-1] =  -0.3743486151220660;
    xtab[26-1] =  -0.3294609198374864;
    xtab[27-1] =  -0.2838154539022487;
    xtab[28-1] =  -0.2375172033464168;
    xtab[29-1] =  -0.1906726556261428;
    xtab[30-1] =  -0.1433895546989752;
    xtab[31-1] =  -0.9577665320919751E-01;
    xtab[32-1] =  -0.4794346235317186E-01;
    xtab[33-1] =    0.000000000000000;
    xtab[34-1] =   0.4794346235317186E-01;
    xtab[35-1] =   0.9577665320919751E-01;
    xtab[36-1] =   0.1433895546989752;
    xtab[37-1] =   0.1906726556261428;
    xtab[38-1] =   0.2375172033464168;
    xtab[39-1] =   0.2838154539022487;
    xtab[40-1] =   0.3294609198374864;
    xtab[41-1] =   0.3743486151220660;
    xtab[42-1] =   0.4183752966234090;
    xtab[43-1] =   0.4614397015691450;
    xtab[44-1] =   0.5034427804550069;
    xtab[45-1] =   0.5442879248622271;
    xtab[46-1] =   0.5838811896604873;
    xtab[47-1] =   0.6221315090854003;
    xtab[48-1] =   0.6589509061936252;
    xtab[49-1] =   0.6942546952139916;
    xtab[50-1] =   0.7279616763294247;
    xtab[51-1] =   0.7599943224419998;
    xtab[52-1] =   0.7902789574921218;
    xtab[53-1] =   0.8187459259226514;
    xtab[54-1] =   0.8453297528999303;
    xtab[55-1] =   0.8699692949264071;
    xtab[56-1] =   0.8926078805047389;
    xtab[57-1] =   0.9131934405428462;
    xtab[58-1] =   0.9316786282287494;
    xtab[59-1] =   0.9480209281684076;
    xtab[60-1] =   0.9621827547180553;
    xtab[61-1] =   0.9741315398335512;
    xtab[62-1] =   0.9838398121870350;
    xtab[63-1] =   0.9912852761768016;
    xtab[64-1] =   0.9964509480618492;
    xtab[65-1] =   0.9993260970754129;

    weight[ 1-1] =   0.1729258251300218E-02;
    weight[ 2-1] =   0.4021524172003703E-02;
    weight[ 3-1] =   0.6307942578971821E-02;
    weight[ 4-1] =   0.8580148266881443E-02;
    weight[ 5-1] =   0.1083267878959798E-01;
    weight[ 6-1] =   0.1306031163999490E-01;
    weight[ 7-1] =   0.1525791214644825E-01;
    weight[ 8-1] =   0.1742042199767025E-01;
    weight[ 9-1] =   0.1954286583675005E-01;
    weight[10-1] =   0.2162036128493408E-01;
    weight[11-1] =   0.2364812969128723E-01;
    weight[12-1] =   0.2562150693803776E-01;
    weight[13-1] =   0.2753595408845034E-01;
    weight[14-1] =   0.2938706778931066E-01;
    weight[15-1] =   0.3117059038018911E-01;
    weight[16-1] =   0.3288241967636860E-01;
    weight[17-1] =   0.3451861839854901E-01;
    weight[18-1] =   0.3607542322556527E-01;
    weight[19-1] =   0.3754925344825770E-01;
    weight[20-1] =   0.3893671920405121E-01;
    weight[21-1] =   0.4023462927300549E-01;
    weight[22-1] =   0.4143999841724028E-01;
    weight[23-1] =   0.4255005424675579E-01;
    weight[24-1] =   0.4356224359580051E-01;
    weight[25-1] =   0.4447423839508296E-01;
    weight[26-1] =   0.4528394102630023E-01;
    weight[27-1] =   0.4598948914665173E-01;
    weight[28-1] =   0.4658925997223349E-01;
    weight[29-1] =   0.4708187401045461E-01;
    weight[30-1] =   0.4746619823288551E-01;
    weight[31-1] =   0.4774134868124067E-01;
    weight[32-1] =   0.4790669250049590E-01;
    weight[33-1] =   0.4796184939446662E-01;
    weight[34-1] =   0.4790669250049590E-01;
    weight[35-1] =   0.4774134868124067E-01;
    weight[36-1] =   0.4746619823288551E-01;
    weight[37-1] =   0.4708187401045461E-01;
    weight[38-1] =   0.4658925997223349E-01;
    weight[39-1] =   0.4598948914665173E-01;
    weight[40-1] =   0.4528394102630023E-01;
    weight[41-1] =   0.4447423839508296E-01;
    weight[42-1] =   0.4356224359580051E-01;
    weight[43-1] =   0.4255005424675579E-01;
    weight[44-1] =   0.4143999841724028E-01;
    weight[45-1] =   0.4023462927300549E-01;
    weight[46-1] =   0.3893671920405121E-01;
    weight[47-1] =   0.3754925344825770E-01;
    weight[48-1] =   0.3607542322556527E-01;
    weight[49-1] =   0.3451861839854901E-01;
    weight[50-1] =   0.3288241967636860E-01;
    weight[51-1] =   0.3117059038018911E-01;
    weight[52-1] =   0.2938706778931066E-01;
    weight[53-1] =   0.2753595408845034E-01;
    weight[54-1] =   0.2562150693803776E-01;
    weight[55-1] =   0.2364812969128723E-01;
    weight[56-1] =   0.2162036128493408E-01;
    weight[57-1] =   0.1954286583675005E-01;
    weight[58-1] =   0.1742042199767025E-01;
    weight[59-1] =   0.1525791214644825E-01;
    weight[60-1] =   0.1306031163999490E-01;
    weight[61-1] =   0.1083267878959798E-01;
    weight[62-1] =   0.8580148266881443E-02;
    weight[63-1] =   0.6307942578971821E-02;
    weight[64-1] =   0.4021524172003703E-02;
    weight[65-1] =   0.1729258251300218E-02;
  }
  else if ( order == 127 )
  {
    xtab[  1-1] =  -0.99982213041530614629963254927125E+00;
    xtab[  2-1] =  -0.99906293435531189513828920479421E+00;
    xtab[  3-1] =  -0.99769756618980462107441703193392E+00;
    xtab[  4-1] =  -0.99572655135202722663543337085008E+00;
    xtab[  5-1] =  -0.99315104925451714736113079489080E+00;
    xtab[  6-1] =  -0.98997261459148415760778669967548E+00;
    xtab[  7-1] =  -0.98619317401693166671043833175407E+00;
    xtab[  8-1] =  -0.98181502080381411003346312451200E+00;
    xtab[  9-1] =  -0.97684081234307032681744391886221E+00;
    xtab[ 10-1] =  -0.97127356816152919228894689830512E+00;
    xtab[ 11-1] =  -0.96511666794529212109082507703391E+00;
    xtab[ 12-1] =  -0.95837384942523877114910286998060E+00;
    xtab[ 13-1] =  -0.95104920607788031054790764659636E+00;
    xtab[ 14-1] =  -0.94314718462481482734544963026201E+00;
    xtab[ 15-1] =  -0.93467258232473796857363487794906E+00;
    xtab[ 16-1] =  -0.92563054405623384912746466814259E+00;
    xtab[ 17-1] =  -0.91602655919146580931308861741716E+00;
    xtab[ 18-1] =  -0.90586645826182138280246131760282E+00;
    xtab[ 19-1] =  -0.89515640941708370896904382642451E+00;
    xtab[ 20-1] =  -0.88390291468002656994525794802849E+00;
    xtab[ 21-1] =  -0.87211280599856071141963753428864E+00;
    xtab[ 22-1] =  -0.85979324109774080981203134414483E+00;
    xtab[ 23-1] =  -0.84695169913409759845333931085437E+00;
    xtab[ 24-1] =  -0.83359597615489951437955716480123E+00;
    xtab[ 25-1] =  -0.81973418036507867415511910167470E+00;
    xtab[ 26-1] =  -0.80537472720468021466656079404644E+00;
    xtab[ 27-1] =  -0.79052633423981379994544995252740E+00;
    xtab[ 28-1] =  -0.77519801587020238244496276354566E+00;
    xtab[ 29-1] =  -0.75939907785653667155666366659810E+00;
    xtab[ 30-1] =  -0.74313911167095451292056688997595E+00;
    xtab[ 31-1] =  -0.72642798867407268553569290153270E+00;
    xtab[ 32-1] =  -0.70927585412210456099944463906757E+00;
    xtab[ 33-1] =  -0.69169312100770067015644143286666E+00;
    xtab[ 34-1] =  -0.67369046373825048534668253831602E+00;
    xtab[ 35-1] =  -0.65527881165548263027676505156852E+00;
    xtab[ 36-1] =  -0.63646934240029724134760815684175E+00;
    xtab[ 37-1] =  -0.61727347512685828385763916340822E+00;
    xtab[ 38-1] =  -0.59770286357006522938441201887478E+00;
    xtab[ 39-1] =  -0.57776938897061258000325165713764E+00;
    xtab[ 40-1] =  -0.55748515286193223292186190687872E+00;
    xtab[ 41-1] =  -0.53686246972339756745816636353452E+00;
    xtab[ 42-1] =  -0.51591385950424935727727729906662E+00;
    xtab[ 43-1] =  -0.49465204002278211739494017368636E+00;
    xtab[ 44-1] =  -0.47308991924540524164509989939699E+00;
    xtab[ 45-1] =  -0.45124058745026622733189858020729E+00;
    xtab[ 46-1] =  -0.42911730928019337626254405355418E+00;
    xtab[ 47-1] =  -0.40673351568978256340867288124339E+00;
    xtab[ 48-1] =  -0.38410279579151693577907781452239E+00;
    xtab[ 49-1] =  -0.36123888860586970607092484346723E+00;
    xtab[ 50-1] =  -0.33815567472039850137600027657095E+00;
    xtab[ 51-1] =  -0.31486716786289498148601475374890E+00;
    xtab[ 52-1] =  -0.29138750639370562079451875284568E+00;
    xtab[ 53-1] =  -0.26773094472238862088834352027938E+00;
    xtab[ 54-1] =  -0.24391184465391785797071324453138E+00;
    xtab[ 55-1] =  -0.21994466666968754245452337866940E+00;
    xtab[ 56-1] =  -0.19584396114861085150428162519610E+00;
    xtab[ 57-1] =  -0.17162435953364216500834492248954E+00;
    xtab[ 58-1] =  -0.14730056544908566938932929319807E+00;
    xtab[ 59-1] =  -0.12288734577408297172603365288567E+00;
    xtab[ 60-1] =  -0.98399521677698970751091751509101E-01;
    xtab[ 61-1] =  -0.73851959621048545273440409360569E-01;
    xtab[ 62-1] =  -0.49259562331926630315379321821927E-01;
    xtab[ 63-1] =  -0.24637259757420944614897071846088E-01;
    xtab[ 64-1] =   0.00000000000000000000000000000000E+00;
    xtab[ 65-1] =   0.24637259757420944614897071846088E-01;
    xtab[ 66-1] =   0.49259562331926630315379321821927E-01;
    xtab[ 67-1] =   0.73851959621048545273440409360569E-01;
    xtab[ 68-1] =   0.98399521677698970751091751509101E-01;
    xtab[ 69-1] =   0.12288734577408297172603365288567E+00;
    xtab[ 70-1] =   0.14730056544908566938932929319807E+00;
    xtab[ 71-1] =   0.17162435953364216500834492248954E+00;
    xtab[ 72-1] =   0.19584396114861085150428162519610E+00;
    xtab[ 73-1] =   0.21994466666968754245452337866940E+00;
    xtab[ 74-1] =   0.24391184465391785797071324453138E+00;
    xtab[ 75-1] =   0.26773094472238862088834352027938E+00;
    xtab[ 76-1] =   0.29138750639370562079451875284568E+00;
    xtab[ 77-1] =   0.31486716786289498148601475374890E+00;
    xtab[ 78-1] =   0.33815567472039850137600027657095E+00;
    xtab[ 79-1] =   0.36123888860586970607092484346723E+00;
    xtab[ 80-1] =   0.38410279579151693577907781452239E+00;
    xtab[ 81-1] =   0.40673351568978256340867288124339E+00;
    xtab[ 82-1] =   0.42911730928019337626254405355418E+00;
    xtab[ 83-1] =   0.45124058745026622733189858020729E+00;
    xtab[ 84-1] =   0.47308991924540524164509989939699E+00;
    xtab[ 85-1] =   0.49465204002278211739494017368636E+00;
    xtab[ 86-1] =   0.51591385950424935727727729906662E+00;
    xtab[ 87-1] =   0.53686246972339756745816636353452E+00;
    xtab[ 88-1] =   0.55748515286193223292186190687872E+00;
    xtab[ 89-1] =   0.57776938897061258000325165713764E+00;
    xtab[ 90-1] =   0.59770286357006522938441201887478E+00;
    xtab[ 91-1] =   0.61727347512685828385763916340822E+00;
    xtab[ 92-1] =   0.63646934240029724134760815684175E+00;
    xtab[ 93-1] =   0.65527881165548263027676505156852E+00;
    xtab[ 94-1] =   0.67369046373825048534668253831602E+00;
    xtab[ 95-1] =   0.69169312100770067015644143286666E+00;
    xtab[ 96-1] =   0.70927585412210456099944463906757E+00;
    xtab[ 97-1] =   0.72642798867407268553569290153270E+00;
    xtab[ 98-1] =   0.74313911167095451292056688997595E+00;
    xtab[ 99-1] =   0.75939907785653667155666366659810E+00;
    xtab[100-1] =   0.77519801587020238244496276354566E+00;
    xtab[101-1] =   0.79052633423981379994544995252740E+00;
    xtab[102-1] =   0.80537472720468021466656079404644E+00;
    xtab[103-1] =   0.81973418036507867415511910167470E+00;
    xtab[104-1] =   0.83359597615489951437955716480123E+00;
    xtab[105-1] =   0.84695169913409759845333931085437E+00;
    xtab[106-1] =   0.85979324109774080981203134414483E+00;
    xtab[107-1] =   0.87211280599856071141963753428864E+00;
    xtab[108-1] =   0.88390291468002656994525794802849E+00;
    xtab[109-1] =   0.89515640941708370896904382642451E+00;
    xtab[110-1] =   0.90586645826182138280246131760282E+00;
    xtab[111-1] =   0.91602655919146580931308861741716E+00;
    xtab[112-1] =   0.92563054405623384912746466814259E+00;
    xtab[113-1] =   0.93467258232473796857363487794906E+00;
    xtab[114-1] =   0.94314718462481482734544963026201E+00;
    xtab[115-1] =   0.95104920607788031054790764659636E+00;
    xtab[116-1] =   0.95837384942523877114910286998060E+00;
    xtab[117-1] =   0.96511666794529212109082507703391E+00;
    xtab[118-1] =   0.97127356816152919228894689830512E+00;
    xtab[119-1] =   0.97684081234307032681744391886221E+00;
    xtab[120-1] =   0.98181502080381411003346312451200E+00;
    xtab[121-1] =   0.98619317401693166671043833175407E+00;
    xtab[122-1] =   0.98997261459148415760778669967548E+00;
    xtab[123-1] =   0.99315104925451714736113079489080E+00;
    xtab[124-1] =   0.99572655135202722663543337085008E+00;
    xtab[125-1] =   0.99769756618980462107441703193392E+00;
    xtab[126-1] =   0.99906293435531189513828920479421E+00;
    xtab[127-1] =   0.99982213041530614629963254927125E+00;

    weight[  1-1] =   0.45645726109586654495731936146574E-03;
    weight[  2-1] =   0.10622766869538486959954760554099E-02;
    weight[  3-1] =   0.16683488125171936761028811985672E-02;
    weight[  4-1] =   0.22734860707492547802810838362671E-02;
    weight[  5-1] =   0.28772587656289004082883197417581E-02;
    weight[  6-1] =   0.34792893810051465908910894094105E-02;
    weight[  7-1] =   0.40792095178254605327114733456293E-02;
    weight[  8-1] =   0.46766539777779034772638165662478E-02;
    weight[  9-1] =   0.52712596565634400891303815906251E-02;
    weight[ 10-1] =   0.58626653903523901033648343751367E-02;
    weight[ 11-1] =   0.64505120486899171845442463868748E-02;
    weight[ 12-1] =   0.70344427036681608755685893032552E-02;
    weight[ 13-1] =   0.76141028256526859356393930849227E-02;
    weight[ 14-1] =   0.81891404887415730817235884718726E-02;
    weight[ 15-1] =   0.87592065795403145773316804234385E-02;
    weight[ 16-1] =   0.93239550065309714787536985834029E-02;
    weight[ 17-1] =   0.98830429087554914716648010899606E-02;
    weight[ 18-1] =   0.10436130863141005225673171997668E-01;
    weight[ 19-1] =   0.10982883090068975788799657376065E-01;
    weight[ 20-1] =   0.11522967656921087154811609734510E-01;
    weight[ 21-1] =   0.12056056679400848183529562144697E-01;
    weight[ 22-1] =   0.12581826520465013101514365424172E-01;
    weight[ 23-1] =   0.13099957986718627426172681912499E-01;
    weight[ 24-1] =   0.13610136522139249906034237533759E-01;
    weight[ 25-1] =   0.14112052399003395774044161633613E-01;
    weight[ 26-1] =   0.14605400905893418351737288078952E-01;
    weight[ 27-1] =   0.15089882532666922992635733981431E-01;
    weight[ 28-1] =   0.15565203152273955098532590262975E-01;
    weight[ 29-1] =   0.16031074199309941802254151842763E-01;
    weight[ 30-1] =   0.16487212845194879399346060358146E-01;
    weight[ 31-1] =   0.16933342169871654545878815295200E-01;
    weight[ 32-1] =   0.17369191329918731922164721250350E-01;
    weight[ 33-1] =   0.17794495722974774231027912900351E-01;
    weight[ 34-1] =   0.18208997148375106468721469154479E-01;
    weight[ 35-1] =   0.18612443963902310429440419898958E-01;
    weight[ 36-1] =   0.19004591238555646611148901044533E-01;
    weight[ 37-1] =   0.19385200901246454628112623489471E-01;
    weight[ 38-1] =   0.19754041885329183081815217323169E-01;
    weight[ 39-1] =   0.20110890268880247225644623956287E-01;
    weight[ 40-1] =   0.20455529410639508279497065713301E-01;
    weight[ 41-1] =   0.20787750081531811812652137291250E-01;
    weight[ 42-1] =   0.21107350591688713643523847921658E-01;
    weight[ 43-1] =   0.21414136912893259295449693233545E-01;
    weight[ 44-1] =   0.21707922796373466052301324695331E-01;
    weight[ 45-1] =   0.21988529885872983756478409758807E-01;
    weight[ 46-1] =   0.22255787825930280235631416460158E-01;
    weight[ 47-1] =   0.22509534365300608085694429903050E-01;
    weight[ 48-1] =   0.22749615455457959852242553240982E-01;
    weight[ 49-1] =   0.22975885344117206754377437838947E-01;
    weight[ 50-1] =   0.23188206663719640249922582981729E-01;
    weight[ 51-1] =   0.23386450514828194170722043496950E-01;
    weight[ 52-1] =   0.23570496544381716050033676844306E-01;
    weight[ 53-1] =   0.23740233018760777777714726703424E-01;
    weight[ 54-1] =   0.23895556891620665983864481754172E-01;
    weight[ 55-1] =   0.24036373866450369675132086026456E-01;
    weight[ 56-1] =   0.24162598453819584716522917710986E-01;
    weight[ 57-1] =   0.24274154023278979833195063936748E-01;
    weight[ 58-1] =   0.24370972849882214952813561907241E-01;
    weight[ 59-1] =   0.24452996155301467956140198471529E-01;
    weight[ 60-1] =   0.24520174143511508275183033290175E-01;
    weight[ 61-1] =   0.24572466031020653286354137335186E-01;
    weight[ 62-1] =   0.24609840071630254092545634003360E-01;
    weight[ 63-1] =   0.24632273575707679066033370218017E-01;
    weight[ 64-1] =   0.24639752923961094419579417477503E-01;
    weight[ 65-1] =   0.24632273575707679066033370218017E-01;
    weight[ 66-1] =   0.24609840071630254092545634003360E-01;
    weight[ 67-1] =   0.24572466031020653286354137335186E-01;
    weight[ 68-1] =   0.24520174143511508275183033290175E-01;
    weight[ 69-1] =   0.24452996155301467956140198471529E-01;
    weight[ 70-1] =   0.24370972849882214952813561907241E-01;
    weight[ 71-1] =   0.24274154023278979833195063936748E-01;
    weight[ 72-1] =   0.24162598453819584716522917710986E-01;
    weight[ 73-1] =   0.24036373866450369675132086026456E-01;
    weight[ 74-1] =   0.23895556891620665983864481754172E-01;
    weight[ 75-1] =   0.23740233018760777777714726703424E-01;
    weight[ 76-1] =   0.23570496544381716050033676844306E-01;
    weight[ 77-1] =   0.23386450514828194170722043496950E-01;
    weight[ 78-1] =   0.23188206663719640249922582981729E-01;
    weight[ 79-1] =   0.22975885344117206754377437838947E-01;
    weight[ 80-1] =   0.22749615455457959852242553240982E-01;
    weight[ 81-1] =   0.22509534365300608085694429903050E-01;
    weight[ 82-1] =   0.22255787825930280235631416460158E-01;
    weight[ 83-1] =   0.21988529885872983756478409758807E-01;
    weight[ 84-1] =   0.21707922796373466052301324695331E-01;
    weight[ 85-1] =   0.21414136912893259295449693233545E-01;
    weight[ 86-1] =   0.21107350591688713643523847921658E-01;
    weight[ 87-1] =   0.20787750081531811812652137291250E-01;
    weight[ 88-1] =   0.20455529410639508279497065713301E-01;
    weight[ 89-1] =   0.20110890268880247225644623956287E-01;
    weight[ 90-1] =   0.19754041885329183081815217323169E-01;
    weight[ 91-1] =   0.19385200901246454628112623489471E-01;
    weight[ 92-1] =   0.19004591238555646611148901044533E-01;
    weight[ 93-1] =   0.18612443963902310429440419898958E-01;
    weight[ 94-1] =   0.18208997148375106468721469154479E-01;
    weight[ 95-1] =   0.17794495722974774231027912900351E-01;
    weight[ 96-1] =   0.17369191329918731922164721250350E-01;
    weight[ 97-1] =   0.16933342169871654545878815295200E-01;
    weight[ 98-1] =   0.16487212845194879399346060358146E-01;
    weight[ 99-1] =   0.16031074199309941802254151842763E-01;
    weight[100-1] =   0.15565203152273955098532590262975E-01;
    weight[101-1] =   0.15089882532666922992635733981431E-01;
    weight[102-1] =   0.14605400905893418351737288078952E-01;
    weight[103-1] =   0.14112052399003395774044161633613E-01;
    weight[104-1] =   0.13610136522139249906034237533759E-01;
    weight[105-1] =   0.13099957986718627426172681912499E-01;
    weight[106-1] =   0.12581826520465013101514365424172E-01;
    weight[107-1] =   0.12056056679400848183529562144697E-01;
    weight[108-1] =   0.11522967656921087154811609734510E-01;
    weight[109-1] =   0.10982883090068975788799657376065E-01;
    weight[110-1] =   0.10436130863141005225673171997668E-01;
    weight[111-1] =   0.98830429087554914716648010899606E-02;
    weight[112-1] =   0.93239550065309714787536985834029E-02;
    weight[113-1] =   0.87592065795403145773316804234385E-02;
    weight[114-1] =   0.81891404887415730817235884718726E-02;
    weight[115-1] =   0.76141028256526859356393930849227E-02;
    weight[116-1] =   0.70344427036681608755685893032552E-02;
    weight[117-1] =   0.64505120486899171845442463868748E-02;
    weight[118-1] =   0.58626653903523901033648343751367E-02;
    weight[119-1] =   0.52712596565634400891303815906251E-02;
    weight[120-1] =   0.46766539777779034772638165662478E-02;
    weight[121-1] =   0.40792095178254605327114733456293E-02;
    weight[122-1] =   0.34792893810051465908910894094105E-02;
    weight[123-1] =   0.28772587656289004082883197417581E-02;
    weight[124-1] =   0.22734860707492547802810838362671E-02;
    weight[125-1] =   0.16683488125171936761028811985672E-02;
    weight[126-1] =   0.10622766869538486959954760554099E-02;
    weight[127-1] =   0.45645726109586654495731936146574E-03;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 through 33, 63, 64, 65 or 127.\n" );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void legendre_set_cos ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_COS: Gauss-Legendre rules for COS(X)*F(X) on [-PI/2,PI/2].
  
  Discussion:
  
    The integration interval is [ -PI/2, PI/2 ].
  
    The weight function is w(x-1] = cos(x).
  
    The integral to approximate:
  
      Integral ( -PI/2 <= X <= PI/2 ) COS(X) * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
  
  Discussion:
  
    The same rule can be used to approximate
  
      Integral ( 0 <= X <= PI ) SIN(X) * F(X) dX
  
    as
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) + PI/2 )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    03 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Gwynne Evans,
    Practical Numerical Integration,
    Wiley, 1993, QA299.3E93, page 310.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1, 2, 4, 8 or 16.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] = 0.0E+00;

    weight[1-1] = 2.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = - 0.68366739008990304094E+00;
    xtab[2-1] =   0.68366739008990304094E+00;

    weight[1-1] = 1.0E+00;
    weight[2-1] = 1.0E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = - 1.1906765638948557415E+00;
    xtab[2-1] = - 0.43928746686001514756E+00;
    xtab[3-1] =   0.43928746686001514756E+00;
    xtab[4-1] =   1.1906765638948557415E+00;

    weight[1-1] = 0.22407061812762016065E+00;
    weight[2-1] = 0.77592938187237983935E+00;
    weight[3-1] = 0.77592938187237983935E+00;
    weight[4-1] = 0.22407061812762016065E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = - 1.4414905401823575701E+00;
    xtab[2-1] = - 1.1537256454567275850E+00;
    xtab[3-1] = - 0.74346864787549244989E+00;
    xtab[4-1] = - 0.25649650741623123020E+00;
    xtab[5-1] =   0.25649650741623123020E+00;
    xtab[6-1] =   0.74346864787549244989E+00;
    xtab[7-1] =   1.1537256454567275850E+00;
    xtab[8-1] =   1.4414905401823575701E+00;

    weight[1-1] = 0.027535633513767011149E+00;
    weight[2-1] = 0.14420409203022750950E+00;
    weight[3-1] = 0.33626447785280459621E+00;
    weight[4-1] = 0.49199579660320088314E+00;
    weight[5-1] = 0.49199579660320088314E+00;
    weight[6-1] = 0.33626447785280459621E+00;
    weight[7-1] = 0.14420409203022750950E+00;
    weight[8-1] = 0.027535633513767011149E+00;
  }
  else if ( order == 16 )
  {
    xtab[ 1-1] = - 1.5327507132362304779E+00;
    xtab[ 2-1] = - 1.4446014873666514608E+00;
    xtab[ 3-1] = - 1.3097818904452936698E+00;
    xtab[ 4-1] = - 1.1330068786005003695E+00;
    xtab[ 5-1] = - 0.92027786206637096497E+00;
    xtab[ 6-1] = - 0.67861108097560545347E+00;
    xtab[ 7-1] = - 0.41577197673418943962E+00;
    xtab[ 8-1] = - 0.14003444424696773778E+00;
    xtab[ 9-1] =   0.14003444424696773778E+00;
    xtab[10-1] =   0.41577197673418943962E+00;
    xtab[11-1] =   0.67861108097560545347E+00;
    xtab[12-1] =   0.92027786206637096497E+00;
    xtab[13-1] =   1.1330068786005003695E+00;
    xtab[14-1] =   1.3097818904452936698E+00;
    xtab[15-1] =   1.4446014873666514608E+00;
    xtab[16-1] =   1.5327507132362304779E+00;

    weight[ 1-1] = 0.0024194677567615628193E+00;
    weight[ 2-1] = 0.014115268156854008264E+00;
    weight[ 3-1] = 0.040437893946503669410E+00;
    weight[ 4-1] = 0.083026647573217742131E+00;
    weight[ 5-1] = 0.13834195526951273359E+00;
    weight[ 6-1] = 0.19741148870253455567E+00;
    weight[ 7-1] = 0.24763632094635522403E+00;
    weight[ 8-1] = 0.27661095764826050408E+00;
    weight[ 9-1] = 0.27661095764826050408E+00;
    weight[10-1] = 0.24763632094635522403E+00;
    weight[11-1] = 0.19741148870253455567E+00;
    weight[12-1] = 0.13834195526951273359E+00;
    weight[13-1] = 0.083026647573217742131E+00;
    weight[14-1] = 0.040437893946503669410E+00;
    weight[15-1] = 0.014115268156854008264E+00;
    weight[16-1] = 0.0024194677567615628193E+00;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET_COS - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void legendre_set_cos2 ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_COS2: Gauss-Legendre rules for COS(X)*F(X) on [0,PI/2].
  
  Discussion:
  
    The integration interval is [ 0, PI/2 ].
  
    The weight function is w(x-1] = cos ( x ).
  
    The integral to approximate:
  
      Integral ( 0 <= X <= PI/2 ) COS(X) * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
  
    The same rule can be used to approximate
  
      Integral ( 0 <= X <= PI/2 ) SIN(X) * F(X) dX
  
    as
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( PI/2 - XTAb[I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    03 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Gwynne Evans,
    Practical Numerical Integration,
    Wiley, 1993, QA299.3E93, page 311.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 2, 4, 8 or 16.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 2 )
  {
    xtab[1-1] = 0.26587388056307823382E+00;
    xtab[2-1] = 1.0351526093171315182E+00;

    weight[1-1] = 0.60362553280827113087E+00;
    weight[2-1] = 0.39637446719172886913E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = 0.095669389196858636773E+00;
    xtab[2-1] = 0.45240902327067096554E+00;
    xtab[3-1] = 0.93185057672024082424E+00;
    xtab[4-1] = 1.3564439599666466230E+00;

    weight[ 1-1] = 0.23783071419515504517E+00;
    weight[ 2-1] = 0.40265695523581253512E+00;
    weight[ 3-1] = 0.28681737948564715225E+00;
    weight[ 4-1] = 0.072694951083385267446E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = 0.029023729768913933432E+00;
    xtab[2-1] = 0.14828524404581819442E+00;
    xtab[3-1] = 0.34531111151664787488E+00;
    xtab[4-1] = 0.59447696797658360178E+00;
    xtab[5-1] = 0.86538380686123504827E+00;
    xtab[6-1] = 1.1263076093187456632E+00;
    xtab[7-1] = 1.3470150460281258016E+00;
    xtab[8-1] = 1.5015603622059195568E+00;

    weight[ 1-1] = 0.073908998095117384985E+00;
    weight[ 2-1] = 0.16002993702338006099E+00;
    weight[ 3-1] = 0.21444434341803549108E+00;
    weight[ 4-1] = 0.21979581268851903339E+00;
    weight[ 5-1] = 0.17581164478209568886E+00;
    weight[ 6-1] = 0.10560448025308322171E+00;
    weight[ 7-1] = 0.042485497299217201089E+00;
    weight[ 8-1] = 0.0079192864405519178899E+00;
  }
  else if ( order == 16 )
  {
    xtab[ 1-1] = 0.0080145034906295973494E+00;
    xtab[ 2-1] = 0.041893031354246254797E+00;
    xtab[ 3-1] = 0.10149954486757579459E+00;
    xtab[ 4-1] = 0.18463185923836617507E+00;
    xtab[ 5-1] = 0.28826388487760574589E+00;
    xtab[ 6-1] = 0.40870579076464794191E+00;
    xtab[ 7-1] = 0.54176054986913847463E+00;
    xtab[ 8-1] = 0.68287636658719416893E+00;
    xtab[ 9-1] = 0.82729287620416833520E+00;
    xtab[10-1] = 0.97018212594829367065E+00;
    xtab[11-1] = 1.1067865150286247873E+00;
    xtab[12-1] = 1.2325555697227748824E+00;
    xtab[13-1] = 1.3432821921580721861E+00;
    xtab[14-1] = 1.4352370549295032923E+00;
    xtab[15-1] = 1.5052970876794669248E+00;
    xtab[16-1] = 1.5510586944086135769E+00;

    weight[ 1-1] = 0.020528714977215248902E+00;
    weight[ 2-1] = 0.046990919853597958123E+00;
    weight[ 3-1] = 0.071441021312218541698E+00;
    weight[ 4-1] = 0.092350338329243052271E+00;
    weight[ 5-1] = 0.10804928026816236935E+00;
    weight[ 6-1] = 0.11698241243306261791E+00;
    weight[ 7-1] = 0.11812395361762037649E+00;
    weight[ 8-1] = 0.11137584940420091049E+00;
    weight[ 9-1] = 0.097778236145946543110E+00;
    weight[10-1] = 0.079418758985944482077E+00;
    weight[11-1] = 0.059039620053768691402E+00;
    weight[12-1] = 0.039458876783728165671E+00;
    weight[13-1] = 0.022987785677206847531E+00;
    weight[14-1] = 0.011010405600421536861E+00;
    weight[15-1] = 0.0038123928030499915653E+00;
    weight[16-1] = 0.00065143375461266656171E+00;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET_COS2 - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void legendre_set_log ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_LOG sets a Gauss-Legendre rule for - LOG(X) * F(X) on [0,1].
  
  Discussion:
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w(x-1] = -log(x);
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) - LOG(X) * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    03 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Gwynne Evans,
    Practical Numerical Integration,
    Wiley, 1993, QA299.3E93, page 309.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 through 8, or 16.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] = 0.25E+00;

    weight[1-1] = 1.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = 0.112008806166976182957205488948E+00;
    xtab[2-1] = 0.602276908118738102757080225338E+00;

    weight[1-1] = 0.718539319030384440665510200891E+00;
    weight[2-1] = 0.281460680969615559334489799109E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = 0.0638907930873254049961166031363E+00;
    xtab[2-1] = 0.368997063715618765546197645857E+00;
    xtab[3-1] = 0.766880303938941455423682659817E+00;

    weight[1-1] = 0.513404552232363325129300497567E+00;
    weight[2-1] = 0.391980041201487554806287180966E+00;
    weight[3-1] = 0.0946154065661491200644123214672E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = 0.0414484801993832208033213101564E+00;
    xtab[2-1] = 0.245274914320602251939675759523E+00;
    xtab[3-1] = 0.556165453560275837180184354376E+00;
    xtab[4-1] = 0.848982394532985174647849188085E+00;

    weight[1-1] = 0.383464068145135124850046522343E+00;
    weight[2-1] = 0.386875317774762627336008234554E+00;
    weight[3-1] = 0.190435126950142415361360014547E+00;
    weight[4-1] = 0.0392254871299598324525852285552E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = 0.0291344721519720533037267621154E+00;
    xtab[2-1] = 0.173977213320897628701139710829E+00;
    xtab[3-1] = 0.411702520284902043174931924646E+00;
    xtab[4-1] = 0.677314174582820380701802667998E+00;
    xtab[5-1] = 0.894771361031008283638886204455E+00;

    weight[1-1] = 0.297893471782894457272257877888E+00;
    weight[2-1] = 0.349776226513224180375071870307E+00;
    weight[3-1] = 0.234488290044052418886906857943E+00;
    weight[4-1] = 0.0989304595166331469761807114404E+00;
    weight[5-1] = 0.0189115521431957964895826824218E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = 0.0216340058441169489956958558537E+00;
    xtab[2-1] = 0.129583391154950796131158505009E+00;
    xtab[3-1] = 0.314020449914765508798248188420E+00;
    xtab[4-1] = 0.538657217351802144548941893993E+00;
    xtab[5-1] = 0.756915337377402852164544156139E+00;
    xtab[6-1] = 0.922668851372120237333873231507E+00;

    weight[1-1] = 0.238763662578547569722268303330E+00;
    weight[2-1] = 0.308286573273946792969383109211E+00;
    weight[3-1] = 0.245317426563210385984932540188E+00;
    weight[4-1] = 0.142008756566476685421345576030E+00;
    weight[5-1] = 0.0554546223248862900151353549662E+00;
    weight[6-1] = 0.0101689586929322758869351162755E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = 0.0167193554082585159416673609320E+00;
    xtab[2-1] = 0.100185677915675121586885031757E+00;
    xtab[3-1] = 0.246294246207930599046668547239E+00;
    xtab[4-1] = 0.433463493257033105832882482601E+00;
    xtab[5-1] = 0.632350988047766088461805812245E+00;
    xtab[6-1] = 0.811118626740105576526226796782E+00;
    xtab[7-1] = 0.940848166743347721760134113379E+00;

    weight[1-1] = 0.196169389425248207525427377585E+00;
    weight[2-1] = 0.270302644247272982145271719533E+00;
    weight[3-1] = 0.239681873007690948308072785252E+00;
    weight[4-1] = 0.165775774810432906560869687736E+00;
    weight[5-1] = 0.0889432271376579644357238403458E+00;
    weight[6-1] = 0.0331943043565710670254494111034E+00;
    weight[7-1] = 0.00593278701512592399918517844468E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = 0.0133202441608924650122526725243E+00;
    xtab[2-1] = 0.0797504290138949384098277291424E+00;
    xtab[3-1] = 0.197871029326188053794476159516E+00;
    xtab[4-1] = 0.354153994351909419671463603538E+00;
    xtab[5-1] = 0.529458575234917277706149699996E+00;
    xtab[6-1] = 0.701814529939099963837152670310E+00;
    xtab[7-1] = 0.849379320441106676048309202301E+00;
    xtab[8-1] = 0.953326450056359788767379678514E+00;

    weight[1-1] = 0.164416604728002886831472568326E+00;
    weight[2-1] = 0.237525610023306020501348561960E+00;
    weight[3-1] = 0.226841984431919126368780402936E+00;
    weight[4-1] = 0.175754079006070244988056212006E+00;
    weight[5-1] = 0.112924030246759051855000442086E+00;
    weight[6-1] = 0.0578722107177820723985279672940E+00;
    weight[7-1] = 0.0209790737421329780434615241150E+00;
    weight[8-1] = 0.00368640710402761901335232127647E+00;
  }
  else if ( order == 16 )
  {
    xtab[ 1-1] = 0.00389783448711591592405360527037E+00;
    xtab[ 2-1] = 0.0230289456168732398203176309848E+00;
    xtab[ 3-1] = 0.0582803983062404123483532298394E+00;
    xtab[ 4-1] = 0.108678365091054036487713613051E+00;
    xtab[ 5-1] = 0.172609454909843937760843776232E+00;
    xtab[ 6-1] = 0.247937054470578495147671753047E+00;
    xtab[ 7-1] = 0.332094549129917155984755859320E+00;
    xtab[ 8-1] = 0.422183910581948600115088366745E+00;
    xtab[ 9-1] = 0.515082473381462603476277704052E+00;
    xtab[10-1] = 0.607556120447728724086384921709E+00;
    xtab[11-1] = 0.696375653228214061156318166581E+00;
    xtab[12-1] = 0.778432565873265405203868167732E+00;
    xtab[13-1] = 0.850850269715391083233822761319E+00;
    xtab[14-1] = 0.911086857222271905418818994060E+00;
    xtab[15-1] = 0.957025571703542157591520509383E+00;
    xtab[16-1] = 0.987047800247984476758697436516E+00;

    weight[ 1-1] = 0.0607917100435912328511733871235E+00;
    weight[ 2-1] = 0.102915677517582144387691736210E+00;
    weight[ 3-1] = 0.122355662046009193557547513197E+00;
    weight[ 4-1] = 0.127569246937015988717042209239E+00;
    weight[ 5-1] = 0.123013574600070915423123365137E+00;
    weight[ 6-1] = 0.111847244855485722621848903429E+00;
    weight[ 7-1] = 0.0965963851521243412529681650802E+00;
    weight[ 8-1] = 0.0793566643514731387824416770520E+00;
    weight[ 9-1] = 0.0618504945819652070951360853113E+00;
    weight[10-1] = 0.0454352465077266686288299526509E+00;
    weight[11-1] = 0.0310989747515818064092528627927E+00;
    weight[12-1] = 0.0194597659273608420780860268669E+00;
    weight[13-1] = 0.0107762549632055256455393162159E+00;
    weight[14-1] = 0.00497254289008764171250524951646E+00;
    weight[15-1] = 0.00167820111005119451503546419059E+00;
    weight[16-1] = 0.000282353764668436321778085987413E+00;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET_LOG - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void legendre_set_sqrtx_01 ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_SQRTX_01 sets Gauss-Legendre rules for SQRT(X)*F(X) on [0,1].
  
  Discussion:
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w(x) = sqrt ( x ).
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) SQRT ( X ) * F(X) dX =
      Integral ( 0 <= Y <= 1 ) 2 * Y**2 * F(Y**2) dY.
      (using Y = SQRT(X) )
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    06 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  int i;
  int order2;

  order2 = 2 * order + 1;

  double xtab2[order2];
  double weight2[order2];

  legendre_set ( order2, xtab2, weight2 );

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = pow ( xtab2[order+1+i], 2 );
  }
  for ( i = 0; i < order; i++ )
  {
    weight[i] = 2.0 * weight2[order+1+i] * xtab[i];
  }

  return;
}
/******************************************************************************/

void legendre_set_sqrtx2_01 ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_SQRTX2_01: Gauss-Legendre rules for F(X)/SQRT(X) on [0,1].
  
  Discussion:
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w(x) = 1 / sqrt ( x ).
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) F(X) / SQRT ( X ) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    06 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  int i;
  int order2;

  order2 = 2 * order + 1;

   double xtab2[order2];
   double weight2[order2];

  legendre_set ( order2, xtab2, weight2 );

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = pow ( xtab2[order+1+i], 2 );
  }
  for ( i = 0; i < order; i++ )
  {
    weight[i] = 2.0 * weight2[order+1+i];
  }

  return;
}
/******************************************************************************/

void legendre_set_x0_01 ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_X0_01 sets a Gauss-Legendre rule for F(X) on [0,1].
  
  Discussion:
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w(x-1] = 1.0.
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    03 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 8.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] =   0.5E+00;

    weight[1-1] = 1.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = 0.2113248654E+00;
    xtab[2-1] = 0.7886751346E+00;

    weight[1-1] = 0.5E+00;
    weight[2-1] = 0.5E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = 0.1127016654E+00;
    xtab[2-1] = 0.5000000000E+00;
    xtab[3-1] = 0.8872983346E+00;

    weight[1-1] = 5.0E+00 / 18.0E+00;
    weight[2-1] = 8.0E+00 / 18.0E+00;
    weight[3-1] = 5.0E+00 / 18.0E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = 0.0694318442E+00;
    xtab[2-1] = 0.3300094782E+00;
    xtab[3-1] = 0.6699905218E+00;
    xtab[4-1] = 0.9305681558E+00;

    weight[1-1] = 0.1739274226E+00;
    weight[2-1] = 0.3260725774E+00;
    weight[3-1] = 0.3260725774E+00;
    weight[4-1] = 0.1739274226E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = 0.0469100770E+00;
    xtab[2-1] = 0.2307653449E+00;
    xtab[3-1] = 0.5000000000E+00;
    xtab[4-1] = 0.7692346551E+00;
    xtab[5-1] = 0.9530899230E+00;

    weight[1-1] = 0.1184634425E+00;
    weight[2-1] = 0.2393143352E+00;
    weight[3-1] = 0.2844444444E+00;
    weight[4-1] = 0.2393143352E+00;
    weight[5-1] = 0.1184634425E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = 0.0337652429E+00;
    xtab[2-1] = 0.1693953068E+00;
    xtab[3-1] = 0.3806904070E+00;
    xtab[4-1] = 0.6193095930E+00;
    xtab[5-1] = 0.8306046932E+00;
    xtab[6-1] = 0.9662347571E+00;

    weight[1-1] = 0.0856622462E+00;
    weight[2-1] = 0.1803807865E+00;
    weight[3-1] = 0.2339569673E+00;
    weight[4-1] = 0.2339569673E+00;
    weight[5-1] = 0.1803807865E+00;
    weight[6-1] = 0.0856622462E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = 0.0254460438E+00;
    xtab[2-1] = 0.1292344072E+00;
    xtab[3-1] = 0.2970774243E+00;
    xtab[4-1] = 0.5000000000E+00;
    xtab[5-1] = 0.7029225757E+00;
    xtab[6-1] = 0.8707655928E+00;
    xtab[7-1] = 0.9745539562E+00;

    weight[1-1] = 0.0647424831E+00;
    weight[2-1] = 0.1398526957E+00;
    weight[3-1] = 0.1909150253E+00;
    weight[4-1] = 0.2089795918E+00;
    weight[5-1] = 0.1909150253E+00;
    weight[6-1] = 0.1398526957E+00;
    weight[7-1] = 0.0647424831E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = 0.0198550718E+00;
    xtab[2-1] = 0.1016667613E+00;
    xtab[3-1] = 0.2372337950E+00;
    xtab[4-1] = 0.4082826788E+00;
    xtab[5-1] = 0.5917173212E+00;
    xtab[6-1] = 0.7627662050E+00;
    xtab[7-1] = 0.8983332387E+00;
    xtab[8-1] = 0.9801449282E+00;

    weight[1-1] = 0.0506142681E+00;
    weight[2-1] = 0.1111905172E+00;
    weight[3-1] = 0.1568533229E+00;
    weight[4-1] = 0.1813418917E+00;
    weight[5-1] = 0.1813418917E+00;
    weight[6-1] = 0.1568533229E+00;
    weight[7-1] = 0.1111905172E+00;
    weight[8-1] = 0.0506142681E+00;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET_X0_01 - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void legendre_set_x1 ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_X1 sets a Gauss-Legendre rule for ( 1 + X ) * F(X) on [-1,1].
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x-1] = 1 + x.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) ( 1 + X ) * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    03 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 9.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] =  0.333333333333333333333333333333E+00;

    weight[1-1] = 2.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = -0.289897948556635619639456814941E+00;
    xtab[2-1] =  0.689897948556635619639456814941E+00;

    weight[1-1] =  0.727834473024091322422523991699E+00;
    weight[2-1] =  1.27216552697590867757747600830E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = -0.575318923521694112050483779752E+00;
    xtab[2-1] =  0.181066271118530578270147495862E+00;
    xtab[3-1] =  0.822824080974592105208907712461E+00;

    weight[1-1] =  0.279307919605816490135525088716E+00;
    weight[2-1] =  0.916964425438344986775682378225E+00;
    weight[3-1] =  0.803727654955838523088792533058E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = -0.720480271312438895695825837750E+00;
    xtab[2-1] = -0.167180864737833640113395337326E+00;
    xtab[3-1] =  0.446313972723752344639908004629E+00;
    xtab[4-1] =  0.885791607770964635613757614892E+00;

    weight[1-1] =  0.124723883800032328695500588386E+00;
    weight[2-1] =  0.519390190432929763305824811559E+00;
    weight[3-1] =  0.813858272041085443165617903743E+00;
    weight[4-1] =  0.542027653725952464833056696312E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = -0.802929828402347147753002204224E+00;
    xtab[2-1] = -0.390928546707272189029229647442E+00;
    xtab[3-1] =  0.124050379505227711989974959990E+00;
    xtab[4-1] =  0.603973164252783654928415726409E+00;
    xtab[5-1] =  0.920380285897062515318386619813E+00;

    weight[1-1] =  0.0629916580867691047411692662740E+00;
    weight[2-1] =  0.295635480290466681402532877367E+00;
    weight[3-1] =  0.585547948338679234792151477424E+00;
    weight[4-1] =  0.668698552377478261966702492391E+00;
    weight[5-1] =  0.387126360906606717097443886545E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = -0.853891342639482229703747931639E+00;
    xtab[2-1] = -0.538467724060109001833766720231E+00;
    xtab[3-1] = -0.117343037543100264162786683611E+00;
    xtab[4-1] =  0.326030619437691401805894055838E+00;
    xtab[5-1] =  0.703842800663031416300046295008E+00;
    xtab[6-1] =  0.941367145680430216055899446174E+00;

    weight[1-1] =  0.0349532072544381270240692132496E+00;
    weight[2-1] =  0.175820662202035902032706497222E+00;
    weight[3-1] =  0.394644603562621056482338042193E+00;
    weight[4-1] =  0.563170215152795712476307356284E+00;
    weight[5-1] =  0.542169988926074467362761586552E+00;
    weight[6-1] =  0.289241322902034734621817304499E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = -0.887474878926155707068695617935E+00;
    xtab[2-1] = -0.639518616526215270024840114382E+00;
    xtab[3-1] = -0.294750565773660725252184459658E+00;
    xtab[4-1] =  0.0943072526611107660028971153047E+00;
    xtab[5-1] =  0.468420354430821063046421216613E+00;
    xtab[6-1] =  0.770641893678191536180719525865E+00;
    xtab[7-1] =  0.955041227122575003782349000858E+00;

    weight[1-1] =  0.0208574488112296163587654972151E+00;
    weight[2-1] =  0.109633426887493901777324193433E+00;
    weight[3-1] =  0.265538785861965879934591955055E+00;
    weight[4-1] =  0.428500262783494679963649011999E+00;
    weight[5-1] =  0.509563589198353307674937943100E+00;
    weight[6-1] =  0.442037032763498409684482945478E+00;
    weight[7-1] =  0.223869453693964204606248453720E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = -0.910732089420060298533757956283E+00;
    xtab[2-1] = -0.711267485915708857029562959544E+00;
    xtab[3-1] = -0.426350485711138962102627520502E+00;
    xtab[4-1] = -0.0903733696068532980645444599064E+00;
    xtab[5-1] =  0.256135670833455395138292079035E+00;
    xtab[6-1] =  0.571383041208738483284917464837E+00;
    xtab[7-1] =  0.817352784200412087992517083851E+00;
    xtab[8-1] =  0.964440169705273096373589797925E+00;

    weight[1-1] =  0.0131807657689951954189692640444E+00;
    weight[2-1] =  0.0713716106239448335742111888042E+00;
    weight[3-1] =  0.181757278018795592332221684383E+00;
    weight[4-1] =  0.316798397969276640481632757440E+00;
    weight[5-1] =  0.424189437743720042818124385645E+00;
    weight[6-1] =  0.450023197883549464687088394417E+00;
    weight[7-1] =  0.364476094545494505382889847132E+00;
    weight[8-1] =  0.178203217446223725304862478136E+00;
  }
  else if ( order == 9 )
  {
    xtab[1-1] = -0.927484374233581078117671398464E+00;
    xtab[2-1] = -0.763842042420002599615429776011E+00;
    xtab[3-1] = -0.525646030370079229365386614293E+00;
    xtab[4-1] = -0.236234469390588049278459503207E+00;
    xtab[5-1] =  0.0760591978379781302337137826389E+00;
    xtab[6-1] =  0.380664840144724365880759065541E+00;
    xtab[7-1] =  0.647766687674009436273648507855E+00;
    xtab[8-1] =  0.851225220581607910728163628088E+00;
    xtab[9-1] =  0.971175180702246902734346518378E+00;

    weight[1-1] =  0.00872338834309252349019620448007E+00;
    weight[2-1] =  0.0482400171391415162069086091476E+00;
    weight[3-1] =  0.127219285964216005046760427743E+00;
    weight[4-1] =  0.233604781180660442262926091607E+00;
    weight[5-1] =  0.337433287379681397577000079834E+00;
    weight[6-1] =  0.401235236773473158616600898930E+00;
    weight[7-1] =  0.394134968689382820640692081477E+00;
    weight[8-1] =  0.304297020437232650320317215016E+00;
    weight[9-1] =  0.145112014093119485838598391765E+00;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET_X1 - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of ORDER = %d\n", order );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void legendre_set_x1_01 ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_X1_01 sets a Gauss-Legendre rule for X * F(X) on [0,1].
  
  Discussion:
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w(x-1] = x.
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) X * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    10 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 8.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] =   0.6666666667E+00;

    weight[1-1] = 0.5000000000E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = 0.3550510257E+00;
    xtab[2-1] = 0.8449489743E+00;

    weight[1-1] = 0.1819586183E+00;
    weight[2-1] = 0.3180413817E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = 0.2123405382E+00;
    xtab[2-1] = 0.5905331356E+00;
    xtab[3-1] = 0.9114120405E+00;

    weight[1-1] = 0.0698269799E+00;
    weight[2-1] = 0.2292411064E+00;
    weight[3-1] = 0.2009319137E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = 0.1397598643E+00;
    xtab[2-1] = 0.4164095676E+00;
    xtab[3-1] = 0.7231569864E+00;
    xtab[4-1] = 0.9428958039E+00;

    weight[1-1] = 0.0311809710E+00;
    weight[2-1] = 0.1298475476E+00;
    weight[3-1] = 0.2034645680E+00;
    weight[4-1] = 0.1355069134E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = 0.0985350858E+00;
    xtab[2-1] = 0.3045357266E+00;
    xtab[3-1] = 0.5620251898E+00;
    xtab[4-1] = 0.8019865821E+00;
    xtab[5-1] = 0.9601901429E+00;

    weight[1-1] = 0.0157479145E+00;
    weight[2-1] = 0.0739088701E+00;
    weight[3-1] = 0.1463869871E+00;
    weight[4-1] = 0.1671746381E+00;
    weight[5-1] = 0.0967815902E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = 0.0730543287E+00;
    xtab[2-1] = 0.2307661380E+00;
    xtab[3-1] = 0.4413284812E+00;
    xtab[4-1] = 0.6630153097E+00;
    xtab[5-1] = 0.8519214003E+00;
    xtab[6-1] = 0.9706835728E+00;

    weight[1-1] = 0.0087383108E+00;
    weight[2-1] = 0.0439551656E+00;
    weight[3-1] = 0.0986611509E+00;
    weight[4-1] = 0.1407925538E+00;
    weight[5-1] = 0.1355424972E+00;
    weight[6-1] = 0.0723103307E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = 0.0562625605E+00;
    xtab[2-1] = 0.1802406917E+00;
    xtab[3-1] = 0.3526247171E+00;
    xtab[4-1] = 0.5471536263E+00;
    xtab[5-1] = 0.7342101772E+00;
    xtab[6-1] = 0.8853209468E+00;
    xtab[7-1] = 0.9775206136E+00;

    weight[1-1] = 0.0052143622E+00;
    weight[2-1] = 0.0274083567E+00;
    weight[3-1] = 0.0663846965E+00;
    weight[4-1] = 0.1071250657E+00;
    weight[5-1] = 0.1273908973E+00;
    weight[6-1] = 0.1105092582E+00;
    weight[7-1] = 0.0559673634E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = 0.0446339553E+00;
    xtab[2-1] = 0.1443662570E+00;
    xtab[3-1] = 0.2868247571E+00;
    xtab[4-1] = 0.4548133152E+00;
    xtab[5-1] = 0.6280678354E+00;
    xtab[6-1] = 0.7856915206E+00;
    xtab[7-1] = 0.9086763921E+00;
    xtab[8-1] = 0.9822200849E+00;

    weight[1-1] = 0.0032951914E+00;
    weight[2-1] = 0.0178429027E+00;
    weight[3-1] = 0.0454393195E+00;
    weight[4-1] = 0.0791995995E+00;
    weight[5-1] = 0.1060473594E+00;
    weight[6-1] = 0.1125057995E+00;
    weight[7-1] = 0.0911190236E+00;
    weight[8-1] = 0.0445508044E+00;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET_X1_01 - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void legendre_set_x2 ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_X2 sets Gauss-Legendre rules for ( 1 + X )^2*F(X) on [-1,1].
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x-1] = ( 1 + x )**2.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) ( 1 + X )**2 * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    02 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 9.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] =  0.5E+00;

    weight[1-1] =  2.66666666666666666666666666666E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = -0.0883036880224505775998524725910E+00;
    xtab[2-1] =  0.754970354689117244266519139258E+00;

    weight[1-1] =  0.806287056638603444666851075928E+00;
    weight[2-1] =  1.86037961002806322199981559074E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = -0.410004419776996766244796955168E+00;
    xtab[2-1] =  0.305992467923296230556472913192E+00;
    xtab[3-1] =  0.854011951853700535688324041976E+00;

    weight[1-1] =  0.239605624068645584091811926047E+00;
    weight[2-1] =  1.16997015407892817602809616291E+00;
    weight[3-1] =  1.25709088851909290654675857771E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = -0.591702835793545726606755921586E+00;
    xtab[2-1] = -0.0340945902087350046811467387661E+00;
    xtab[3-1] =  0.522798524896275389882037174551E+00;
    xtab[4-1] =  0.902998901106005341405865485802E+00;

    weight[1-1] =  0.0828179259993445222751812523731E+00;
    weight[2-1] =  0.549071097383384602539010760334E+00;
    weight[3-1] =  1.14767031839371367238662411421E+00;
    weight[4-1] =  0.887107324890223869465850539752E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = -0.702108425894032836232448374820E+00;
    xtab[2-1] = -0.268666945261773544694327777841E+00;
    xtab[3-1] =  0.220227225868961343518209179230E+00;
    xtab[4-1] =  0.653039358456608553790815164028E+00;
    xtab[5-1] =  0.930842120163569816951085142737E+00;

    weight[1-1] =  0.0329106016247920636689299329544E+00;
    weight[2-1] =  0.256444805783695354037991444453E+00;
    weight[3-1] =  0.713601289772720001490035944563E+00;
    weight[4-1] =  1.00959169519929190423066348132E+00;
    weight[5-1] =  0.654118274286167343239045863379E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = -0.773611232355123732602532012021E+00;
    xtab[2-1] = -0.431362254623427837535325249187E+00;
    xtab[3-1] = -0.0180728263295041680220798103354E+00;
    xtab[4-1] =  0.395126163954217534500188844163E+00;
    xtab[5-1] =  0.736872116684029732026178298518E+00;
    xtab[6-1] =  0.948190889812665614490712786006E+00;

    weight[1-1] =  0.0146486064549543818622276447204E+00;
    weight[2-1] =  0.125762377479560410622810097040E+00;
    weight[3-1] =  0.410316569036929681761034600615E+00;
    weight[4-1] =  0.756617493988329628546336413760E+00;
    weight[5-1] =  0.859011997894245060846045458784E+00;
    weight[6-1] =  0.500309621812647503028212451747E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = -0.822366333126005527278634734418E+00;
    xtab[2-1] = -0.547034493182875002223997992852E+00;
    xtab[3-1] = -0.200043026557985860387937545780E+00;
    xtab[4-1] =  0.171995710805880507163425502299E+00;
    xtab[5-1] =  0.518891747903884926692601716998E+00;
    xtab[6-1] =  0.793821941703901970495546427988E+00;
    xtab[7-1] =  0.959734452453198985538996625765E+00;

    weight[1-1] =  0.00714150426951365443207221475404E+00;
    weight[2-1] =  0.0653034050584375560578544725498E+00;
    weight[3-1] =  0.235377690316228918725962815880E+00;
    weight[4-1] =  0.505171029671130381676271523850E+00;
    weight[5-1] =  0.733870426238362032891332767175E+00;
    weight[6-1] =  0.725590596901489156295739839779E+00;
    weight[7-1] =  0.394212014211504966587433032679E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = -0.857017929919813794402037235698E+00;
    xtab[2-1] = -0.631543407166567521509503573952E+00;
    xtab[3-1] = -0.339104543648722903660229021109E+00;
    xtab[4-1] = -0.0111941563689783438801237300122E+00;
    xtab[5-1] =  0.316696017045595559454075475675E+00;
    xtab[6-1] =  0.609049663022520165351466780939E+00;
    xtab[7-1] =  0.834198765028697794599267293239E+00;
    xtab[8-1] =  0.967804480896157932935972899807E+00;

    weight[1-1] =  0.00374814227227757804631954025851E+00;
    weight[2-1] =  0.0357961737041152639660521680263E+00;
    weight[3-1] =  0.137974910241879862433949246199E+00;
    weight[4-1] =  0.326515411108352185491692769217E+00;
    weight[5-1] =  0.547577467373226177976217604887E+00;
    weight[6-1] =  0.682278153375510121675529810121E+00;
    weight[7-1] =  0.614544746137780998436053880546E+00;
    weight[8-1] =  0.318231662453524478640851647411E+00;
  }
  else if ( order == 9 )
  {
    xtab[1-1] = -0.882491728426548422828684254270E+00;
    xtab[2-1] = -0.694873684026474640346360850039E+00;
    xtab[3-1] = -0.446537143480670863635920316400E+00;
    xtab[4-1] = -0.159388112702326252531544826624E+00;
    xtab[5-1] =  0.141092709224374414981503995427E+00;
    xtab[6-1] =  0.428217823321559204544020866175E+00;
    xtab[7-1] =  0.676480966471850715860378175342E+00;
    xtab[8-1] =  0.863830940812464825046988286026E+00;
    xtab[9-1] =  0.973668228805771018909618924364E+00;

    weight[1-1] =  0.00209009877215570354392734918986E+00;
    weight[2-1] =  0.0205951891648697848186537272448E+00;
    weight[3-1] =  0.0832489326348178964194106978875E+00;
    weight[4-1] =  0.210746247220398685903797568021E+00;
    weight[5-1] =  0.388325022916052063676224499399E+00;
    weight[6-1] =  0.554275165518437673725822282791E+00;
    weight[7-1] =  0.621388553284444032628761363828E+00;
    weight[8-1] =  0.523916296267173054255512857631E+00;
    weight[9-1] =  0.262081160888317771694556320674E+00;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET_X2 - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of ORDER = %d\n", order );
    exit ( 1 );
  }
  return;
}
/******************************************************************************/

void legendre_set_x2_01 ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LEGENDRE_SET_X2_01 sets a Gauss-Legendre rule for X**2 * F(X) on [0,1].
  
  Discussion:
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w(x-1] = x**2.
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) X*X * F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    01 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 8.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] =   0.75E+00;

    weight[1-1] = 1.0E+00 / 3.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] = 0.4558481560E+00;
    xtab[2-1] = 0.8774851773E+00;

    weight[1-1] = 0.1007858821E+00;
    weight[2-1] = 0.2325474513E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = 0.2949977901E+00;
    xtab[2-1] = 0.6529962340E+00;
    xtab[3-1] = 0.9270059759E+00;

    weight[1-1] = 0.0299507030E+00;
    weight[2-1] = 0.1462462693E+00;
    weight[3-1] = 0.1571363611E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] = 0.2041485821E+00;
    xtab[2-1] = 0.4829527049E+00;
    xtab[3-1] = 0.7613992624E+00;
    xtab[4-1] = 0.9514994506E+00;

    weight[1-1] = 0.0103522408E+00;
    weight[2-1] = 0.0686338872E+00;
    weight[3-1] = 0.1434587898E+00;
    weight[4-1] = 0.1108884156E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] = 0.1489457871E+00;
    xtab[2-1] = 0.3656665274E+00;
    xtab[3-1] = 0.6101136129E+00;
    xtab[4-1] = 0.8265196792E+00;
    xtab[5-1] = 0.9654210601E+00;

    weight[1-1] = 0.0041138252E+00;
    weight[2-1] = 0.0320556007E+00;
    weight[3-1] = 0.0892001612E+00;
    weight[4-1] = 0.1261989619E+00;
    weight[5-1] = 0.0817647843E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] = 0.1131943838E+00;
    xtab[2-1] = 0.2843188727E+00;
    xtab[3-1] = 0.4909635868E+00;
    xtab[4-1] = 0.6975630820E+00;
    xtab[5-1] = 0.8684360583E+00;
    xtab[6-1] = 0.9740954449E+00;

    weight[1-1] = 0.0018310758E+00;
    weight[2-1] = 0.0157202972E+00;
    weight[3-1] = 0.0512895711E+00;
    weight[4-1] = 0.0945771867E+00;
    weight[5-1] = 0.1073764997E+00;
    weight[6-1] = 0.0625387027E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = 0.0888168334E+00;
    xtab[2-1] = 0.2264827534E+00;
    xtab[3-1] = 0.3999784867E+00;
    xtab[4-1] = 0.5859978554E+00;
    xtab[5-1] = 0.7594458740E+00;
    xtab[6-1] = 0.8969109709E+00;
    xtab[7-1] = 0.9798672262E+00;

    weight[1-1] = 0.0008926880E+00;
    weight[2-1] = 0.0081629256E+00;
    weight[3-1] = 0.0294222113E+00;
    weight[4-1] = 0.0631463787E+00;
    weight[5-1] = 0.0917338033E+00;
    weight[6-1] = 0.0906988246E+00;
    weight[7-1] = 0.0492765018E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] = 0.0714910350E+00;
    xtab[2-1] = 0.1842282964E+00;
    xtab[3-1] = 0.3304477282E+00;
    xtab[4-1] = 0.4944029218E+00;
    xtab[5-1] = 0.6583480085E+00;
    xtab[6-1] = 0.8045248315E+00;
    xtab[7-1] = 0.9170993825E+00;
    xtab[8-1] = 0.9839022404E+00;

    weight[1-1] = 0.0004685178E+00;
    weight[2-1] = 0.0044745217E+00;
    weight[3-1] = 0.0172468638E+00;
    weight[4-1] = 0.0408144264E+00;
    weight[5-1] = 0.0684471834E+00;
    weight[6-1] = 0.0852847692E+00;
    weight[7-1] = 0.0768180933E+00;
    weight[8-1] = 0.0397789578E+00;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LEGENDRE_SET_X2_01 - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void lobatto_compute ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    LOBATTO_COMPUTE computes a Lobatto quadrature rule.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x) = 1.0.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) )
  
    The quadrature rule will integrate exactly all polynomials up to
    X**(2*N-3).
  
    The Lobatto rule is distinguished by the fact that both endpoints
    (-1 and 1) are always abscissas of the rule.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    04 February 2007
  
  Author:
  
    Original MATLAB version by Greg von Winckel.
    C version by John Burkardt.
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
    Spectral Methods in Fluid Dynamics,
    Springer, 1993,
    ISNB13: 978-3540522058,
    LC: QA377.S676.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int N, the order of the rule.  N must be at least 2.
  
    Output, double X[N], W[N], the abscissas and weights
    of the rule.
*/
{
  int i;
  int j;
  double pi = 3.141592653589793;
  double test;
  double tolerance;

  if ( n < 2 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LOBATTO_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", n );
    fprintf ( stderr, "  ORDER must be at least 2.\n" );
    exit ( 1 );
  }
  tolerance = 100.0 * r8_epsilon ( );
/*
  Initial estimate for the abscissas is the Chebyshev-Gauss-Lobatto nodes.
*/
  for ( i = 0; i < n; i++ )
  {
    x[i] = cos ( pi * ( double ) ( i ) / ( double ) ( n - 1 ) );
  }

   double xold[n];
   double p[n*n];

  do
  {
    for ( i = 0; i < n; i++ )
    {
      xold[i] = x[i];
    }
    for ( i = 0; i < n; i++ )
    {
      p[i+0*n] = 1.0;
    }
    for ( i = 0; i < n; i++ )
    {
      p[i+1*n] = x[i];
    }

    for ( j = 2; j <= n-1; j++ )
    {
      for ( i = 0; i < n; i++)
      {
        p[i+j*n] = ( ( double ) ( 2 * j - 1 ) * x[i] * p[i+(j-1)*n]
                   + ( double ) (   - j + 1 ) *        p[i+(j-2)*n] )
                   / ( double ) (     j     );
      }
    }

    for ( i = 0; i < n; i++ )
    {
      x[i] = xold[i] - ( x[i] * p[i+(n-1)*n] - p[i+(n-2)*n] )
           / ( ( double ) ( n ) * p[i+(n-1)*n] );
    }

    test = 0.0;
    for ( i = 0; i < n; i++ )
    {
      test = r8_max ( test, r8_abs ( x[i] - xold[i] ) );
    }

  } while ( tolerance < test );

  r8vec_reverse ( n, x );

  for ( i = 0; i < n; i++ )
  {
    w[i] = 2.0 / ( ( double ) ( ( n - 1 ) * n ) * pow ( p[i+(n-1)*n], 2 ) );
  }

  return;
}
/******************************************************************************/

void lobatto_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    LOBATTO_SET sets abscissas and weights for Lobatto quadrature.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x-1] = 1.0.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
  
    The quadrature rule will integrate exactly all polynomials up to
    X**(2*ORDER-3).
  
    The Lobatto rule is distinguished by the fact that both endpoints
    (-1 and 1) are always abscissas of the rule.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    30 April 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 2 and 20.
  
    Output, double XTAB[ORDER], the abscissas for the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are positive, symmetric and should sum to 2.
*/
{
  if ( order == 2 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =    1.0E+00;

    weight[1-1] =  1.0E+00;
    weight[2-1] =  1.0E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =    0.0E+00;
    xtab[3-1] =    1.0E+00;

    weight[1-1] =  1.0 / 3.0E+00;
    weight[2-1] =  4.0 / 3.0E+00;
    weight[3-1] =  1.0 / 3.0E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.447213595499957939281834733746E+00;
    xtab[3-1] =    0.447213595499957939281834733746E+00;
    xtab[4-1] =    1.0E+00;

    weight[1-1] =  1.0E+00 / 6.0E+00;
    weight[2-1] =  5.0E+00 / 6.0E+00;
    weight[3-1] =  5.0E+00 / 6.0E+00;
    weight[4-1] =  1.0E+00 / 6.0E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.654653670707977143798292456247E+00;
    xtab[3-1] =    0.0E+00;
    xtab[4-1] =    0.654653670707977143798292456247E+00;
    xtab[5-1] =    1.0E+00;

    weight[1-1] =  9.0E+00 / 90.0E+00;
    weight[2-1] = 49.0E+00 / 90.0E+00;
    weight[3-1] = 64.0E+00 / 90.0E+00;
    weight[4-1] = 49.0E+00 / 90.0E+00;
    weight[5-1] =  9.0E+00 / 90.0E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.765055323929464692851002973959E+00;
    xtab[3-1] =  - 0.285231516480645096314150994041E+00;
    xtab[4-1] =    0.285231516480645096314150994041E+00;
    xtab[5-1] =    0.765055323929464692851002973959E+00;
    xtab[6-1] =    1.0E+00;

    weight[1-1] =  0.066666666666666666666666666667E+00;
    weight[2-1] =  0.378474956297846980316612808212E+00;
    weight[3-1] =  0.554858377035486353016720525121E+00;
    weight[4-1] =  0.554858377035486353016720525121E+00;
    weight[5-1] =  0.378474956297846980316612808212E+00;
    weight[6-1] =  0.066666666666666666666666666667E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.830223896278566929872032213967E+00;
    xtab[3-1] =  - 0.468848793470714213803771881909E+00;
    xtab[4-1] =    0.0E+00;
    xtab[5-1] =    0.468848793470714213803771881909E+00;
    xtab[6-1] =    0.830223896278566929872032213967E+00;
    xtab[7-1] =    1.0E+00;

    weight[1-1] =  0.476190476190476190476190476190E-01;
    weight[2-1] =  0.276826047361565948010700406290E+00;
    weight[3-1] =  0.431745381209862623417871022281E+00;
    weight[4-1] =  0.487619047619047619047619047619E+00;
    weight[5-1] =  0.431745381209862623417871022281E+00;
    weight[6-1] =  0.276826047361565948010700406290E+00;
    weight[7-1] =  0.476190476190476190476190476190E-01;
  }
  else if ( order == 8 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.871740148509606615337445761221E+00;
    xtab[3-1] =  - 0.591700181433142302144510731398E+00;
    xtab[4-1] =  - 0.209299217902478868768657260345E+00;
    xtab[5-1] =    0.209299217902478868768657260345E+00;
    xtab[6-1] =    0.591700181433142302144510731398E+00;
    xtab[7-1] =    0.871740148509606615337445761221E+00;
    xtab[8-1] =    1.0E+00;

    weight[1-1] =  0.357142857142857142857142857143E-01;
    weight[2-1] =  0.210704227143506039382991065776E+00;
    weight[3-1] =  0.341122692483504364764240677108E+00;
    weight[4-1] =  0.412458794658703881567052971402E+00;
    weight[5-1] =  0.412458794658703881567052971402E+00;
    weight[6-1] =  0.341122692483504364764240677108E+00;
    weight[7-1] =  0.210704227143506039382991065776E+00;
    weight[8-1] =  0.357142857142857142857142857143E-01;
  }
  else if ( order == 9 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.899757995411460157312345244418E+00;
    xtab[3-1] =  - 0.677186279510737753445885427091E+00;
    xtab[4-1] =  - 0.363117463826178158710752068709E+00;
    xtab[5-1] =    0.0E+00;
    xtab[6-1] =    0.363117463826178158710752068709E+00;
    xtab[7-1] =    0.677186279510737753445885427091E+00;
    xtab[8-1] =    0.899757995411460157312345244418E+00;
    xtab[9-1] =    1.0E+00;

    weight[1-1] =  0.277777777777777777777777777778E-01;
    weight[2-1] =  0.165495361560805525046339720029E+00;
    weight[3-1] =  0.274538712500161735280705618579E+00;
    weight[4-1] =  0.346428510973046345115131532140E+00;
    weight[5-1] =  0.371519274376417233560090702948E+00;
    weight[6-1] =  0.346428510973046345115131532140E+00;
    weight[7-1] =  0.274538712500161735280705618579E+00;
    weight[8-1] =  0.165495361560805525046339720029E+00;
    weight[9-1] =  0.277777777777777777777777777778E-01;
  }
  else if ( order == 10 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.919533908166458813828932660822E+00;
    xtab[3-1] =  - 0.738773865105505075003106174860E+00;
    xtab[4-1] =  - 0.477924949810444495661175092731E+00;
    xtab[5-1] =  - 0.165278957666387024626219765958E+00;
    xtab[6-1] =    0.165278957666387024626219765958E+00;
    xtab[7-1] =    0.477924949810444495661175092731E+00;
    xtab[8-1] =    0.738773865105505075003106174860E+00;
    xtab[9-1] =    0.919533908166458813828932660822E+00;
    xtab[10-1] =   1.0E+00;

    weight[1-1] =  0.222222222222222222222222222222E-01;
    weight[2-1] =  0.133305990851070111126227170755E+00;
    weight[3-1] =  0.224889342063126452119457821731E+00;
    weight[4-1] =  0.292042683679683757875582257374E+00;
    weight[5-1] =  0.327539761183897456656510527917E+00;
    weight[6-1] =  0.327539761183897456656510527917E+00;
    weight[7-1] =  0.292042683679683757875582257374E+00;
    weight[8-1] =  0.224889342063126452119457821731E+00;
    weight[9-1] =  0.133305990851070111126227170755E+00;
    weight[10-1] = 0.222222222222222222222222222222E-01;
  }
  else if ( order == 11 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.934001430408059134332274136099E+00;
    xtab[3-1] =  - 0.784483473663144418622417816108E+00;
    xtab[4-1] =  - 0.565235326996205006470963969478E+00;
    xtab[5-1] =  - 0.295758135586939391431911515559E+00;
    xtab[6-1] =    0.0E+00;
    xtab[7-1] =    0.295758135586939391431911515559E+00;
    xtab[8-1] =    0.565235326996205006470963969478E+00;
    xtab[9-1] =    0.784483473663144418622417816108E+00;
    xtab[10-1] =   0.934001430408059134332274136099E+00;
    xtab[11-1] =   1.0E+00;

    weight[1-1] =  0.181818181818181818181818181818E-01;
    weight[2-1] =  0.109612273266994864461403449580E+00;
    weight[3-1] =  0.187169881780305204108141521899E+00;
    weight[4-1] =  0.248048104264028314040084866422E+00;
    weight[5-1] =  0.286879124779008088679222403332E+00;
    weight[6-1] =  0.300217595455690693785931881170E+00;
    weight[7-1] =  0.286879124779008088679222403332E+00;
    weight[8-1] =  0.248048104264028314040084866422E+00;
    weight[9-1] =  0.187169881780305204108141521899E+00;
    weight[10-1] = 0.109612273266994864461403449580E+00;
    weight[11-1] = 0.181818181818181818181818181818E-01;
  }
  else if ( order == 12 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.944899272222882223407580138303E+00;
    xtab[3-1] =  - 0.819279321644006678348641581717E+00;
    xtab[4-1] =  - 0.632876153031869677662404854444E+00;
    xtab[5-1] =  - 0.399530940965348932264349791567E+00;
    xtab[6-1] =  - 0.136552932854927554864061855740E+00;
    xtab[7-1] =    0.136552932854927554864061855740E+00;
    xtab[8-1] =    0.399530940965348932264349791567E+00;
    xtab[9-1] =    0.632876153031869677662404854444E+00;
    xtab[10-1] =   0.819279321644006678348641581717E+00;
    xtab[11-1] =   0.944899272222882223407580138303E+00;
    xtab[12-1] =   1.0E+00;

    weight[1-1] =  0.151515151515151515151515151515E-01;
    weight[2-1] =  0.916845174131961306683425941341E-01;
    weight[3-1] =  0.157974705564370115164671062700E+00;
    weight[4-1] =  0.212508417761021145358302077367E+00;
    weight[5-1] =  0.251275603199201280293244412148E+00;
    weight[6-1] =  0.271405240910696177000288338500E+00;
    weight[7-1] =  0.271405240910696177000288338500E+00;
    weight[8-1] =  0.251275603199201280293244412148E+00;
    weight[9-1] =  0.212508417761021145358302077367E+00;
    weight[10-1] = 0.157974705564370115164671062700E+00;
    weight[11-1] = 0.916845174131961306683425941341E-01;
    weight[12-1] = 0.151515151515151515151515151515E-01;
  }
  else if ( order == 13 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.953309846642163911896905464755E+00;
    xtab[3-1] =  - 0.846347564651872316865925607099E+00;
    xtab[4-1] =  - 0.686188469081757426072759039566E+00;
    xtab[5-1] =  - 0.482909821091336201746937233637E+00;
    xtab[6-1] =  - 0.249286930106239992568673700374E+00;
    xtab[7-1] =    0.0E+00;
    xtab[8-1] =    0.249286930106239992568673700374E+00;
    xtab[9-1] =    0.482909821091336201746937233637E+00;
    xtab[10-1] =   0.686188469081757426072759039566E+00;
    xtab[11-1] =   0.846347564651872316865925607099E+00;
    xtab[12-1] =   0.953309846642163911896905464755E+00;
    xtab[13-1] =   1.0E+00;

    weight[1-1] =  0.128205128205128205128205128205E-01;
    weight[2-1] =  0.778016867468189277935889883331E-01;
    weight[3-1] =  0.134981926689608349119914762589E+00;
    weight[4-1] =  0.183646865203550092007494258747E+00;
    weight[5-1] =  0.220767793566110086085534008379E+00;
    weight[6-1] =  0.244015790306676356458578148360E+00;
    weight[7-1] =  0.251930849333446736044138641541E+00;
    weight[8-1] =  0.244015790306676356458578148360E+00;
    weight[9-1] =  0.220767793566110086085534008379E+00;
    weight[10-1] = 0.183646865203550092007494258747E+00;
    weight[11-1] = 0.134981926689608349119914762589E+00;
    weight[12-1] = 0.778016867468189277935889883331E-01;
    weight[13-1] = 0.128205128205128205128205128205E-01;
  }
  else if ( order == 14 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.959935045267260901355100162015E+00;
    xtab[3-1] =  - 0.867801053830347251000220202908E+00;
    xtab[4-1] =  - 0.728868599091326140584672400521E+00;
    xtab[5-1] =  - 0.550639402928647055316622705859E+00;
    xtab[6-1] =  - 0.342724013342712845043903403642E+00;
    xtab[7-1] =  - 0.116331868883703867658776709736E+00;
    xtab[8-1] =    0.116331868883703867658776709736E+00;
    xtab[9-1] =    0.342724013342712845043903403642E+00;
    xtab[10-1] =   0.550639402928647055316622705859E+00;
    xtab[11-1] =   0.728868599091326140584672400521E+00;
    xtab[12-1] =   0.867801053830347251000220202908E+00;
    xtab[13-1] =   0.959935045267260901355100162015E+00;
    xtab[14-1] =   1.0E+00;

    weight[1-1] =  0.109890109890109890109890109890E-01;
    weight[2-1] =  0.668372844976812846340706607461E-01;
    weight[3-1] =  0.116586655898711651540996670655E+00;
    weight[4-1] =  0.160021851762952142412820997988E+00;
    weight[5-1] =  0.194826149373416118640331778376E+00;
    weight[6-1] =  0.219126253009770754871162523954E+00;
    weight[7-1] =  0.231612794468457058889628357293E+00;
    weight[8-1] =  0.231612794468457058889628357293E+00;
    weight[9-1] =  0.219126253009770754871162523954E+00;
    weight[10-1] = 0.194826149373416118640331778376E+00;
    weight[11-1] = 0.160021851762952142412820997988E+00;
    weight[12-1] = 0.116586655898711651540996670655E+00;
    weight[13-1] = 0.668372844976812846340706607461E-01;
    weight[14-1] = 0.109890109890109890109890109890E-01;
  }
  else if ( order == 15 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.965245926503838572795851392070E+00;
    xtab[3-1] =  - 0.885082044222976298825401631482E+00;
    xtab[4-1] =  - 0.763519689951815200704118475976E+00;
    xtab[5-1] =  - 0.606253205469845711123529938637E+00;
    xtab[6-1] =  - 0.420638054713672480921896938739E+00;
    xtab[7-1] =  - 0.215353955363794238225679446273E+00;
    xtab[8-1] =    0.0E+00;
    xtab[9-1] =    0.215353955363794238225679446273E+00;
    xtab[10-1] =   0.420638054713672480921896938739E+00;
    xtab[11-1] =   0.606253205469845711123529938637E+00;
    xtab[12-1] =   0.763519689951815200704118475976E+00;
    xtab[13-1] =   0.885082044222976298825401631482E+00;
    xtab[14-1] =   0.965245926503838572795851392070E+00;
    xtab[15-1] =   1.0E+00;

    weight[1-1] =  0.952380952380952380952380952381E-02;
    weight[2-1] =  0.580298930286012490968805840253E-01;
    weight[3-1] =  0.101660070325718067603666170789E+00;
    weight[4-1] =  0.140511699802428109460446805644E+00;
    weight[5-1] =  0.172789647253600949052077099408E+00;
    weight[6-1] =  0.196987235964613356092500346507E+00;
    weight[7-1] =  0.211973585926820920127430076977E+00;
    weight[8-1] =  0.217048116348815649514950214251E+00;
    weight[9-1] =  0.211973585926820920127430076977E+00;
    weight[10-1] = 0.196987235964613356092500346507E+00;
    weight[11-1] = 0.172789647253600949052077099408E+00;
    weight[12-1] = 0.140511699802428109460446805644E+00;
    weight[13-1] = 0.101660070325718067603666170789E+00;
    weight[14-1] = 0.580298930286012490968805840253E-01;
    weight[15-1] = 0.952380952380952380952380952381E-02;
  }
  else if ( order == 16 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.969568046270217932952242738367E+00;
    xtab[3-1] =  - 0.899200533093472092994628261520E+00;
    xtab[4-1] =  - 0.792008291861815063931088270963E+00;
    xtab[5-1] =  - 0.652388702882493089467883219641E+00;
    xtab[6-1] =  - 0.486059421887137611781890785847E+00;
    xtab[7-1] =  - 0.299830468900763208098353454722E+00;
    xtab[8-1] =  - 0.101326273521949447843033005046E+00;
    xtab[9-1] =    0.101326273521949447843033005046E+00;
    xtab[10-1] =   0.299830468900763208098353454722E+00;
    xtab[11-1] =   0.486059421887137611781890785847E+00;
    xtab[12-1] =   0.652388702882493089467883219641E+00;
    xtab[13-1] =   0.792008291861815063931088270963E+00;
    xtab[14-1] =   0.899200533093472092994628261520E+00;
    xtab[15-1] =   0.969568046270217932952242738367E+00;
    xtab[16-1] =   1.0E+00;

    weight[1-1] =  0.833333333333333333333333333333E-02;
    weight[2-1] =  0.508503610059199054032449195655E-01;
    weight[3-1] =  0.893936973259308009910520801661E-01;
    weight[4-1] =  0.124255382132514098349536332657E+00;
    weight[5-1] =  0.154026980807164280815644940485E+00;
    weight[6-1] =  0.177491913391704125301075669528E+00;
    weight[7-1] =  0.193690023825203584316913598854E+00;
    weight[8-1] =  0.201958308178229871489199125411E+00;
    weight[9-1] =  0.201958308178229871489199125411E+00;
    weight[10-1] = 0.193690023825203584316913598854E+00;
    weight[11-1] = 0.177491913391704125301075669528E+00;
    weight[12-1] = 0.154026980807164280815644940485E+00;
    weight[13-1] = 0.124255382132514098349536332657E+00;
    weight[14-1] = 0.893936973259308009910520801661E-01;
    weight[15-1] = 0.508503610059199054032449195655E-01;
    weight[16-1] = 0.833333333333333333333333333333E-02;
  }
  else if ( order == 17 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.973132176631418314156979501874E+00;
    xtab[3-1] =  - 0.910879995915573595623802506398E+00;
    xtab[4-1] =  - 0.815696251221770307106750553238E+00;
    xtab[5-1] =  - 0.691028980627684705394919357372E+00;
    xtab[6-1] =  - 0.541385399330101539123733407504E+00;
    xtab[7-1] =  - 0.372174433565477041907234680735E+00;
    xtab[8-1] =  - 0.189511973518317388304263014753E+00;
    xtab[9-1] =    0.0E+00;
    xtab[10-1] =   0.189511973518317388304263014753E+00;
    xtab[11-1] =   0.372174433565477041907234680735E+00;
    xtab[12-1] =   0.541385399330101539123733407504E+00;
    xtab[13-1] =   0.691028980627684705394919357372E+00;
    xtab[14-1] =   0.815696251221770307106750553238E+00;
    xtab[15-1] =   0.910879995915573595623802506398E+00;
    xtab[16-1] =   0.973132176631418314156979501874E+00;
    xtab[17-1] =   1.0E+00;

    weight[1-1] =  0.735294117647058823529411764706E-02;
    weight[2-1] =  0.449219405432542096474009546232E-01;
    weight[3-1] =  0.791982705036871191902644299528E-01;
    weight[4-1] =  0.110592909007028161375772705220E+00;
    weight[5-1] =  0.137987746201926559056201574954E+00;
    weight[6-1] =  0.160394661997621539516328365865E+00;
    weight[7-1] =  0.177004253515657870436945745363E+00;
    weight[8-1] =  0.187216339677619235892088482861E+00;
    weight[9-1] =  0.190661874753469433299407247028E+00;
    weight[10-1] = 0.187216339677619235892088482861E+00;
    weight[11-1] = 0.177004253515657870436945745363E+00;
    weight[12-1] = 0.160394661997621539516328365865E+00;
    weight[13-1] = 0.137987746201926559056201574954E+00;
    weight[14-1] = 0.110592909007028161375772705220E+00;
    weight[15-1] = 0.791982705036871191902644299528E-01;
    weight[16-1] = 0.449219405432542096474009546232E-01;
    weight[17-1] = 0.735294117647058823529411764706E-02;
  }
  else if ( order == 18 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.976105557412198542864518924342E+00;
    xtab[3-1] =  - 0.920649185347533873837854625431E+00;
    xtab[4-1] =  - 0.835593535218090213713646362328E+00;
    xtab[5-1] =  - 0.723679329283242681306210365302E+00;
    xtab[6-1] =  - 0.588504834318661761173535893194E+00;
    xtab[7-1] =  - 0.434415036912123975342287136741E+00;
    xtab[8-1] =  - 0.266362652878280984167665332026E+00;
    xtab[9-1] =  - 0.897490934846521110226450100886E-01;
    xtab[10-1] =   0.897490934846521110226450100886E-01;
    xtab[11-1] =   0.266362652878280984167665332026E+00;
    xtab[12-1] =   0.434415036912123975342287136741E+00;
    xtab[13-1] =   0.588504834318661761173535893194E+00;
    xtab[14-1] =   0.723679329283242681306210365302E+00;
    xtab[15-1] =   0.835593535218090213713646362328E+00;
    xtab[16-1] =   0.920649185347533873837854625431E+00;
    xtab[17-1] =   0.976105557412198542864518924342E+00;
    xtab[18-1] =   1.0E+00;

    weight[1-1] =  0.653594771241830065359477124183E-02;
    weight[2-1] =  0.399706288109140661375991764101E-01;
    weight[3-1] =  0.706371668856336649992229601678E-01;
    weight[4-1] =  0.990162717175028023944236053187E-01;
    weight[5-1] =  0.124210533132967100263396358897E+00;
    weight[6-1] =  0.145411961573802267983003210494E+00;
    weight[7-1] =  0.161939517237602489264326706700E+00;
    weight[8-1] =  0.173262109489456226010614403827E+00;
    weight[9-1] =  0.179015863439703082293818806944E+00;
    weight[10-1] = 0.179015863439703082293818806944E+00;
    weight[11-1] = 0.173262109489456226010614403827E+00;
    weight[12-1] = 0.161939517237602489264326706700E+00;
    weight[13-1] = 0.145411961573802267983003210494E+00;
    weight[14-1] = 0.124210533132967100263396358897E+00;
    weight[15-1] = 0.990162717175028023944236053187E-01;
    weight[16-1] = 0.706371668856336649992229601678E-01;
    weight[17-1] = 0.399706288109140661375991764101E-01;
    weight[18-1] = 0.653594771241830065359477124183E-02;
  }
  else if ( order == 19 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.978611766222080095152634063110E+00;
    xtab[3-1] =  - 0.928901528152586243717940258797E+00;
    xtab[4-1] =  - 0.852460577796646093085955970041E+00;
    xtab[5-1] =  - 0.751494202552613014163637489634E+00;
    xtab[6-1] =  - 0.628908137265220497766832306229E+00;
    xtab[7-1] =  - 0.488229285680713502777909637625E+00;
    xtab[8-1] =  - 0.333504847824498610298500103845E+00;
    xtab[9-1] =  - 0.169186023409281571375154153445E+00;
    xtab[10-1] =   0.0E+00;
    xtab[11-1] =   0.169186023409281571375154153445E+00;
    xtab[12-1] =   0.333504847824498610298500103845E+00;
    xtab[13-1] =   0.488229285680713502777909637625E+00;
    xtab[14-1] =   0.628908137265220497766832306229E+00;
    xtab[15-1] =   0.751494202552613014163637489634E+00;
    xtab[16-1] =   0.852460577796646093085955970041E+00;
    xtab[17-1] =   0.928901528152586243717940258797E+00;
    xtab[18-1] =   0.978611766222080095152634063110E+00;
    xtab[19-1] =   1.0E+00;

    weight[1-1] =  0.584795321637426900584795321637E-02;
    weight[2-1] =  0.357933651861764771154255690351E-01;
    weight[3-1] =  0.633818917626297368516956904183E-01;
    weight[4-1] =  0.891317570992070844480087905562E-01;
    weight[5-1] =  0.112315341477305044070910015464E+00;
    weight[6-1] =  0.132267280448750776926046733910E+00;
    weight[7-1] =  0.148413942595938885009680643668E+00;
    weight[8-1] =  0.160290924044061241979910968184E+00;
    weight[9-1] =  0.167556584527142867270137277740E+00;
    weight[10-1] = 0.170001919284827234644672715617E+00;
    weight[11-1] = 0.167556584527142867270137277740E+00;
    weight[12-1] = 0.160290924044061241979910968184E+00;
    weight[13-1] = 0.148413942595938885009680643668E+00;
    weight[14-1] = 0.132267280448750776926046733910E+00;
    weight[15-1] = 0.112315341477305044070910015464E+00;
    weight[16-1] = 0.891317570992070844480087905562E-01;
    weight[17-1] = 0.633818917626297368516956904183E-01;
    weight[18-1] = 0.357933651861764771154255690351E-01;
    weight[19-1] = 0.584795321637426900584795321637E-02;
  }
  else if ( order == 20 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.980743704893914171925446438584E+00;
    xtab[3-1] =  - 0.935934498812665435716181584931E+00;
    xtab[4-1] =  - 0.866877978089950141309847214616E+00;
    xtab[5-1] =  - 0.775368260952055870414317527595E+00;
    xtab[6-1] =  - 0.663776402290311289846403322971E+00;
    xtab[7-1] =  - 0.534992864031886261648135961829E+00;
    xtab[8-1] =  - 0.392353183713909299386474703816E+00;
    xtab[9-1] =  - 0.239551705922986495182401356927E+00;
    xtab[10-1] = - 0.805459372388218379759445181596E-01;
    xtab[11-1] =   0.805459372388218379759445181596E-01;
    xtab[12-1] =   0.239551705922986495182401356927E+00;
    xtab[13-1] =   0.392353183713909299386474703816E+00;
    xtab[14-1] =   0.534992864031886261648135961829E+00;
    xtab[15-1] =   0.663776402290311289846403322971E+00;
    xtab[16-1] =   0.775368260952055870414317527595E+00;
    xtab[17-1] =   0.866877978089950141309847214616E+00;
    xtab[18-1] =   0.935934498812665435716181584931E+00;
    xtab[19-1] =   0.980743704893914171925446438584E+00;
    xtab[20-1] =   1.0E+00;

    weight[1-1] =  0.526315789473684210526315789474E-02;
    weight[2-1] =  0.322371231884889414916050281173E-01;
    weight[3-1] =  0.571818021275668260047536271732E-01;
    weight[4-1] =  0.806317639961196031447768461137E-01;
    weight[5-1] =  0.101991499699450815683781205733E+00;
    weight[6-1] =  0.120709227628674725099429705002E+00;
    weight[7-1] =  0.136300482358724184489780792989E+00;
    weight[8-1] =  0.148361554070916825814713013734E+00;
    weight[9-1] =  0.156580102647475487158169896794E+00;
    weight[10-1] = 0.160743286387845749007726726449E+00;
    weight[11-1] = 0.160743286387845749007726726449E+00;
    weight[12-1] = 0.156580102647475487158169896794E+00;
    weight[13-1] = 0.148361554070916825814713013734E+00;
    weight[14-1] = 0.136300482358724184489780792989E+00;
    weight[15-1] = 0.120709227628674725099429705002E+00;
    weight[16-1] = 0.101991499699450815683781205733E+00;
    weight[17-1] = 0.806317639961196031447768461137E-01;
    weight[18-1] = 0.571818021275668260047536271732E-01;
    weight[19-1] = 0.322371231884889414916050281173E-01;
    weight[20-1] = 0.526315789473684210526315789474E-02;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "LOBATTO_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are between 2 and 20.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

double log_gamma ( double x )

/******************************************************************************/
/*
  Purpose:
  
    LOG_GAMMA calculates the natural logarithm of GAMMA(X).
  
  Discussion:
  
    The method uses Stirling's approximation, and is accurate to about
    12 decimal places.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    27 April 2006
  
  Author:
  
    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
    C version by John Burkardt.
  
  Reference:
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
  Parameters:
  
    Input, double X, the evaluation point.  The routine
    will fail if GAMMA(X) is not positive.  X should be greater than 0.
  
    Output, double LOG_GAMMA, the natural logarithm of the
    gamma function of X.
*/
{
  int i;
  int k;
  int m;
  double p;
  double pi = 3.141592653589793;
  double x2;
  double y;
  double z;

  if ( x < 0.5 )
  {
    m = 1;
    x2 = 1.0 - x;
  }
  else
  {
    m = 0;
    x2 = x;
  }

  k = -1;

  for ( ; ; )
  {
    k = k + 1;

    if ( 6.0 < x2 + ( double ) k )
    {
      break;
    }
  }

  z = x2 + ( double ) k;

  y = ( z - 0.5 ) * log ( z ) - z + 0.9189385332047 +
       ( ( ( ( (
       - 4146.0   / z / z
       + 1820.0 ) / z / z
       - 1287.0 ) / z / z
       + 1716.0 ) / z / z
       - 6006.0 ) / z / z
       + 180180.0 ) / z / 2162160.0;

  if ( 0 < k )
  {
    for ( i = 1; i <= k; i++ )
    {
      y = y - log ( x2 + ( double ) ( k - i ) );
    }

  }

  if ( m != 0 )
  {
    p = pi / sin ( pi * ( 1.0 - x2 ) );

    if ( p <= 0.0 )
    {
      fprintf ( stderr, "\n" );
      fprintf ( stderr, "LOG_GAMMA - Fatal error!\n" );
      exit ( 1 );
    }
    else
    {
      y = log ( p ) - y;
    }
  }

  return y;
}
/******************************************************************************/

void moulton_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    MOULTON_SET sets weights for Adams-Moulton quadrature.
  
  Discussion:
  
    Adams-Moulton quadrature formulas are normally used in solving
    ordinary differential equations, and are not suitable for general
    quadrature computations.  However, an Adams-Moulton formula is
    equivalent to approximating the integral of F(Y(X)) between X(M)
    and X(M+1), using an implicit formula that relies on known values
    of F(Y(X)) at X(M-N+1) through X(M), plus the unknown value at X(M+1).
  
    Suppose the unknown function is denoted by Y(X), with derivative F(Y(X)),
    and that approximate values of the function are known at a series of
    X values, which we write as X(1), X(2), ..., X(M).  We write the value
    Y(X(1)) as Y(1) and so on.
  
    Then the solution of the ODE Y' = F(X,Y) at the next point X(M+1) is
    computed by:
  
      Y(M+1-1] = Y(M) + Integral ( X(M) < X < X(M+1) ) F(Y(X)) dX
             = Y(M) + H * Sum ( 1 <= I <= N ) W(I) * F(Y(M+2-I))
             approximately.
  
    Note that this formula is implicit, since the unknown value Y(M+1)
    appears on the right hand side.  Hence, in ODE applications, this
    equation must be solved via a nonlinear equation solver.  For
    quadrature problems, where the function to be integrated is known
    beforehand, this is not a problem, and the calculation is explicit.
  
    In the documentation that follows, we replace F(Y(X)) by F(X).
  
  
    The Adams-Moulton formulas require equally spaced data.
  
    Here is how the formula is applied in the case with non-unit spacing:
  
      Integral ( A <= X <= A+H ) F(X) dX =
      H * Sum ( 1 <= I <= ORDER ) weight[I) * F ( A - (I-2)*H ),
      approximately.
  
    The integration interval is [ 0, 1 ].
  
    The weight function is w(x-1] = 1.0.
  
    The integral to approximate:
  
      Integral ( 0 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) weight[I) * F ( 2 - I )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    30 April 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Jean Lapidus, John Seinfeld,
    Numerical Solution of Ordinary Differential Equations,
    Academic Press, 1971.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.  ORDER must be
    between 1 and 10 or 12, 14, 16, 18 or 20.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    weight[1) is the weight at X = 1, weight[2) the weight at X = 0,
    and so on.  The weights are rational.  The weights are not
    symmetric, and some weights may be negative.
    They should sum to 1.
*/
{
  double d;
  int i;

  if ( order == 1 )
  {
    weight[1-1] =  1.0;
  }
  else if ( order == 2 )
  {
    d = 2.0;

    weight[1-1] =  1.0 / d;
    weight[2-1] =  1.0 / d;
  }
  else if ( order == 3 )
  {
    d = 12.0;

    weight[1-1] =    5.0 / d;
    weight[2-1] =    8.0 / d;
    weight[3-1] =  - 1.0 / d;
  }
  else if ( order == 4 )
  {
    d = 24.0;

    weight[1-1] =    9.0 / d;
    weight[2-1] =   19.0 / d;
    weight[3-1] =  - 5.0 / d;
    weight[4-1] =    1.0 / d;
  }
  else if ( order == 5 )
  {
    d = 720.0;

    weight[1-1] =    251.0 / d;
    weight[2-1] =    646.0 / d;
    weight[3-1] =  - 264.0 / d;
    weight[4-1] =    106.0 / d;
    weight[5-1] =   - 19.0 / d;
  }
  else if ( order == 6 )
  {
    d = 1440.0;

    weight[1-1] =    475.0 / d;
    weight[2-1] =   1427.0 / d;
    weight[3-1] =  - 798.0 / d;
    weight[4-1] =    482.0 / d;
    weight[5-1] =  - 173.0 / d;
    weight[6-1] =     27.0 / d;
  }
  else if ( order == 7 )
  {
    d = 60480.0;

    weight[1-1] =    19087.0 / d;
    weight[2-1] =    65112.0 / d;
    weight[3-1] =  - 46461.0 / d;
    weight[4-1] =    37504.0 / d;
    weight[5-1] =  - 20211.0 / d;
    weight[6-1] =     6312.0 / d;
    weight[7-1] =    - 863.0 / d;
  }
  else if ( order == 8 )
  {
    d = 120960.0;

    weight[1-1] =    36799.0 / d;
    weight[2-1] =   139849.0 / d;
    weight[3-1] = - 121797.0 / d;
    weight[4-1] =   123133.0 / d;
    weight[5-1] =  - 88547.0 / d;
    weight[6-1] =    41499.0 / d;
    weight[7-1] =  - 11351.0 / d;
    weight[8-1] =     1375.0 / d;
  }
  else if ( order == 9 )
  {
    d = 3628800.0;

    weight[1-1] =   1070017.0 / d;
    weight[2-1] =   4467094.0 / d;
    weight[3-1] = - 4604594.0 / d;
    weight[4-1] =   5595358.0 / d;
    weight[5-1] = - 5033120.0 / d;
    weight[6-1] =   3146338.0 / d;
    weight[7-1] = - 1291214.0 / d;
    weight[8-1] =    312874.0 / d;
    weight[9-1] =   - 33953.0 / d;
  }
  else if ( order == 10 )
  {
    d = 7257600.0;

    weight[1-1] =    2082753.0 / d;
    weight[2-1] =    9449717.0 / d;
    weight[3-1] = - 11271304.0 / d;
    weight[4-1] =   16002320.0 / d;
    weight[5-1] = - 17283646.0 / d;
    weight[6-1] =   13510082.0 / d;
    weight[7-1] =  - 7394032.0 / d;
    weight[8-1] =    2687864.0 / d;
    weight[9-1] =   - 583435.0 / d;
    weight[10-1] =     57281.0 / d;
  }
  else if ( order == 12 )
  {
    d = 958003200.0E+00;

    weight[1-1] =    262747265.0E+00 / d;
    weight[2-1] =   1374799219.0E+00 / d;
    weight[3-1] =  -2092490673.0E+00 / d;
    weight[4-1] =   3828828885.0E+00 / d;
    weight[5-1] =  -5519460582.0E+00 / d;
    weight[6-1] =   6043521486.0E+00 / d;
    weight[7-1] =  -4963166514.0E+00 / d;
    weight[8-1] =   3007739418.0E+00 / d;
    weight[9-1] =  -1305971115.0E+00 / d;
    weight[10-1] =   384709327.0E+00 / d;
    weight[11-1] =   -68928781.0E+00 / d;
    weight[12-1] =     5675265.0E+00 / d;
  }
  else if ( order == 14 )
  {
    d = 5230697472000.0E+00;

    weight[1-1] =    1382741929621.0E+00 / d;
    weight[2-1] =    8153167962181.0E+00 / d;
    weight[3-1] =  -15141235084110.0E+00 / d;
    weight[4-1] =   33928990133618.0E+00 / d;
    weight[5-1] =  -61188680131285.0E+00 / d;
    weight[6-1] =   86180228689563.0E+00 / d;
    weight[7-1] =  -94393338653892.0E+00 / d;
    weight[8-1] =   80101021029180.0E+00 / d;
    weight[9-1] =  -52177910882661.0E+00 / d;
    weight[10-1] =  25620259777835.0E+00 / d;
    weight[11-1] =  -9181635605134.0E+00 / d;
    weight[12-1] =   2268078814386.0E+00 / d;
    weight[13-1] =   -345457086395.0E+00 / d;
    weight[14-1] =     24466579093.0E+00 / d;
  }
  else if ( order == 16 )
  {
    d = 62768369664000.0E+00;

    weight[1-1] =     16088129229375.0E+00 / d;
    weight[2-1] =    105145058757073.0E+00 / d;
    weight[3-1] =   -230992163723849.0E+00 / d;
    weight[4-1] =    612744541065337.0E+00 / d;
    weight[5-1] =  -1326978663058069.0E+00 / d;
    weight[6-1] =   2285168598349733.0E+00 / d;
    weight[7-1] =  -3129453071993581.0E+00 / d;
    weight[8-1] =   3414941728852893.0E+00 / d;
    weight[9-1] =  -2966365730265699.0E+00 / d;
    weight[10-1] =  2039345879546643.0E+00 / d;
    weight[11-1] = -1096355235402331.0E+00 / d;
    weight[12-1] =   451403108933483.0E+00 / d;
    weight[13-1] =  -137515713789319.0E+00 / d;
    weight[14-1] =    29219384284087.0E+00 / d;
    weight[15-1] =    -3867689367599.0E+00 / d;
    weight[16-1] =      240208245823.0E+00 / d;
  }
  else if ( order == 18 )
  {
    d = 64023737057280000.0E+00;

    weight[1-1] =      15980174332775873.0E+00 / d;
    weight[2-1] =     114329243705491117.0E+00 / d;
    weight[3-1] =    -290470969929371220.0E+00 / d;
    weight[4-1] =     890337710266029860.0E+00 / d;
    weight[5-1] =   -2250854333681641520.0E+00 / d;
    weight[6-1] =    4582441343348851896.0E+00 / d;
    weight[7-1] =   -7532171919277411636.0E+00 / d;
    weight[8-1] =   10047287575124288740.0E+00 / d;
    weight[9-1] =  -10910555637627652470.0E+00 / d;
    weight[10-1] =   9644799218032932490.0E+00 / d;
    weight[11-1] =  -6913858539337636636.0E+00 / d;
    weight[12-1] =   3985516155854664396.0E+00 / d;
    weight[13-1] =  -1821304040326216520.0E+00 / d;
    weight[14-1] =    645008976643217360.0E+00 / d;
    weight[15-1] =   -170761422500096220.0E+00 / d;
    weight[16-1] =     31816981024600492.0E+00 / d;
    weight[17-1] =     -3722582669836627.0E+00 / d;
    weight[18-1] =       205804074290625.0E+00 / d;
  }
  else if ( order == 20 )
  {
    d = 102181884343418880000.0E+00;

    weight[1-1] =       24919383499187492303.0E+00 / d;
    weight[2-1] =      193280569173472261637.0E+00 / d;
    weight[3-1] =     -558160720115629395555.0E+00 / d;
    weight[4-1] =     1941395668950986461335.0E+00 / d;
    weight[5-1] =    -5612131802364455926260.0E+00 / d;
    weight[6-1] =    13187185898439270330756.0E+00 / d;
    weight[7-1] =   -25293146116627869170796.0E+00 / d;
    weight[8-1] =    39878419226784442421820.0E+00 / d;
    weight[9-1] =   -51970649453670274135470.0E+00 / d;
    weight[10-1] =   56154678684618739939910.0E+00 / d;
    weight[11-1] =  -50320851025594566473146.0E+00 / d;
    weight[12-1] =   37297227252822858381906.0E+00 / d;
    weight[13-1] =  -22726350407538133839300.0E+00 / d;
    weight[14-1] =   11268210124987992327060.0E+00 / d;
    weight[15-1] =   -4474886658024166985340.0E+00 / d;
    weight[16-1] =    1389665263296211699212.0E+00 / d;
    weight[17-1] =    -325187970422032795497.0E+00 / d;
    weight[18-1] =      53935307402575440285.0E+00 / d;
    weight[19-1] =      -5652892248087175675.0E+00 / d;
    weight[20-1] =        281550972898020815.0E+00 / d;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "MOULTON_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 through 10, 12, 14, 16, 18, or 20.\n" );
    exit ( 1 );
  }
  for ( i = 0; i < order; i++ )
  {
    xtab[i] = ( double ) ( 1 - i );
  }
  return;
}
/******************************************************************************/

void nc_compute ( int order, double x_min, double x_max, double xtab[],
  double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    NC_COMPUTE computes a Newton-Cotes quadrature rule.
  
  Discussion:
  
    For the interval [X_MIN,X_MAX], the Newton-Cotes quadrature rule
    estimates
  
      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
  
    using ORDER abscissas XTAB(I) and a weight vector WEIGHT(I):
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) ).
  
    For the CLOSED rule, the abscissas include the end points.
    For the OPEN rule, the abscissas do not include the end points.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 May 2007
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
  
    Input, double X_MIN, X_MAX, the endpoints of the interval.
  
    Input, double XTAB[ORDER], the abscissas.
  
    Output, double WEIGHT[ORDER], the weights.
*/
{
  int i;
  int j;
  int k;
  double yvala;
  double yvalb;

  double diftab[order];

  for ( i = 1; i <= order; i++ )
  {
/*
  Compute the Lagrange basis polynomial which is 1 at XTAB(I),
  and zero at the other nodes.
*/
    for ( j = 1; j <= order; j++ )
    {
      diftab[j-1] = 0.0;
    }
    diftab[i-1] = 1.0;

    for ( j = 2; j <= order; j++ )
    {
      for ( k = j; k <= order; k++ )
      {
        diftab[order+j-k-1] = ( diftab[order+j-k-2] - diftab[order+j-k-1] )
          / ( xtab[order-k] - xtab[order+j-k-1] );
      }
    }
    for ( j = 1; j <= order-1; j++ )
    {
      for ( k = 1; k <= order - j; k++ )
      {
        diftab[order-k-1] = diftab[order-k-1] - xtab[order-k-j] *
          diftab[order-k];
      }
    }
/*
  Evaluate the antiderivative of the polynomial at the left and
  right endpoints.
*/
    yvala = diftab[order-1] / ( double ) ( order );
    for ( j = order-2; 0 <= j; j-- )
    {
      yvala = yvala * x_min + diftab[j] / ( double ) ( j + 1 );
    }
    yvala = yvala * x_min;

    yvalb = diftab[order-1] / ( double ) ( order );
    for ( j = order-2; 0 <= j; j-- )
    {
      yvalb = yvalb * x_max + diftab[j] / ( double ) ( j + 1 );
    }
    yvalb = yvalb * x_max;

    weight[i-1] = yvalb - yvala;
  }

  return;
}
/******************************************************************************/

void ncc_compute ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    NCC_COMPUTE computes a Newton-Cotes closed quadrature rule.
  
  Discussion:
  
    For the interval [X_MIN,X_MAX], the Newton-Cotes open quadrature rule
    estimates
  
      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
  
    using ORDER equally spaced abscissas XTAB(I) and a weight vector
    WEIGHT(I):
  
      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
  
    For the CLOSED rule, the abscissas include the end points.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 May 2007
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
*/
{
  int i;
  double x_max = 1.0;
  double x_min = -1.0;
/*
  Compute a closed quadrature rule.
*/
  if ( order == 1 )
  {
    xtab[0] = ( x_max + x_min ) / 2.0;
    weight[0] = x_max - x_min;
  }
  else
  {
    for ( i = 0; i < order; i++ )
    {
      xtab[i] = ( ( double ) ( order - i - 1 ) * x_min
                + ( double ) (         i     ) * x_max )
                / ( double ) ( order     - 1 );
    }
    nc_compute ( order, x_min, x_max, xtab, weight );
  }
  return;
}
/******************************************************************************/

void ncc_set ( int order, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    NCC_SET sets abscissas and weights for closed Newton-Cotes quadrature.
  
  Discussion:
  
    The closed Newton-Cotes rules use equally spaced abscissas, and
    hence may be used with tabulated function data.
  
    The rules are called "closed" because they include the endpoints.
    As a favor, we include an order 1 rule, the midpoint rule, even
    though this does not satisfy the requirement that the endpoints
    be included//
  
    The higher order rules involve negative weights.  These can produce
    loss of accuracy due to the subtraction of large, nearly equal quantities.
  
    ORDER = 1 is the midpoint rule (and is not really an NCC rule.)
    ORDER = 2 is the trapezoidal rule.
    ORDER = 3 is Simpson's rule.
    ORDER = 4 is Simpson's 3/8 rule.
    ORDER = 5 is Bode's rule.
  
    The Kopal reference for ORDER = 12 lists
      W(6) = 15494566.0D+00 / 43545600.0D+00
    but this results in a set of coeffients that don't add up to 2.
    The correct value is
      W(6) = 15493566.0D+00 / 43545600.0.
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x) = 1.0.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
  
    In Mathematica, the closed Newton-Cotes weights
    can be computed by:
  
      << NumericalMath`NewtonCotes`
      NewtonCotesWeights [ order, -1, 1, QuadratureType -> Closed ]
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    08 May 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Johnson,
    Quarterly Journal of Mathematics,
    Volume 46, Number 52, 1915.
  
    Zdenek Kopal,
    Numerical Analysis,
    John Wiley, 1955,
    LC: QA297.K6.
  
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Cambridge University Press, 1999,
    ISBN: 0-521-64314-7,
    LC: QA76.95.W65.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3,
    LC: QA47.M315.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 21.
  
    Output, double X[ORDER], the abscissas of the rule.
    The abscissas are uniformly spaced in the interval, and, except
    for the rule of order 1, include -1 and 1.
  
    Output, double W[ORDER], the weights of the rule.
    The weights are symmetric, rational, and should sum to 2.0.
    Some weights may be negative.
*/
{
  if ( order == 1 )
  {
/*
  2
*/
    x[0] = 0.00000000000000000000;

    w[0] = 2.00000000000000000000;
  }
  else if ( order == 2 )
  {
/*
  1
  1
*/
    x[0] = -1.00000000000000000000;
    x[1] =  1.00000000000000000000;

    w[0] = 1.00000000000000000000;
    w[1] = 1.00000000000000000000;
  }
  else if ( order == 3 )
  {
/*
  1 / 3
  4 / 3
  1 / 3
*/
    x[0] = -1.00000000000000000000;
    x[1] =  0.00000000000000000000;
    x[2] =  1.00000000000000000000;

    w[0] = 0.33333333333333333333;
    w[1] = 1.33333333333333333333;
    w[2] = 0.33333333333333333333;
  }
  else if ( order == 4 )
  {
/*
  1 / 4
  3 / 4
  3 / 4
  1 / 4
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.33333333333333333333;
    x[2] =  0.33333333333333333333;
    x[3] =  1.00000000000000000000;

    w[0] = 0.25000000000000000000;
    w[1] = 0.75000000000000000000;
    w[2] = 0.75000000000000000000;
    w[3] = 0.25000000000000000000;
  }
  else if ( order == 5 )
  {
/*
   7 / 45
  32 / 45
  12 / 45
  32 / 45
   7 / 45
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.50000000000000000000;
    x[2] =  0.00000000000000000000;
    x[3] =  0.50000000000000000000;
    x[4] =  1.00000000000000000000;

    w[0] = 0.15555555555555555556;
    w[1] = 0.71111111111111111111;
    w[2] = 0.26666666666666666667;
    w[3] = 0.71111111111111111111;
    w[4] = 0.15555555555555555556;
  }
  else if ( order == 6 )
  {
/*
  19 / 144
  75 / 144
  50 / 144
  50 / 144
  75 / 144
  19 / 144
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.60000000000000000000;
    x[2] = -0.20000000000000000000;
    x[3] =  0.20000000000000000000;
    x[4] =  0.60000000000000000000;
    x[5] =  1.00000000000000000000;

    w[0] = 0.13194444444444444444;
    w[1] = 0.52083333333333333333;
    w[2] = 0.34722222222222222222;
    w[3] = 0.34722222222222222222;
    w[4] = 0.52083333333333333333;
    w[5] = 0.13194444444444444444;
  }
  else if ( order == 7 )
  {
/*
   41 / 420
  216 / 420
   27 / 420
  272 / 420
   27 / 420
  216 / 420
   41 / 420
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.66666666666666666667;
    x[2] = -0.33333333333333333333;
    x[3] =  0.00000000000000000000;
    x[4] =  0.33333333333333333333;
    x[5] =  0.66666666666666666667;
    x[6] =  1.00000000000000000000;

    w[0] = 0.097619047619047619048;
    w[1] = 0.51428571428571428571;
    w[2] = 0.064285714285714285714;
    w[3] = 0.64761904761904761905;
    w[4] = 0.064285714285714285714;
    w[5] = 0.51428571428571428571;
    w[6] = 0.097619047619047619048;
  }
  else if ( order == 8 )
  {
/*
   751 / 8640
  3577 / 8640
  1323 / 8640
  2989 / 8640
  2989 / 8640
  1323 / 8640
  3577 / 8640
   751 / 8640
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.71428571428571428571;
    x[2] = -0.42857142857142857143;
    x[3] = -0.14285714285714285714;
    x[4] =  0.14285714285714285714;
    x[5] =  0.42857142857142857143;
    x[6] =  0.71428571428571428571;
    x[7] =  1.00000000000000000000;

    w[0] = 0.086921296296296296296;
    w[1] = 0.41400462962962962963;
    w[2] = 0.15312500000000000000;
    w[3] = 0.34594907407407407407;
    w[4] = 0.34594907407407407407;
    w[5] = 0.15312500000000000000;
    w[6] = 0.41400462962962962963;
    w[7] = 0.086921296296296296296;
  }
  else if ( order == 9 )
  {
/*
    989 / 14175
   5888 / 14175
   -928 / 14175
  10496 / 14175
  -4540 / 14175
  10496 / 14175
   -928 / 14175
   5888 / 14175
    989 / 14175
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.75000000000000000000;
    x[2] = -0.50000000000000000000;
    x[3] = -0.25000000000000000000;
    x[4] =  0.00000000000000000000;
    x[5] =  0.25000000000000000000;
    x[6] =  0.50000000000000000000;
    x[7] =  0.75000000000000000000;
    x[8] =  1.00000000000000000000;

    w[0] =  0.069770723104056437390;
    w[1] =  0.41537918871252204586;
    w[2] = -0.065467372134038800705;
    w[3] =  0.74045855379188712522;
    w[4] = -0.32028218694885361552;
    w[5] =  0.74045855379188712522;
    w[6] = -0.065467372134038800705;
    w[7] =  0.41537918871252204586;
    w[8] =  0.069770723104056437390;
  }
  else if ( order == 10 )
  {
/*
   2857 / 44800
  15741 / 44800
   1080 / 44800
  19344 / 44800
   5778 / 44800
   5778 / 44800
  19344 / 44800
   1080 / 44800
  15741 / 44800
   2857 / 44800
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.77777777777777777778;
    x[2] = -0.55555555555555555556;
    x[3] = -0.33333333333333333333;
    x[4] = -0.11111111111111111111;
    x[5] =  0.11111111111111111111;
    x[6] =  0.33333333333333333333;
    x[7] =  0.55555555555555555556;
    x[8] =  0.77777777777777777778;
    x[9] =  1.00000000000000000000;

    w[0] =  0.063772321428571428571;
    w[1] =  0.35136160714285714286;
    w[2] =  0.024107142857142857143;
    w[3] =  0.43178571428571428571;
    w[4] =  0.12897321428571428571;
    w[5] =  0.12897321428571428571;
    w[6] =  0.43178571428571428571;
    w[7] =  0.024107142857142857143;
    w[8] =  0.35136160714285714286;
    w[9] =  0.063772321428571428571;
  }
  else if ( order == 11 )
  {
/*
     16067 / 299376
    106300 / 299376
   - 48525 / 299376
    272400 / 299376
  - 260550 / 299376
    427368 / 299376
  - 260550 / 299376
    272400 / 299376
   - 48525 / 299376
    106300 / 299376
     16067 / 299376
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.80000000000000000000;
    x[2] = -0.60000000000000000000;
    x[3] = -0.40000000000000000000;
    x[4] = -0.20000000000000000000;
    x[5] =  0.00000000000000000000;
    x[6] =  0.20000000000000000000;
    x[7] =  0.40000000000000000000;
    x[8] =  0.60000000000000000000;
    x[9] =  0.80000000000000000000;
    x[10] = 1.00000000000000000000;

    w[0] =  0.053668296723852279408;
    w[1] =  0.35507188284966062744;
    w[2] = -0.16208714125380792047;
    w[3] =  0.90989257655924322591;
    w[4] = -0.87031024531024531025;
    w[5] =  1.4275292608625941959;
    w[6] = -0.87031024531024531025;
    w[7] =  0.90989257655924322591;
    w[8] = -0.16208714125380792047;
    w[9] =  0.35507188284966062744;
    w[10] = 0.053668296723852279408;
  }
  else if ( order == 12 )
  {
/*
     2171465 / 43545600
    13486539 / 43545600
   - 3237113 / 43545600
    25226685 / 43545600
   - 9595542 / 43545600
    15493566 / 43545600
    15493566 / 43545600
   - 9595542 / 43545600
    25226685 / 43545600
   - 3237113 / 43545600
    13486539 / 43545600
     2171465 / 43545600
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.81818181818181818182;
    x[2] = -0.63636363636363636364;
    x[3] = -0.45454545454545454545;
    x[4] = -0.27272727272727272727;
    x[5] = -0.090909090909090909091;
    x[6] =  0.090909090909090909091;
    x[7] =  0.27272727272727272727;
    x[8] =  0.45454545454545454545;
    x[9] =  0.63636363636363636364;
    x[10] = 0.81818181818181818182;
    x[11] = 1.00000000000000000000;

    w[0] =   0.049866461823927101705;
    w[1] =   0.30971071704144620811;
    w[2] =  -0.074338463587595532040;
    w[3] =   0.57931650958994708995;
    w[4] =  -0.22035617835097001764;
    w[5] =   0.35580095348324514991;
    w[6] =   0.35580095348324514991;
    w[7] =  -0.22035617835097001764;
    w[8] =   0.57931650958994708995;
    w[9] =  -0.074338463587595532040;
    w[10] =  0.30971071704144620811;
    w[11] =  0.049866461823927101705;
  }
  else if ( order == 13 )
  {
/*
      1364651 / 31531500
      9903168 / 31531500
    - 7587864 / 31531500
     35725120 / 31531500
   - 51491295 / 31531500
     87516288 / 31531500
   - 87797136 / 31531500
     87516288 / 31531500
   - 51491295 / 31531500
     35725120 / 31531500
    - 7587864 / 31531500
      9903168 / 31531500
      1364651 / 31531500
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.83333333333333333333;
    x[2] = -0.66666666666666666667;
    x[3] = -0.50000000000000000000;
    x[4] = -0.33333333333333333333;
    x[5] = -0.16666666666666666667;
    x[6] =  0.00000000000000000000;
    x[7] =  0.16666666666666666667;
    x[8] =  0.33333333333333333333;
    x[9] =  0.50000000000000000000;
    x[10] = 0.66666666666666666667;
    x[11] = 0.83333333333333333333;
    x[12] = 1.00000000000000000000;

    w[0] =   0.043278974993260707546;
    w[1] =   0.31407221350078492936;
    w[2] =  -0.24064392750107035821;
    w[3] =   1.1329977958549387121;
    w[4] =  -1.6330112744398458684;
    w[5] =   2.7755193378050520908;
    w[6] =  -2.7844262404262404262;
    w[7] =   2.7755193378050520908;
    w[8] =  -1.6330112744398458684;
    w[9] =   1.1329977958549387121;
    w[10] = -0.24064392750107035821;
    w[11] =  0.31407221350078492936;
    w[12] =  0.043278974993260707546;
  }
  else if ( order == 14 )
  {
/*
      6137698213 / 150885504000
     42194238652 / 150885504000
   - 23361540993 / 150885504000
    116778274403 / 150885504000
  - 113219777650 / 150885504000
    154424590209 / 150885504000
   - 32067978834 / 150885504000
   - 32067978834 / 150885504000
    154424590209 / 150885504000
  - 113219777650 / 150885504000
    116778274403 / 150885504000
   - 23361540993 / 150885504000
     42194238652 / 150885504000
      6137698213 / 150885504000
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.84615384615384615385;
    x[2] = -0.69230769230769230769;
    x[3] = -0.53846153846153846154;
    x[4] = -0.38461538461538461538;
    x[5] = -0.23076923076923076923;
    x[6] = -0.076923076923076923077;
    x[7] =  0.076923076923076923077;
    x[8] =  0.23076923076923076923;
    x[9] =  0.38461538461538461538;
    x[10] = 0.53846153846153846154;
    x[11] = 0.69230769230769230769;
    x[12] = 0.84615384615384615385;
    x[13] = 1.00000000000000000000;

    w[0] =   0.040669438210247155353;
    w[1] =   0.27975217053157074652;
    w[2] =  -0.15542374057682837445;
    w[3] =   0.77579230848776566369;
    w[4] =  -0.75384763266423526013;
    w[5] =   1.0273523591123107492;
    w[6] =  -0.21429490310083068020;
    w[7] =  -0.21429490310083068020;
    w[8] =   1.0273523591123107492;
    w[9] =  -0.75384763266423526013;
    w[10] =  0.77579230848776566369;
    w[11] = -0.15542374057682837445;
    w[12] =  0.27975217053157074652;
    w[13] =  0.040669438210247155353;
  }
  else if ( order == 15 )
  {
/*
       90241897 / 2501928000
      710986864 / 2501928000
    - 770720657 / 2501928000
     3501442784 / 2501928000
   - 6625093363 / 2501928000
    12630121616 / 2501928000
  - 16802270373 / 2501928000
    19534438464 / 2501928000
  - 16802270373 / 2501928000
    12630121616 / 2501928000
   - 6625093363 / 2501928000
     3501442784 / 2501928000
    - 770720657 / 2501928000
      710986864 / 2501928000
       90241897 / 2501928000
*/
      x[0] = -1.00000000000000000000;
      x[1] = -0.85714285714285714286;
      x[2] = -0.71428571428571428571;
      x[3] = -0.57142857142857142857;
      x[4] = -0.42857142857142857143;
      x[5] = -0.28571428571428571429;
      x[6] = -0.14285714285714285714;
      x[7] =  0.00000000000000000000;
      x[8] =  0.14285714285714285714;
      x[9] =  0.28571428571428571429;
      x[10] = 0.42857142857142857143;
      x[11] = 0.57142857142857142857;
      x[12] = 0.71428571428571428571;
      x[13] = 0.85714285714285714286;
      x[14] = 1.00000000000000000000;
  
      w[0] =   0.036068942431596752584;
      w[1] =   0.28417558938546592868;
      w[2] =  -0.30805069410470645039;
      w[3] =   1.3994978208805369299;
      w[4] =  -2.6479952112930507992;
      w[5] =   5.0481555088715582543;
      w[6] =  -6.7157289790113864188;
      w[7] =   7.8077540456799716059;
      w[8] =  -6.7157289790113864188;
      w[9] =   5.0481555088715582543;
      w[10] = -2.6479952112930507992;
      w[11] =  1.3994978208805369299;
      w[12] = -0.30805069410470645039;
      w[13] =  0.28417558938546592868;
      w[14] =  0.036068942431596752584;
    }
    else if ( order == 16 )
    {
/*
     105930069 / 3099672576
     796661595 / 3099672576
   - 698808195 / 3099672576
    3143332755 / 3099672576
  - 4688522055 / 3099672576
    7385654007 / 3099672576
  - 6000998415 / 3099672576
    3056422815 / 3099672576
    3056422815 / 3099672576
  - 6000998415 / 3099672576
    7385654007 / 3099672576
  - 4688522055 / 3099672576
    3143332755 / 3099672576
   - 698808195 / 3099672576
     796661595 / 3099672576
     105930069 / 3099672576
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.86666666666666666667;
    x[2] = -0.73333333333333333333;
    x[3] = -0.60000000000000000000;
    x[4] = -0.46666666666666666667;
    x[5] = -0.33333333333333333333;
    x[6] = -0.20000000000000000000;
    x[7] = -0.066666666666666666667;
    x[8] =  0.066666666666666666667;
    x[9] =  0.20000000000000000000;
    x[10] = 0.33333333333333333333;
    x[11] = 0.46666666666666666667;
    x[12] = 0.60000000000000000000;
    x[13] = 0.73333333333333333333;
    x[14] = 0.86666666666666666667;
    x[15] = 1.00000000000000000000;

    w[0] =   0.034174599543251887002;
    w[1] =   0.25701475735481036820;
    w[2] =  -0.22544581011901045383;
    w[3] =   1.0140854164204471124;
    w[4] =  -1.5125862296882804695;
    w[5] =   2.3827206990135980091;
    w[6] =  -1.9360104229924960952;
    w[7] =   0.98604699046767964179;
    w[8] =   0.98604699046767964179;
    w[9] =  -1.9360104229924960952;
    w[10] =  2.3827206990135980091;
    w[11] = -1.5125862296882804695;
    w[12] =  1.0140854164204471124;
    w[13] = -0.22544581011901045383;
    w[14] =  0.25701475735481036820;
    w[15] =  0.034174599543251887002;
  }
  else if ( order == 17 )
  {
/*
       15043611773 / 488462349375
      127626606592 / 488462349375
    - 179731134720 / 488462349375
      832211855360 / 488462349375
   - 1929498607520 / 488462349375
     4177588893696 / 488462349375
   - 6806534407936 / 488462349375
     9368875018240 / 488462349375
  - 10234238972220 / 488462349375
     9368875018240 / 488462349375
   - 6806534407936 / 488462349375
     4177588893696 / 488462349375
   - 1929498607520 / 488462349375
      832211855360 / 488462349375
    - 179731134720 / 488462349375
      127626606592 / 488462349375
       15043611773 / 488462349375
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.87500000000000000000;
    x[2] = -0.75000000000000000000;
    x[3] = -0.62500000000000000000;
    x[4] = -0.50000000000000000000;
    x[5] = -0.37500000000000000000;
    x[6] = -0.25000000000000000000;
    x[7] = -0.12500000000000000000;
    x[8] =  0.00000000000000000000;
    x[9] =  0.12500000000000000000;
    x[10] = 0.25000000000000000000;
    x[11] = 0.37500000000000000000;
    x[12] = 0.50000000000000000000;
    x[13] = 0.62500000000000000000;
    x[14] = 0.75000000000000000000;
    x[15] = 0.87500000000000000000;
    x[16] = 1.00000000000000000000;

    w[0]  =   0.030797894233299012495;
    w[1]  =   0.26128238288028031086;
    w[2]  =  -0.36795289329867605622;
    w[3]  =   1.7037379778090086905;
    w[4]  =  -3.9501480717783930427;
    w[5]  =   8.5525299934402953388;
    w[6]  = -13.934614237197880038;
    w[7]  =  19.180342211078732848;
    w[8]  = -20.951950514333334128;
    w[9] =   19.180342211078732848;
    w[10] = -13.934614237197880038;
    w[11] =   8.5525299934402953388;
    w[12] =  -3.9501480717783930427;
    w[13] =   1.7037379778090086905;
    w[14] =  -0.36795289329867605622;
    w[15] =   0.26128238288028031086;
    w[16] =   0.030797894233299012495;
  }
  else if ( order == 18 )
  {
/*
       55294720874657 / 1883051089920000
      450185515446285 / 1883051089920000
    - 542023437008852 / 1883051089920000
     2428636525764260 / 1883051089920000
   - 4768916800123440 / 1883051089920000
     8855416648684984 / 1883051089920000
  - 10905371859796660 / 1883051089920000
    10069615750132836 / 1883051089920000
   - 3759785974054070 / 1883051089920000
   - 3759785974054070 / 1883051089920000
    10069615750132836 / 1883051089920000
  - 10905371859796660 / 1883051089920000
     8855416648684984 / 1883051089920000
   - 4768916800123440 / 1883051089920000
     2428636525764260 / 1883051089920000
    - 542023437008852 / 1883051089920000
      450185515446285 / 1883051089920000
       55294720874657 / 1883051089920000
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.88235294117647058824;
    x[2] = -0.76470588235294117647;
    x[3] = -0.64705882352941176471;
    x[4] = -0.52941176470588235294;
    x[5] = -0.41176470588235294118;
    x[6] = -0.29411764705882352941;
    x[7] = -0.17647058823529411765;
    x[8] = -0.058823529411764705882;
    x[9] =  0.058823529411764705882;
    x[10] = 0.17647058823529411765;
    x[11] = 0.29411764705882352941;
    x[12] = 0.41176470588235294118;
    x[13] = 0.52941176470588235294;
    x[14] = 0.64705882352941176471;
    x[15] = 0.76470588235294117647;
    x[16] = 0.88235294117647058824;
    x[17] = 1.00000000000000000000;

    w[0] =   0.029364429446790078519;
    w[1] =   0.23907238516051669677;
    w[2] =  -0.28784319231183443641;
    w[3] =   1.2897348026109258587;
    w[4] =  -2.5325477495812627261;
    w[5] =   4.7026959045817496499;
    w[6] =  -5.7913308450170443690;
    w[7] =   5.3475000248456540826;
    w[8] =  -1.9966457597354948350;
    w[9] =  -1.9966457597354948350;
    w[10] =  5.3475000248456540826;
    w[11] = -5.7913308450170443690;
    w[12] =  4.7026959045817496499;
    w[13] = -2.5325477495812627261;
    w[14] =  1.2897348026109258587;
    w[15] = -0.28784319231183443641;
    w[16] =  0.23907238516051669677;
    w[17] =  0.029364429446790078519;
  }
  else if ( order == 19 )
  {
/*
       203732352169 / 7604556960000
      1848730221900 / 7604556960000
    - 3212744374395 / 7604556960000
     15529830312096 / 7604556960000
   - 42368630685840 / 7604556960000
    103680563465808 / 7604556960000
  - 198648429867720 / 7604556960000
    319035784479840 / 7604556960000
  - 419127951114198 / 7604556960000
    461327344340680 / 7604556960000
  - 419127951114198 / 7604556960000
    319035784479840 / 7604556960000
  - 198648429867720 / 7604556960000
    103680563465808 / 7604556960000
   - 42368630685840 / 7604556960000
     15529830312096 / 7604556960000
    - 3212744374395 / 7604556960000
      1848730221900 / 7604556960000
       203732352169 / 7604556960000
*/
    x[0] = -1.00000000000000000000;
    x[1] = -0.88888888888888888889;
    x[2] = -0.77777777777777777778;
    x[3] = -0.66666666666666666667;
    x[4] = -0.55555555555555555556;
    x[5] = -0.44444444444444444444;
    x[6] = -0.33333333333333333333;
    x[7] = -0.22222222222222222222;
    x[8] = -0.11111111111111111111;
    x[9] =  0.00000000000000000000;
    x[10] = 0.11111111111111111111;
    x[11] = 0.22222222222222222222;
    x[12] = 0.33333333333333333333;
    x[13] = 0.44444444444444444444;
    x[14] = 0.55555555555555555556;
    x[15] = 0.66666666666666666667;
    x[16] = 0.77777777777777777778;
    x[17] = 0.88888888888888888889;
    x[18] = 1.00000000000000000000;

    w[0] =    0.026790824664820447344;
    w[1] =    0.24310820888374278151;
    w[2] =   -0.42247620621346493274;
    w[3] =    2.0421742376029227612;
    w[4] =   -5.5714791681749728126;
    w[5] =   13.634004454324976218;
    w[6] =  -26.122288374274995239;
    w[7] =   41.953237533490708445;
    w[8] =  -55.115367445968607749;
    w[9] =   60.664591871329740161;
    w[10] = -55.115367445968607749;
    w[11] =  41.953237533490708445;
    w[12] = -26.122288374274995239;
    w[13] =  13.634004454324976218;
    w[14] =  -5.5714791681749728126;
    w[15] =   2.0421742376029227612;
    w[16] =  -0.42247620621346493274;
    w[17] =   0.24310820888374278151;
    w[18] =   0.026790824664820447344;
  }
  else if ( order == 20 )
  {
/*
       69028763155644023 / 2688996956405760000
      603652082270808125 / 2688996956405760000
    - 926840515700222955 / 2688996956405760000
     4301581538450500095 / 2688996956405760000
  - 10343692234243192788 / 2688996956405760000
    22336420328479961316 / 2688996956405760000
  - 35331888421114781580 / 2688996956405760000
    43920768370565135580 / 2688996956405760000
  - 37088370261379851390 / 2688996956405760000
    15148337305921759574 / 2688996956405760000
    15148337305921759574 / 2688996956405760000
  - 37088370261379851390 / 2688996956405760000
    43920768370565135580 / 2688996956405760000
  - 35331888421114781580 / 2688996956405760000
    22336420328479961316 / 2688996956405760000
  - 10343692234243192788 / 2688996956405760000
     4301581538450500095 / 2688996956405760000
    - 926840515700222955 / 2688996956405760000
      603652082270808125 / 2688996956405760000
       69028763155644023 / 2688996956405760000
*/
    x[0] =  -1.00000000000000000000;
    x[1] =  -0.89473684210526315789;
    x[2] =  -0.78947368421052631579;
    x[3] =  -0.68421052631578947368;
    x[4] =  -0.57894736842105263158;
    x[5] =  -0.47368421052631578947;
    x[6] =  -0.36842105263157894737;
    x[7] =  -0.26315789473684210526;
    x[8] =  -0.15789473684210526316;
    x[9] =  -0.052631578947368421053;
    x[10] =  0.052631578947368421053;
    x[11] =  0.15789473684210526316;
    x[12] =  0.26315789473684210526;
    x[13] =  0.36842105263157894737;
    x[14] =  0.47368421052631578947;
    x[15] =  0.57894736842105263158;
    x[16] =  0.68421052631578947368;
    x[17] =  0.78947368421052631579;
    x[18] =  0.89473684210526315789;
    x[19] =  1.00000000000000000000;

    w[0] =    0.025670822345560078100;
    w[1] =    0.22448968595251886556;
    w[2] =   -0.34467890099030890987;
    w[3] =    1.5996974366978074270;
    w[4] =   -3.8466730910952978835;
    w[5] =    8.3065993344729824120;
    w[6] =  -13.139430424771119113;
    w[7] =   16.333513604742678295;
    w[8] =  -13.792641220001198577;
    w[9] =    5.6334527526463774045;
    w[10] =   5.6334527526463774045;
    w[11] = -13.792641220001198577;
    w[12] =  16.333513604742678295;
    w[13] = -13.139430424771119113;
    w[14] =   8.3065993344729824120;
    w[15] =  -3.8466730910952978835;
    w[16] =   1.5996974366978074270;
    w[17] =  -0.34467890099030890987;
    w[18] =   0.22448968595251886556;
    w[19] =   0.025670822345560078100;
  }
  else if ( order == 21 )
  {
    x[0] =  -1.00000000000000000000;
    x[1] =  -0.90000000000000000000;
    x[2] =  -0.80000000000000000000;
    x[3] =  -0.70000000000000000000;
    x[4] =  -0.60000000000000000000;
    x[5] =  -0.50000000000000000000;
    x[6] =  -0.40000000000000000000;
    x[7] =  -0.30000000000000000000;
    x[8] =  -0.20000000000000000000;
    x[9] =  -0.10000000000000000000;
    x[10] =  0.00000000000000000000;
    x[11] =  0.10000000000000000000;
    x[12] =  0.20000000000000000000;
    x[13] =  0.30000000000000000000;
    x[14] =  0.40000000000000000000;
    x[15] =  0.50000000000000000000;
    x[16] =  0.60000000000000000000;
    x[17] =  0.70000000000000000000;
    x[18] =  0.80000000000000000000;
    x[19] =  0.90000000000000000000;
    x[20] =  1.00000000000000000000;

    w[0] =     0.023650546498063206389;
    w[1] =     0.22827543528921394997;
    w[2] =    -0.47295674102285392846;
    w[3] =     2.4123737869637513288;
    w[4] =    -7.5420634534306609355;
    w[5] =    20.673596439879602287;
    w[6] =   -45.417631687959024596;
    w[7] =    83.656114844387109207;
    w[8] =  -128.15055898030800930;
    w[9] =   165.59456694494570344;
    w[10] = -180.01073427048578932;
    w[11] =  165.59456694494570344;
    w[12] = -128.15055898030800930;
    w[13] =   83.656114844387109207;
    w[14] =  -45.417631687959024596;
    w[15] =   20.673596439879602287;
    w[16] =   -7.5420634534306609355;
    w[17] =    2.4123737869637513288;
    w[18] =   -0.47295674102285392846;
    w[19] =    0.22827543528921394997;
    w[20] =    0.023650546498063206389;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "NCC_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 through 21.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void nco_compute ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    NCO_COMPUTE computes a Newton-Cotes open quadrature rule.
  
  Discussion:
  
    For the interval [X_MIN,X_MAX], the Newton-Cotes open quadrature rule
    estimates
  
      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
  
    using ORDER equally spaced abscissas XTAB(I) and a weight vector
    WEIGHT(I):
  
      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) ).
  
    For the OPEN rule, the abscissas do not include the end points.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 May 2007
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
  
    Output, double XTAB[ORDER], the abscissas.
  
    Output, double WEIGHT[ORDER], the weights.
*/
{
  int i;
  double x_max = 1.0;
  double x_min = -1.0;

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = ( ( double ) ( order - i     ) * x_min
              + ( double ) (       + i + 1 ) * x_max )
              / ( double ) ( order     + 1 );
  }

  nc_compute ( order, x_min, x_max, xtab, weight );

  return;
}
/******************************************************************************/

void nco_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    NCO_SET sets abscissas and weights for open Newton-Cotes quadrature.
  
  Discussion:
  
    The open Newton-Cotes rules use equally spaced abscissas, and
    hence may be used with equally spaced data.
  
    The rules are called "open" because they do not include the interval
    endpoints.
  
    Most of the rules involve negative weights.  These can produce loss
    of accuracy due to the subtraction of large, nearly equal quantities.
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x-1] = 1.0.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
    In Mathematica, the open Newton-Cotes weights and abscissas
    can be computed by the commands:
  
      << NumericalMath`NewtonCotes`
      NewtonCotesWeights [ order, -1, 1, QuadratureType -> Open ]
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    02 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Wolfram Media / Cambridge University Press, 1999.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 7, and 9.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are rational, symmetric, and should sum to 2.
    Some weights may be negative.
*/
{
  double d;
  int i;

  if ( order == 1 )
  {
    weight[1-1] = 2.0;
  }
  else if ( order == 2 )
  {
    weight[1-1] = 1.0;
    weight[2-1] = 1.0;
  }
  else if ( order == 3 )
  {
    d = 3.0;

    weight[1-1] =   4.0 / d;
    weight[2-1] = - 2.0 / d;
    weight[3-1] =   4.0 / d;
  }
  else if ( order == 4 )
  {
    d = 12.0;

    weight[1-1] = 11.0 / d;
    weight[2-1] =  1.0 / d;
    weight[3-1] =  1.0 / d;
    weight[4-1] = 11.0 / d;
  }
  else if ( order == 5 )
  {
    d = 10.0;

    weight[1-1] =   11.0 / d;
    weight[2-1] = - 14.0 / d;
    weight[3-1] =   26.0 / d;
    weight[4-1] = - 14.0 / d;
    weight[5-1] =   11.0 / d;
  }
  else if ( order == 6 )
  {
    d = 1440.0;

    weight[1-1] =  1222.0 / d;
    weight[2-1] = - 906.0 / d;
    weight[3-1] =  1124.0 / d;
    weight[4-1] =  1124.0 / d;
    weight[5-1] = - 906.0 / d;
    weight[6-1] =  1222.0 / d;
  }
  else if ( order == 7 )
  {
    d = 945.0;

    weight[1-1] =    920.0 / d;
    weight[2-1] = - 1908.0 / d;
    weight[3-1] =   4392.0 / d;
    weight[4-1] = - 4918.0 / d;
    weight[5-1] =   4392.0 / d;
    weight[6-1] = - 1908.0 / d;
    weight[7-1] =    920.0 / d;
  }
  else if ( order == 9 )
  {
    d = 4536.0;

    weight[1-1] =    4045.0 / d;
    weight[2-1] = - 11690.0 / d;
    weight[3-1] =   33340.0 / d;
    weight[4-1] = - 55070.0 / d;
    weight[5-1] =   67822.0 / d;
    weight[6-1] = - 55070.0 / d;
    weight[7-1] =   33340.0 / d;
    weight[8-1] = - 11690.0 / d;
    weight[9-1] =    4045.0 / d;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "NCO_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 to 7, and 9.\n" );
    exit ( 1 );
  }
/*
  Set the abscissas.
*/
  for ( i = 0; i < order; i++ )
  {
    xtab[i] = ( double ) ( 2 * i - order + 1 )
            / ( double ) ( order + 1 );
  }

  return;
}
/******************************************************************************/

void ncoh_compute ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    NCOH_COMPUTE computes a Newton-Cotes "open half" quadrature rule.
  
  Discussion:
  
    For the interval [X_MIN,X_MAX], the Newton-Cotes "open half" quadrature
    rule estimates
  
      Integral ( X_MIN <= X <= X_MAX ) F(X) dX
  
    using ORDER equally spaced abscissas XTAB(I), each of which is
    the midpoint of one of ORDER equal subintervals,
    and a weight vector WEIGHT(I):
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) ).
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    26 May 2007
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
  
    Output, double XTAB[ORDER], the abscissas.
  
    Output, double WEIGHT[ORDER], the weights.
*/
{
  int i;
  double x_max = 1.0;
  double x_min = -1.0;

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = ( ( double ) ( 2 * order - 2 * i - 1 ) * x_min
              + ( double ) (             2 * i + 1 ) * x_max )
              / ( double ) ( 2 * order             );
  }

  nc_compute ( order, x_min, x_max, xtab, weight );

  return;
}
/******************************************************************************/

void ncoh_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    NCOH_SET sets abscissas and weights for "open half" Newton-Cotes rules.
  
  Discussion:
  
    The open Newton-Cotes rules use equally spaced abscissas, and
    hence may be used with equally spaced data.
  
    The rules are called "open" because the abscissas do not include
    the interval endpoints.
  
    For this uncommon type of open Newton-Cotes rule, the abscissas for
    rule N are found by dividing the interval into N equal subintervals,
    and using the midpoint of each subinterval as the abscissa.
  
    Most of the rules involve negative weights.  These can produce loss
    of accuracy due to the subtraction of large, nearly equal quantities.
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x-1] = 1.0.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
  
    In Mathematica, the "open half" Newton-Cotes weights and abscissas
    can be computed by the commands:
  
      << NumericalMath`NewtonCotes`
      NewtonCotesWeights [ order, -1, 1, QuadratureType -> Open ]
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    05 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Stephen Wolfram,
    The Mathematica Book,
    Fourth Edition,
    Wolfram Media / Cambridge University Press, 1999.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 10.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are rational, symmetric, and should sum to 2.
    Some weights may be negative.
*/
{
  double a;
  double b;
  double d;
  int i;

  if ( order == 1 )
  {
    weight[1-1] = 2.0E+00;
  }
  else if ( order == 2 )
  {
    weight[1-1] = 1.0E+00;
    weight[2-1] = 1.0E+00;
  }
  else if ( order == 3 )
  {
    d = 4.0E+00;

    weight[1-1] =   3.0E+00 / d;
    weight[2-1] =   2.0E+00 / d;
    weight[3-1] =   3.0E+00 / d;
  }
  else if ( order == 4 )
  {
    d = 24.0E+00;

    weight[1-1] = 13.0E+00 / d;
    weight[2-1] = 11.0E+00 / d;
    weight[3-1] = 11.0E+00 / d;
    weight[4-1] = 13.0E+00 / d;
  }
  else if ( order == 5 )
  {
    d = 576.0E+00;

    weight[1-1] =  275.0E+00 / d;
    weight[2-1] =  100.0E+00 / d;
    weight[3-1] =  402.0E+00 / d;
    weight[4-1] =  100.0E+00 / d;
    weight[5-1] =  275.0E+00 / d;
  }
  else if ( order == 6 )
  {
    d = 640.0E+00;

    weight[1-1] =   247.0E+00 / d;
    weight[2-1] =   139.0E+00 / d;
    weight[3-1] =   254.0E+00 / d;
    weight[4-1] =   254.0E+00 / d;
    weight[5-1] =   139.0E+00 / d;
    weight[6-1] =   247.0E+00 / d;
  }
  else if ( order == 7 )
  {
    d = 138240.0E+00;

    weight[1-1] =   49490.0E+00 / d;
    weight[2-1] =    1764.0E+00 / d;
    weight[3-1] =  112014.0E+00 / d;
    weight[4-1] =  -50056.0E+00 / d;
    weight[5-1] =  112014.0E+00 / d;
    weight[6-1] =    1764.0E+00 / d;
    weight[7-1] =   49490.0E+00 / d;
  }
  else if ( order == 8 )
  {
    d = 967680.0E+00;

    weight[1-1] =  295627.0E+00 / d;
    weight[2-1] =   71329.0E+00 / d;
    weight[3-1] =  471771.0E+00 / d;
    weight[4-1] =  128953.0E+00 / d;
    weight[5-1] =  128953.0E+00 / d;
    weight[6-1] =  471771.0E+00 / d;
    weight[7-1] =   71329.0E+00 / d;
    weight[8-1] =  295627.0E+00 / d;
  }
  else if ( order == 9 )
  {
    d = 2867200.0E+00;

    weight[1-1] =    832221.0E+00 / d;
    weight[2-1] =   -260808.0E+00 / d;
    weight[3-1] =   2903148.0E+00 / d;
    weight[4-1] =  -3227256.0E+00 / d;
    weight[5-1] =   5239790.0E+00 / d;
    weight[6-1] =  -3227256.0E+00 / d;
    weight[7-1] =   2903148.0E+00 / d;
    weight[8-1] =   -260808.0E+00 / d;
    weight[9-1] =    832221.0E+00 / d;
  }
  else if ( order == 10 )
  {
    d = 18579456.0E+00;

    weight[1-1] =    4751285.0E+00 / d;
    weight[2-1] =    -492755.0E+00 / d;
    weight[3-1] =   12269956.0E+00 / d;
    weight[4-1] =   -6274220.0E+00 / d;
    weight[5-1] =    8325190.0E+00 / d;
    weight[6-1] =    8325190.0E+00 / d;
    weight[7-1] =   -6274220.0E+00 / d;
    weight[8-1] =   12269956.0E+00 / d;
    weight[9-1] =    -492755.0E+00 / d;
    weight[10-1] =   4751285.0E+00 / d;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "NCOH_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 to 10.\n" );
    exit ( 1 );
  }
/*
  Set the abscissas.
*/
  a = -1.0E+00;
  b = +1.0E+00;

  for ( i = 0; i < order; i++ )
  {
    xtab[i] = ( ( double ) ( 2 * order - 2 * i - 1 ) * a
              + ( double ) (             2 * i + 1 ) * b )
              / ( double ) ( 2 * order                   );
  }

  return;
}
/******************************************************************************/

void patterson_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    PATTERSON_SET sets abscissas and weights for Gauss-Patterson quadrature.
  
  Discussion:
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x-1] = 1.0.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
  
    The zeroth rule, of order 1, is the standard Gauss-Legendre rule.
  
    The first rule, of order 3, is the standard Gauss-Legendre rule.
  
    The second rule, of order 7, includes the abscissas of the previous
    rule.
  
    Each subsequent rule is nested in a similar way.  Rules are available
    of orders 1, 3, 7, 15, 31, 63 and 127.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    14 April 2007
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Prem Kythe, Michael Schaeferkotter,
    Handbook of Computational Methods for Integration,
    Chapman and Hall, 2004,
    ISBN: 1-58488-428-2,
    LC: QA299.3.K98.
  
    Thomas Patterson,
    The Optimal Addition of Points to Quadrature Formulae,
    Mathematics of Computation,
    Volume 22, Number 104, October 1968, pages 847-856.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be 1, 3, 7, 15, 31, 63 or 127.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are positive, symmetric and should sum to 2.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] =   0.0;

    weight[1-1] = 2.0;
  }
  else if ( order == 3 )
  {
    xtab[1-1] = -0.77459666924148337704;
    xtab[2-1] =  0.0;
    xtab[3-1] =  0.77459666924148337704;

    weight[1-1] = 0.555555555555555555556;
    weight[2-1] = 0.888888888888888888889;
    weight[3-1] = 0.555555555555555555556;
  }
  else if ( order == 7 )
  {
    xtab[1-1] = -0.96049126870802028342;
    xtab[2-1] = -0.77459666924148337704;
    xtab[3-1] = -0.43424374934680255800;
    xtab[4-1] =  0.0;
    xtab[5-1] =  0.43424374934680255800;
    xtab[6-1] =  0.77459666924148337704;
    xtab[7-1] =  0.96049126870802028342;

    weight[1-1] = 0.104656226026467265194;
    weight[2-1] = 0.268488089868333440729;
    weight[3-1] = 0.401397414775962222905;
    weight[4-1] = 0.450916538658474142345;
    weight[5-1] = 0.401397414775962222905;
    weight[6-1] = 0.268488089868333440729;
    weight[7-1] = 0.104656226026467265194;
  }
  else if ( order == 15 )
  {
    xtab[ 1-1] = -0.99383196321275502221;
    xtab[ 2-1] = -0.96049126870802028342;
    xtab[ 3-1] = -0.88845923287225699889;
    xtab[ 4-1] = -0.77459666924148337704;
    xtab[ 5-1] = -0.62110294673722640294;
    xtab[ 6-1] = -0.43424374934680255800;
    xtab[ 7-1] = -0.22338668642896688163;
    xtab[ 8-1] =  0.0;
    xtab[ 9-1] =  0.22338668642896688163;
    xtab[10-1] =  0.43424374934680255800;
    xtab[11-1] =  0.62110294673722640294;
    xtab[12-1] =  0.77459666924148337704;
    xtab[13-1] =  0.88845923287225699889;
    xtab[14-1] =  0.96049126870802028342;
    xtab[15-1] =  0.99383196321275502221;

    weight[ 1-1] = 0.0170017196299402603390;
    weight[ 2-1] = 0.0516032829970797396969;
    weight[ 3-1] = 0.0929271953151245376859;
    weight[ 4-1] = 0.134415255243784220360;
    weight[ 5-1] = 0.171511909136391380787;
    weight[ 6-1] = 0.200628529376989021034;
    weight[ 7-1] = 0.219156858401587496404;
    weight[ 8-1] = 0.225510499798206687386;
    weight[ 9-1] = 0.219156858401587496404;
    weight[10-1] = 0.200628529376989021034;
    weight[11-1] = 0.171511909136391380787;
    weight[12-1] = 0.134415255243784220360;
    weight[13-1] = 0.0929271953151245376859;
    weight[14-1] = 0.0516032829970797396969;
    weight[15-1] = 0.0170017196299402603390;
  }
  else if ( order == 31 )
  {
    xtab[ 1-1] = -0.99909812496766759766;
    xtab[ 2-1] = -0.99383196321275502221;
    xtab[ 3-1] = -0.98153114955374010687;
    xtab[ 4-1] = -0.96049126870802028342;
    xtab[ 5-1] = -0.92965485742974005667;
    xtab[ 6-1] = -0.88845923287225699889;
    xtab[ 7-1] = -0.83672593816886873550;
    xtab[ 8-1] = -0.77459666924148337704;
    xtab[ 9-1] = -0.70249620649152707861;
    xtab[10-1] = -0.62110294673722640294;
    xtab[11-1] = -0.53131974364437562397;
    xtab[12-1] = -0.43424374934680255800;
    xtab[13-1] = -0.33113539325797683309;
    xtab[14-1] = -0.22338668642896688163;
    xtab[15-1] = -0.11248894313318662575;
    xtab[16-1] =  0.0;
    xtab[17-1] =  0.11248894313318662575;
    xtab[18-1] =  0.22338668642896688163;
    xtab[19-1] =  0.33113539325797683309;
    xtab[20-1] =  0.43424374934680255800;
    xtab[21-1] =  0.53131974364437562397;
    xtab[22-1] =  0.62110294673722640294;
    xtab[23-1] =  0.70249620649152707861;
    xtab[24-1] =  0.77459666924148337704;
    xtab[25-1] =  0.83672593816886873550;
    xtab[26-1] =  0.88845923287225699889;
    xtab[27-1] =  0.92965485742974005667;
    xtab[28-1] =  0.96049126870802028342;
    xtab[29-1] =  0.98153114955374010687;
    xtab[30-1] =  0.99383196321275502221;
    xtab[31-1] =  0.99909812496766759766;

    weight[ 1-1] = 0.00254478079156187441540;
    weight[ 2-1] = 0.00843456573932110624631;
    weight[ 3-1] = 0.0164460498543878109338;
    weight[ 4-1] = 0.0258075980961766535646;
    weight[ 5-1] = 0.0359571033071293220968;
    weight[ 6-1] = 0.0464628932617579865414;
    weight[ 7-1] = 0.0569795094941233574122;
    weight[ 8-1] = 0.0672077542959907035404;
    weight[ 9-1] = 0.0768796204990035310427;
    weight[10-1] = 0.0857559200499903511542;
    weight[11-1] = 0.0936271099812644736167;
    weight[12-1] = 0.100314278611795578771;
    weight[13-1] = 0.105669893580234809744;
    weight[14-1] = 0.109578421055924638237;
    weight[15-1] = 0.111956873020953456880;
    weight[16-1] = 0.112755256720768691607;
    weight[17-1] = 0.111956873020953456880;
    weight[18-1] = 0.109578421055924638237;
    weight[19-1] = 0.105669893580234809744;
    weight[20-1] = 0.100314278611795578771;
    weight[21-1] = 0.0936271099812644736167;
    weight[22-1] = 0.0857559200499903511542;
    weight[23-1] = 0.0768796204990035310427;
    weight[24-1] = 0.0672077542959907035404;
    weight[25-1] = 0.0569795094941233574122;
    weight[26-1] = 0.0464628932617579865414;
    weight[27-1] = 0.0359571033071293220968;
    weight[28-1] = 0.0258075980961766535646;
    weight[29-1] = 0.0164460498543878109338;
    weight[30-1] = 0.00843456573932110624631;
    weight[31-1] = 0.00254478079156187441540;
  }
  else if ( order == 63 )
  {
    xtab[ 1-1] = -0.99987288812035761194;
    xtab[ 2-1] = -0.99909812496766759766;
    xtab[ 3-1] = -0.99720625937222195908;
    xtab[ 4-1] = -0.99383196321275502221;
    xtab[ 5-1] = -0.98868475754742947994;
    xtab[ 6-1] = -0.98153114955374010687;
    xtab[ 7-1] = -0.97218287474858179658;
    xtab[ 8-1] = -0.96049126870802028342;
    xtab[ 9-1] = -0.94634285837340290515;
    xtab[10-1] = -0.92965485742974005667;
    xtab[11-1] = -0.91037115695700429250;
    xtab[12-1] = -0.88845923287225699889;
    xtab[13-1] = -0.86390793819369047715;
    xtab[14-1] = -0.83672593816886873550;
    xtab[15-1] = -0.80694053195021761186;
    xtab[16-1] = -0.77459666924148337704;
    xtab[17-1] = -0.73975604435269475868;
    xtab[18-1] = -0.70249620649152707861;
    xtab[19-1] = -0.66290966002478059546;
    xtab[20-1] = -0.62110294673722640294;
    xtab[21-1] = -0.57719571005204581484;
    xtab[22-1] = -0.53131974364437562397;
    xtab[23-1] = -0.48361802694584102756;
    xtab[24-1] = -0.43424374934680255800;
    xtab[25-1] = -0.38335932419873034692;
    xtab[26-1] = -0.33113539325797683309;
    xtab[27-1] = -0.27774982202182431507;
    xtab[28-1] = -0.22338668642896688163;
    xtab[29-1] = -0.16823525155220746498;
    xtab[30-1] = -0.11248894313318662575;
    xtab[31-1] = -0.056344313046592789972;
    xtab[32-1] =  0.0;
    xtab[33-1] =  0.056344313046592789972;
    xtab[34-1] =  0.11248894313318662575;
    xtab[35-1] =  0.16823525155220746498;
    xtab[36-1] =  0.22338668642896688163;
    xtab[37-1] =  0.27774982202182431507;
    xtab[38-1] =  0.33113539325797683309;
    xtab[39-1] =  0.38335932419873034692;
    xtab[40-1] =  0.43424374934680255800;
    xtab[41-1] =  0.48361802694584102756;
    xtab[42-1] =  0.53131974364437562397;
    xtab[43-1] =  0.57719571005204581484;
    xtab[44-1] =  0.62110294673722640294;
    xtab[45-1] =  0.66290966002478059546;
    xtab[46-1] =  0.70249620649152707861;
    xtab[47-1] =  0.73975604435269475868;
    xtab[48-1] =  0.77459666924148337704;
    xtab[49-1] =  0.80694053195021761186;
    xtab[50-1] =  0.83672593816886873550;
    xtab[51-1] =  0.86390793819369047715;
    xtab[52-1] =  0.88845923287225699889;
    xtab[53-1] =  0.91037115695700429250;
    xtab[54-1] =  0.92965485742974005667;
    xtab[55-1] =  0.94634285837340290515;
    xtab[56-1] =  0.96049126870802028342;
    xtab[57-1] =  0.97218287474858179658;
    xtab[58-1] =  0.98153114955374010687;
    xtab[59-1] =  0.98868475754742947994;
    xtab[60-1] =  0.99383196321275502221;
    xtab[61-1] =  0.99720625937222195908;
    xtab[62-1] =  0.99909812496766759766;
    xtab[63-1] =  0.99987288812035761194;

    weight[ 1-1] = 0.000363221481845530659694;
    weight[ 2-1] = 0.00126515655623006801137;
    weight[ 3-1] = 0.00257904979468568827243;
    weight[ 4-1] = 0.00421763044155885483908;
    weight[ 5-1] = 0.00611550682211724633968;
    weight[ 6-1] = 0.00822300795723592966926;
    weight[ 7-1] = 0.0104982469096213218983;
    weight[ 8-1] = 0.0129038001003512656260;
    weight[ 9-1] = 0.0154067504665594978021;
    weight[10-1] = 0.0179785515681282703329;
    weight[11-1] = 0.0205942339159127111492;
    weight[12-1] = 0.0232314466399102694433;
    weight[13-1] = 0.0258696793272147469108;
    weight[14-1] = 0.0284897547458335486125;
    weight[15-1] = 0.0310735511116879648799;
    weight[16-1] = 0.0336038771482077305417;
    weight[17-1] = 0.0360644327807825726401;
    weight[18-1] = 0.0384398102494555320386;
    weight[19-1] = 0.0407155101169443189339;
    weight[20-1] = 0.0428779600250077344929;
    weight[21-1] = 0.0449145316536321974143;
    weight[22-1] = 0.0468135549906280124026;
    weight[23-1] = 0.0485643304066731987159;
    weight[24-1] = 0.0501571393058995374137;
    weight[25-1] = 0.0515832539520484587768;
    weight[26-1] = 0.0528349467901165198621;
    weight[27-1] = 0.0539054993352660639269;
    weight[28-1] = 0.0547892105279628650322;
    weight[29-1] = 0.0554814043565593639878;
    weight[30-1] = 0.0559784365104763194076;
    weight[31-1] = 0.0562776998312543012726;
    weight[32-1] = 0.0563776283603847173877;
    weight[33-1] = 0.0562776998312543012726;
    weight[34-1] = 0.0559784365104763194076;
    weight[35-1] = 0.0554814043565593639878;
    weight[36-1] = 0.0547892105279628650322;
    weight[37-1] = 0.0539054993352660639269;
    weight[38-1] = 0.0528349467901165198621;
    weight[39-1] = 0.0515832539520484587768;
    weight[40-1] = 0.0501571393058995374137;
    weight[41-1] = 0.0485643304066731987159;
    weight[42-1] = 0.0468135549906280124026;
    weight[43-1] = 0.0449145316536321974143;
    weight[44-1] = 0.0428779600250077344929;
    weight[45-1] = 0.0407155101169443189339;
    weight[46-1] = 0.0384398102494555320386;
    weight[47-1] = 0.0360644327807825726401;
    weight[48-1] = 0.0336038771482077305417;
    weight[49-1] = 0.0310735511116879648799;
    weight[50-1] = 0.0284897547458335486125;
    weight[51-1] = 0.0258696793272147469108;
    weight[52-1] = 0.0232314466399102694433;
    weight[53-1] = 0.0205942339159127111492;
    weight[54-1] = 0.0179785515681282703329;
    weight[55-1] = 0.0154067504665594978021;
    weight[56-1] = 0.0129038001003512656260;
    weight[57-1] = 0.0104982469096213218983;
    weight[58-1] = 0.00822300795723592966926;
    weight[59-1] = 0.00611550682211724633968;
    weight[60-1] = 0.00421763044155885483908;
    weight[61-1] = 0.00257904979468568827243;
    weight[62-1] = 0.00126515655623006801137;
    weight[63-1] = 0.000363221481845530659694;
  }
  else if ( order == 127 )
  {
    xtab[  1-1] = -0.99998243035489159858;
    xtab[  2-1] = -0.99987288812035761194;
    xtab[  3-1] = -0.99959879967191068325;
    xtab[  4-1] = -0.99909812496766759766;
    xtab[  5-1] = -0.99831663531840739253;
    xtab[  6-1] = -0.99720625937222195908;
    xtab[  7-1] = -0.99572410469840718851;
    xtab[  8-1] = -0.99383196321275502221;
    xtab[  9-1] = -0.99149572117810613240;
    xtab[ 10-1] = -0.98868475754742947994;
    xtab[ 11-1] = -0.98537149959852037111;
    xtab[ 12-1] = -0.98153114955374010687;
    xtab[ 13-1] = -0.97714151463970571416;
    xtab[ 14-1] = -0.97218287474858179658;
    xtab[ 15-1] = -0.96663785155841656709;
    xtab[ 16-1] = -0.96049126870802028342;
    xtab[ 17-1] = -0.95373000642576113641;
    xtab[ 18-1] = -0.94634285837340290515;
    xtab[ 19-1] = -0.93832039777959288365;
    xtab[ 20-1] = -0.92965485742974005667;
    xtab[ 21-1] = -0.92034002547001242073;
    xtab[ 22-1] = -0.91037115695700429250;
    xtab[ 23-1] = -0.89974489977694003664;
    xtab[ 24-1] = -0.88845923287225699889;
    xtab[ 25-1] = -0.87651341448470526974;
    xtab[ 26-1] = -0.86390793819369047715;
    xtab[ 27-1] = -0.85064449476835027976;
    xtab[ 28-1] = -0.83672593816886873550;
    xtab[ 29-1] = -0.82215625436498040737;
    xtab[ 30-1] = -0.80694053195021761186;
    xtab[ 31-1] = -0.79108493379984836143;
    xtab[ 32-1] = -0.77459666924148337704;
    xtab[ 33-1] = -0.75748396638051363793;
    xtab[ 34-1] = -0.73975604435269475868;
    xtab[ 35-1] = -0.72142308537009891548;
    xtab[ 36-1] = -0.70249620649152707861;
    xtab[ 37-1] = -0.68298743109107922809;
    xtab[ 38-1] = -0.66290966002478059546;
    xtab[ 39-1] = -0.64227664250975951377;
    xtab[ 40-1] = -0.62110294673722640294;
    xtab[ 41-1] = -0.59940393024224289297;
    xtab[ 42-1] = -0.57719571005204581484;
    xtab[ 43-1] = -0.55449513263193254887;
    xtab[ 44-1] = -0.53131974364437562397;
    xtab[ 45-1] = -0.50768775753371660215;
    xtab[ 46-1] = -0.48361802694584102756;
    xtab[ 47-1] = -0.45913001198983233287;
    xtab[ 48-1] = -0.43424374934680255800;
    xtab[ 49-1] = -0.40897982122988867241;
    xtab[ 50-1] = -0.38335932419873034692;
    xtab[ 51-1] = -0.35740383783153215238;
    xtab[ 52-1] = -0.33113539325797683309;
    xtab[ 53-1] = -0.30457644155671404334;
    xtab[ 54-1] = -0.27774982202182431507;
    xtab[ 55-1] = -0.25067873030348317661;
    xtab[ 56-1] = -0.22338668642896688163;
    xtab[ 57-1] = -0.19589750271110015392;
    xtab[ 58-1] = -0.16823525155220746498;
    xtab[ 59-1] = -0.14042423315256017459;
    xtab[ 60-1] = -0.11248894313318662575;
    xtab[ 61-1] = -0.084454040083710883710;
    xtab[ 62-1] = -0.056344313046592789972;
    xtab[ 63-1] = -0.028184648949745694339;
    xtab[ 64-1] =  0.0;
    xtab[ 65-1] =  0.028184648949745694339;
    xtab[ 66-1] =  0.056344313046592789972;
    xtab[ 67-1] =  0.084454040083710883710;
    xtab[ 68-1] =  0.11248894313318662575;
    xtab[ 69-1] =  0.14042423315256017459;
    xtab[ 70-1] =  0.16823525155220746498;
    xtab[ 71-1] =  0.19589750271110015392;
    xtab[ 72-1] =  0.22338668642896688163;
    xtab[ 73-1] =  0.25067873030348317661;
    xtab[ 74-1] =  0.27774982202182431507;
    xtab[ 75-1] =  0.30457644155671404334;
    xtab[ 76-1] =  0.33113539325797683309;
    xtab[ 77-1] =  0.35740383783153215238;
    xtab[ 78-1] =  0.38335932419873034692;
    xtab[ 79-1] =  0.40897982122988867241;
    xtab[ 80-1] =  0.43424374934680255800;
    xtab[ 81-1] =  0.45913001198983233287;
    xtab[ 82-1] =  0.48361802694584102756;
    xtab[ 83-1] =  0.50768775753371660215;
    xtab[ 84-1] =  0.53131974364437562397;
    xtab[ 85-1] =  0.55449513263193254887;
    xtab[ 86-1] =  0.57719571005204581484;
    xtab[ 87-1] =  0.59940393024224289297;
    xtab[ 88-1] =  0.62110294673722640294;
    xtab[ 89-1] =  0.64227664250975951377;
    xtab[ 90-1] =  0.66290966002478059546;
    xtab[ 91-1] =  0.68298743109107922809;
    xtab[ 92-1] =  0.70249620649152707861;
    xtab[ 93-1] =  0.72142308537009891548;
    xtab[ 94-1] =  0.73975604435269475868;
    xtab[ 95-1] =  0.75748396638051363793;
    xtab[ 96-1] =  0.77459666924148337704;
    xtab[ 97-1] =  0.79108493379984836143;
    xtab[ 98-1] =  0.80694053195021761186;
    xtab[ 99-1] =  0.82215625436498040737;
    xtab[100-1] =  0.83672593816886873550;
    xtab[101-1] =  0.85064449476835027976;
    xtab[102-1] =  0.86390793819369047715;
    xtab[103-1] =  0.87651341448470526974;
    xtab[104-1] =  0.88845923287225699889;
    xtab[105-1] =  0.89974489977694003664;
    xtab[106-1] =  0.91037115695700429250;
    xtab[107-1] =  0.92034002547001242073;
    xtab[108-1] =  0.92965485742974005667;
    xtab[109-1] =  0.93832039777959288365;
    xtab[110-1] =  0.94634285837340290515;
    xtab[111-1] =  0.95373000642576113641;
    xtab[112-1] =  0.96049126870802028342;
    xtab[113-1] =  0.96663785155841656709;
    xtab[114-1] =  0.97218287474858179658;
    xtab[115-1] =  0.97714151463970571416;
    xtab[116-1] =  0.98153114955374010687;
    xtab[117-1] =  0.98537149959852037111;
    xtab[118-1] =  0.98868475754742947994;
    xtab[119-1] =  0.99149572117810613240;
    xtab[120-1] =  0.99383196321275502221;
    xtab[121-1] =  0.99572410469840718851;
    xtab[122-1] =  0.99720625937222195908;
    xtab[123-1] =  0.99831663531840739253;
    xtab[124-1] =  0.99909812496766759766;
    xtab[125-1] =  0.99959879967191068325;
    xtab[126-1] =  0.99987288812035761194;
    xtab[127-1] =  0.99998243035489159858;

    weight[  1-1] = 0.0000505360952078625176247;
    weight[  2-1] = 0.000180739564445388357820;
    weight[  3-1] = 0.000377746646326984660274;
    weight[  4-1] = 0.000632607319362633544219;
    weight[  5-1] = 0.000938369848542381500794;
    weight[  6-1] = 0.00128952408261041739210;
    weight[  7-1] = 0.00168114286542146990631;
    weight[  8-1] = 0.00210881524572663287933;
    weight[  9-1] = 0.00256876494379402037313;
    weight[ 10-1] = 0.00305775341017553113613;
    weight[ 11-1] = 0.00357289278351729964938;
    weight[ 12-1] = 0.00411150397865469304717;
    weight[ 13-1] = 0.00467105037211432174741;
    weight[ 14-1] = 0.00524912345480885912513;
    weight[ 15-1] = 0.00584344987583563950756;
    weight[ 16-1] = 0.00645190005017573692280;
    weight[ 17-1] = 0.00707248999543355546805;
    weight[ 18-1] = 0.00770337523327974184817;
    weight[ 19-1] = 0.00834283875396815770558;
    weight[ 20-1] = 0.00898927578406413572328;
    weight[ 21-1] = 0.00964117772970253669530;
    weight[ 22-1] = 0.0102971169579563555237;
    weight[ 23-1] = 0.0109557333878379016480;
    weight[ 24-1] = 0.0116157233199551347270;
    weight[ 25-1] = 0.0122758305600827700870;
    weight[ 26-1] = 0.0129348396636073734547;
    weight[ 27-1] = 0.0135915710097655467896;
    weight[ 28-1] = 0.0142448773729167743063;
    weight[ 29-1] = 0.0148936416648151820348;
    weight[ 30-1] = 0.0155367755558439824399;
    weight[ 31-1] = 0.0161732187295777199419;
    weight[ 32-1] = 0.0168019385741038652709;
    weight[ 33-1] = 0.0174219301594641737472;
    weight[ 34-1] = 0.0180322163903912863201;
    weight[ 35-1] = 0.0186318482561387901863;
    weight[ 36-1] = 0.0192199051247277660193;
    weight[ 37-1] = 0.0197954950480974994880;
    weight[ 38-1] = 0.0203577550584721594669;
    weight[ 39-1] = 0.0209058514458120238522;
    weight[ 40-1] = 0.0214389800125038672465;
    weight[ 41-1] = 0.0219563663053178249393;
    weight[ 42-1] = 0.0224572658268160987071;
    weight[ 43-1] = 0.0229409642293877487608;
    weight[ 44-1] = 0.0234067774953140062013;
    weight[ 45-1] = 0.0238540521060385400804;
    weight[ 46-1] = 0.0242821652033365993580;
    weight[ 47-1] = 0.0246905247444876769091;
    weight[ 48-1] = 0.0250785696529497687068;
    weight[ 49-1] = 0.0254457699654647658126;
    weight[ 50-1] = 0.0257916269760242293884;
    weight[ 51-1] = 0.0261156733767060976805;
    weight[ 52-1] = 0.0264174733950582599310;
    weight[ 53-1] = 0.0266966229274503599062;
    weight[ 54-1] = 0.0269527496676330319634;
    weight[ 55-1] = 0.0271855132296247918192;
    weight[ 56-1] = 0.0273946052639814325161;
    weight[ 57-1] = 0.0275797495664818730349;
    weight[ 58-1] = 0.0277407021782796819939;
    weight[ 59-1] = 0.0278772514766137016085;
    weight[ 60-1] = 0.0279892182552381597038;
    weight[ 61-1] = 0.0280764557938172466068;
    weight[ 62-1] = 0.0281388499156271506363;
    weight[ 63-1] = 0.0281763190330166021307;
    weight[ 64-1] = 0.0281888141801923586938;
    weight[ 65-1] = 0.0281763190330166021307;
    weight[ 66-1] = 0.0281388499156271506363;
    weight[ 67-1] = 0.0280764557938172466068;
    weight[ 68-1] = 0.0279892182552381597038;
    weight[ 69-1] = 0.0278772514766137016085;
    weight[ 70-1] = 0.0277407021782796819939;
    weight[ 71-1] = 0.0275797495664818730349;
    weight[ 72-1] = 0.0273946052639814325161;
    weight[ 73-1] = 0.0271855132296247918192;
    weight[ 74-1] = 0.0269527496676330319634;
    weight[ 75-1] = 0.0266966229274503599062;
    weight[ 76-1] = 0.0264174733950582599310;
    weight[ 77-1] = 0.0261156733767060976805;
    weight[ 78-1] = 0.0257916269760242293884;
    weight[ 79-1] = 0.0254457699654647658126;
    weight[ 80-1] = 0.0250785696529497687068;
    weight[ 81-1] = 0.0246905247444876769091;
    weight[ 82-1] = 0.0242821652033365993580;
    weight[ 83-1] = 0.0238540521060385400804;
    weight[ 84-1] = 0.0234067774953140062013;
    weight[ 85-1] = 0.0229409642293877487608;
    weight[ 86-1] = 0.0224572658268160987071;
    weight[ 87-1] = 0.0219563663053178249393;
    weight[ 88-1] = 0.0214389800125038672465;
    weight[ 89-1] = 0.0209058514458120238522;
    weight[ 90-1] = 0.0203577550584721594669;
    weight[ 91-1] = 0.0197954950480974994880;
    weight[ 92-1] = 0.0192199051247277660193;
    weight[ 93-1] = 0.0186318482561387901863;
    weight[ 94-1] = 0.0180322163903912863201;
    weight[ 95-1] = 0.0174219301594641737472;
    weight[ 96-1] = 0.0168019385741038652709;
    weight[ 97-1] = 0.0161732187295777199419;
    weight[ 98-1] = 0.0155367755558439824399;
    weight[ 99-1] = 0.0148936416648151820348;
    weight[100-1] = 0.0142448773729167743063;
    weight[101-1] = 0.0135915710097655467896;
    weight[102-1] = 0.0129348396636073734547;
    weight[103-1] = 0.0122758305600827700870;
    weight[104-1] = 0.0116157233199551347270;
    weight[105-1] = 0.0109557333878379016480;
    weight[106-1] = 0.0102971169579563555237;
    weight[107-1] = 0.00964117772970253669530;
    weight[108-1] = 0.00898927578406413572328;
    weight[109-1] = 0.00834283875396815770558;
    weight[110-1] = 0.00770337523327974184817;
    weight[111-1] = 0.00707248999543355546805;
    weight[112-1] = 0.00645190005017573692280;
    weight[113-1] = 0.00584344987583563950756;
    weight[114-1] = 0.00524912345480885912513;
    weight[115-1] = 0.00467105037211432174741;
    weight[116-1] = 0.00411150397865469304717;
    weight[117-1] = 0.00357289278351729964938;
    weight[118-1] = 0.00305775341017553113613;
    weight[119-1] = 0.00256876494379402037313;
    weight[120-1] = 0.00210881524572663287933;
    weight[121-1] = 0.00168114286542146990631;
    weight[122-1] = 0.00128952408261041739210;
    weight[123-1] = 0.000938369848542381500794;
    weight[124-1] = 0.000632607319362633544219;
    weight[125-1] = 0.000377746646326984660274;
    weight[126-1] = 0.000180739564445388357820;
    weight[127-1] = 0.0000505360952078625176247;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "PATTERSON_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal input value of ORDER.\n" );
    fprintf ( stderr, "  Order must be 1, 3, 7, 15, 31, 63, or 127.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

double r8_abs ( double x )

/******************************************************************************/
/*
  Purpose:
  
    R8_ABS returns the absolute value of an R8.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    14 November 2006
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, double X, the quantity whose absolute value is desired.
  
    Output, double R8_ABS, the absolute value of X.
*/
{
  double value;

  if ( 0.0 <= x )
  {
    value = x;
  }
  else
  {
    value = -x;
  }
  return value;
}
/******************************************************************************/

double r8_epsilon ( )

/******************************************************************************/
/*
  Purpose:
  
    R8_EPSILON returns the R8 roundoff unit.
  
  Discussion:
  
    The roundoff unit is a number R which is a power of 2 with the property that,
    to the precision of the computer's arithmetic,
      1 < 1 + R
    but
      1 = ( 1 + R / 2 )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    01 July 2004
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Output, double R8_EPSILON, the double precision round-off unit.
*/
{
  double r;

  r = 1.0;

  while ( 1.0 < ( double ) ( 1.0 + r )  )
  {
    r = r / 2.0;
  }

  return ( 2.0 * r );
}
/******************************************************************************/

double r8_factorial ( int n )

/******************************************************************************/
/*
  Purpose:
  
    R8_FACTORIAL computes the factorial of N, also denoted "N!".
  
  Formula:
  
    factorial ( N ) = N! = product ( 1 <= I <= N ) I
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    16 January 1999
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the argument of the factorial function.
    If N is less than 1, the function value is returned as 1.
  
    Output, double R8_FACTORIAL, the factorial of N.
*/
{
  int i;
  double value;

  value = 1.0;

  for ( i = 1; i <= n; i++ )
  {
    value = value * ( double ) ( i );
  }

  return value;
}
/******************************************************************************/

double r8_factorial2 ( int n )

/******************************************************************************/
/*
  Purpose:
  
    R8_FACTORIAL2 computes the double factorial function N!!
  
  Formula:
  
    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
  
  Example:
  
     N    N!!
  
     0     1
     1     1
     2     2
     3     3
     4     8
     5    15
     6    48
     7   105
     8   384
     9   945
    10  3840
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    22 January 2008
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the argument of the double factorial
    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
  
    Output, double R8_FACTORIAL2, the value of N!!.
*/
{
  int n_copy;
  double value;

  value = 1.0;

  if ( n < 1 )
  {
    return value;
  }

  n_copy = n;

  while ( 1 < n_copy )
  {
    value = value * ( double ) n_copy;
    n_copy = n_copy - 2;
  }

  return value;
}
/******************************************************************************/

double r8_gamma ( double x )

/******************************************************************************/
/*
  Purpose:
  
    R8_GAMMA evaluates Gamma(X) for a real argument.
  
  Discussion:
  
    This routine calculates the gamma function for a real argument X.
  
    Computation is based on an algorithm outlined in reference 1.
    The program uses rational functions that approximate the gamma
    function to at least 20 significant decimal digits.  Coefficients
    for the approximation over the interval (1,2) are unpublished.
    Those for the approximation for 12 <= X are from reference 2.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    18 January 2008
  
  Author:
  
    Original FORTRAN77 version by William Cody, Laura Stoltz.
    C version by John Burkardt.
  
  Reference:
  
    William Cody,
    An Overview of Software Development for Special Functions,
    in Numerical Analysis Dundee, 1975,
    edited by GA Watson,
    Lecture Notes in Mathematics 506,
    Springer, 1976.
  
    John Hart, Ward Cheney, Charles Lawson, Hans Maehly,
    Charles Mesztenyi, John Rice, Henry Thatcher,
    Christoph Witzgall,
    Computer Approximations,
    Wiley, 1968,
    LC: QA297.C64.
  
  Parameters:
  
    Input, double X, the argument of the function.
  
    Output, double R8_GAMMA, the value of the function.
*/
{
/*
  Coefficients for minimax approximation over (12, INF).
*/
  double c[7] = {
   -1.910444077728E-03,
    8.4171387781295E-04,
   -5.952379913043012E-04,
    7.93650793500350248E-04,
   -2.777777777777681622553E-03,
    8.333333333333333331554247E-02,
    5.7083835261E-03 };
  double eps = 2.22E-16;
  double fact;
  double half = 0.5;
  int i;
  int n;
  double one = 1.0;
  double p[8] = {
  -1.71618513886549492533811E+00,
   2.47656508055759199108314E+01,
  -3.79804256470945635097577E+02,
   6.29331155312818442661052E+02,
   8.66966202790413211295064E+02,
  -3.14512729688483675254357E+04,
  -3.61444134186911729807069E+04,
   6.64561438202405440627855E+04 };
  int parity;
  double pi = 3.1415926535897932384626434;
  double q[8] = {
  -3.08402300119738975254353E+01,
   3.15350626979604161529144E+02,
  -1.01515636749021914166146E+03,
  -3.10777167157231109440444E+03,
   2.25381184209801510330112E+04,
   4.75584627752788110767815E+03,
  -1.34659959864969306392456E+05,
  -1.15132259675553483497211E+05 };
  double res;
  double sqrtpi = 0.9189385332046727417803297;
  double sum;
  double twelve = 12.0;
  double two = 2.0;
  double value;
  double xbig = 171.624;
  double xden;
  double xinf = 1.79E+308;
  double xminin = 2.23E-308;
  double xnum;
  double y;
  double y1;
  double ysq;
  double z;
  double zero = 0.0;;

  parity = 0;
  fact = one;
  n = 0;
  y = x;
/*
  Argument is negative.
*/
  if ( y <= zero )
  {
    y = - x;
    y1 = ( double ) ( int ) ( y );
    res = y - y1;

    if ( res != zero )
    {
      if ( y1 != ( double ) ( int ) ( y1 * half ) * two )
      {
        parity = 1;
      }

      fact = - pi / sin ( pi * res );
      y = y + one;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Argument is positive.
*/
  if ( y < eps )
  {
/*
  Argument < EPS.
*/
    if ( xminin <= y )
    {
      res = one / y;
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
  else if ( y < twelve )
  {
    y1 = y;
/*
  0.0 < argument < 1.0.
*/
    if ( y < one )
    {
      z = y;
      y = y + one;
    }
/*
  1.0 < argument < 12.0.
  Reduce argument if necessary.
*/
    else
    {
      n = ( int ) ( y ) - 1;
      y = y - ( double ) ( n );
      z = y - one;
    }
/*
  Evaluate approximation for 1.0 < argument < 2.0.
*/
    xnum = zero;
    xden = one;
    for ( i = 0; i < 8; i++ )
    {
      xnum = ( xnum + p[i] ) * z;
      xden = xden * z + q[i];
    }
    res = xnum / xden + one;
/*
  Adjust result for case  0.0 < argument < 1.0.
*/
    if ( y1 < y )
    {
      res = res / y1;
    }
/*
  Adjust result for case 2.0 < argument < 12.0.
*/
    else if ( y < y1 )
    {
      for ( i = 1; i <= n; i++ )
      {
        res = res * y;
        y = y + one;
      }
    }
  }
  else
  {
/*
  Evaluate for 12.0 <= argument.
*/
    if ( y <= xbig )
    {
      ysq = y * y;
      sum = c[6];
      for ( i = 0; i < 6; i++ )
      {
        sum = sum / ysq + c[i];
      }
      sum = sum / y - y + sqrtpi;
      sum = sum + ( y - half ) * log ( y );
      res = exp ( sum );
    }
    else
    {
      res = xinf;
      value = res;
      return value;
    }
  }
/*
  Final adjustments and return.
*/
  if ( parity )
  {
    res = - res;
  }

  if ( fact != one )
  {
    res = fact / res;
  }

  value = res;

  return value;
}
/******************************************************************************/

double r8_psi ( double xx )

/******************************************************************************/
/*
  Purpose:
  
    R8_PSI evaluates the function Psi(X).
  
  Discussion:
  
    This routine evaluates the logarithmic derivative of the
    Gamma function,
  
      PSI(X) = d/dX ( GAMMA(X) ) / GAMMA(X)
             = d/dX LN ( GAMMA(X) )
  
    for real X, where either
  
      - XMAX1 < X < - XMIN, and X is not a negative integer,
  
    or
  
      XMIN < X.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    09 February 2008
  
  Author:
  
    Original FORTRAN77 version by William Cody.
    C version by John Burkardt.
  
  Reference:
  
    William Cody, Anthony Strecok, Henry Thacher,
    Chebyshev Approximations for the Psi Function,
    Mathematics of Computation,
    Volume 27, Number 121, January 1973, pages 123-127.
  
  Parameters:
  
    Input, double XX, the argument of the function.
  
    Output, double R8_PSI, the value of the function.
*/
{
  double aug;
  double den;
  double four = 4.0;
  double fourth = 0.25;
  double half = 0.5;
  int i;
  int n;
  int nq;
  double one = 1.0;
  double p1[9] = {
   4.5104681245762934160E-03,
   5.4932855833000385356,
   3.7646693175929276856E+02,
   7.9525490849151998065E+03,
   7.1451595818951933210E+04,
   3.0655976301987365674E+05,
   6.3606997788964458797E+05,
   5.8041312783537569993E+05,
   1.6585695029761022321E+05 };
  double p2[7] = {
  -2.7103228277757834192,
  -1.5166271776896121383E+01,
  -1.9784554148719218667E+01,
  -8.8100958828312219821,
  -1.4479614616899842986,
  -7.3689600332394549911E-02,
  -6.5135387732718171306E-21 };
  double piov4 = 0.78539816339744830962;
  double q1[8] = {
   9.6141654774222358525E+01,
   2.6287715790581193330E+03,
   2.9862497022250277920E+04,
   1.6206566091533671639E+05,
   4.3487880712768329037E+05,
   5.4256384537269993733E+05,
   2.4242185002017985252E+05,
   6.4155223783576225996E-08 };
  double q2[6] = {
   4.4992760373789365846E+01,
   2.0240955312679931159E+02,
   2.4736979003315290057E+02,
   1.0742543875702278326E+02,
   1.7463965060678569906E+01,
   8.8427520398873480342E-01 };
  double sgn;
  double three = 3.0;
  double upper;
  double value;
  double w;
  double x;
  double x01 = 187.0;
  double x01d = 128.0;
  double x02 = 6.9464496836234126266E-04;
  double xinf = 1.70E+38;
  double xlarge = 2.04E+15;
  double xmax1 = 3.60E+16;
  double xmin1 = 5.89E-39;
  double xsmall = 2.05E-09;
  double z;
  double zero = 0.0;

  x = xx;
  w = r8_abs ( x );
  aug = zero;
/*
  Check for valid arguments, then branch to appropriate algorithm.
*/
  if ( xmax1 <= - x || w < xmin1 )
  {
    if ( zero < x )
    {
      value = - xinf;
    }
    else
    {
      value = xinf;
    }
    return value;
  }

  if ( x < half )
  {
/*
  X < 0.5, use reflection formula: psi(1-x) = psi(x) + pi * cot(pi*x)
  Use 1/X for PI*COTAN(PI*X)  when  XMIN1 < |X| <= XSMALL.
*/
    if ( w <= xsmall )
    {
      aug = - one / x;
    }
/*
  Argument reduction for cotangent.
*/
    else
    {
      if ( x < zero )
      {
        sgn = piov4;
      }
      else
      {
        sgn = - piov4;
      }

      w = w - ( double ) ( ( int ) ( w ) );
      nq = (int) ( w * four );
      w = four * ( w - ( double ) ( nq ) * fourth );
/*
  W is now related to the fractional part of 4.0 * X.
  Adjust argument to correspond to values in the first
  quadrant and determine the sign.
*/
      n = nq / 2;

      if ( n + n != nq )
      {
        w = one - w;
      }

      z = piov4 * w;

      if ( ( n % 2 ) != 0 )
      {
        sgn = - sgn;
      }
/*
  Determine the final value for  -pi * cotan(pi*x).
*/
      n = ( nq + 1 ) / 2;
      if ( ( n % 2 ) == 0 )
      {
/*
  Check for singularity.
*/
        if ( z == zero )
        {
          if ( zero < x )
          {
            value = -xinf;
          }
          else
          {
            value = xinf;
          }
          return value;
        }
        aug = sgn * ( four / tan ( z ) );
      }
      else
      {
        aug = sgn * ( four * tan ( z ) );
      }
    }
    x = one - x;
  }
/*
  0.5 <= X <= 3.0.
*/
  if ( x <= three )
  {
    den = x;
    upper = p1[0] * x;
    for ( i = 1; i <= 7; i++ )
    {
      den = ( den + q1[i-1] ) * x;
      upper = ( upper + p1[i]) * x;
    }
    den = ( upper + p1[8] ) / ( den + q1[7] );
    x = ( x - x01 / x01d ) - x02;
    value = den * x + aug;
    return value;
  }
/*
  3.0 < X.
*/
  if ( x < xlarge )
  {
    w = one / ( x * x );
    den = w;
    upper = p2[0] * w;
    for ( i = 1; i <= 5; i++ )
    {
      den = ( den + q2[i-1] ) * w;
      upper = ( upper + p2[i] ) * w;
    }
    aug = ( upper + p2[6] ) / ( den + q2[5] ) - half / x + aug;
  }

  value = aug + log ( x );

  return value;
}
/******************************************************************************/

double r8_huge ( )

/******************************************************************************/
/*
  Purpose:
  
    R8_HUGE returns a "huge" R8.
  
  Discussion:
  
    The value returned by this function is NOT required to be the
    maximum representable R8.  This value varies from machine to machine,
    from compiler to compiler, and may cause problems when being printed.
    We simply want a "very large" but non-infinite number.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    06 October 2007
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Output, double R8_HUGE, a "huge" R8 value.
*/
{
  double value;

  value = 1.0E+30;

  return value;
}
/******************************************************************************/

double r8_hyper_2f1 ( double a, double b, double c, double x )

/******************************************************************************/
/*
  Purpose:
  
    R8_HYPER_2F1 evaluates the hypergeometric function 2F1(A,B,C,X).
  
  Discussion:
  
    A minor bug was corrected.  The HW variable, used in several places as
    the "old" value of a quantity being iteratively improved, was not
    being initialized.  JVB, 11 February 2008.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    05 July 2009
  
  Author:
  
    Original FORTRAN77 version by Shanjie Zhang, Jianming Jin.
    C version by John Burkardt.
  
    The original FORTRAN77 version of this routine is copyrighted by
    Shanjie Zhang and Jianming Jin.  However, they give permission to
    incorporate this routine into a user program provided that the copyright
    is acknowledged.
  
  Reference:
  
    Shanjie Zhang, Jianming Jin,
    Computation of Special Functions,
    Wiley, 1996,
    ISBN: 0-471-11963-6,
    LC: QA351.C45
  
  Parameters:
  
    Input, double A, B, C, X, the arguments of the function.
    C must not be equal to a nonpositive integer.
    X < 1.
  
    Output, double R8_HYPER_2F1, the value of the function.
*/
{
  double a0;
  double aa;
  double bb;
  double c0;
  double c1;
  double el = 0.5772156649015329;
  double eps;
  double f0;
  double f1;
  double g0;
  double g1;
  double g2;
  double g3;
  double ga;
  double gabc;
  double gam;
  double gb;
  double gbm;
  double gc;
  double gca;
  double gcab;
  double gcb;
  double gm;
  double hf;
  double hw;
  int j;
  int k;
  int l0;
  int l1;
  int l2;
  int l3;
  int l4;
  int l5;
  int m;
  int nm;
  double pa;
  double pb;
  double pi = 3.141592653589793;
  double r;
  double r0;
  double r1;
  double rm;
  double rp;
  double sm;
  double sp;
  double sp0;
  double x1;

  l0 = ( c == ( int ) ( c ) ) && ( c < 0.0 );
  l1 = ( 1.0 - x < 1.0E-15 ) && ( c - a - b <= 0.0 );
  l2 = ( a == ( int ) ( a ) ) && ( a < 0.0 );
  l3 = ( b == ( int ) ( b ) ) && ( b < 0.0 );
  l4 = ( c - a == ( int ) ( c - a ) ) && ( c - a <= 0.0 );
  l5 = ( c - b == ( int ) ( c - b ) ) && ( c - b <= 0.0 );

  if ( l0 || l1 )
  {
    hf = 0.0;
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "R8_HYPER_2F1 - Fatal error!\n" );
    fprintf ( stderr, "  The hypergeometric series is divergent.\n" );
    exit ( 1 );
  }

  if ( 0.95 < x )
  {
    eps = 1.0E-08;
  }
  else
  {
    eps = 1.0E-15;
  }

  if ( x == 0.0 || a == 0.0 || b == 0.0 )
  {
    hf = 1.0;
    return hf;
  }
  else if ( 1.0 - x == eps && 0.0 < c - a - b )
  {
    gc = r8_gamma ( c );
    gcab = r8_gamma ( c - a - b );
    gca = r8_gamma ( c - a );
    gcb = r8_gamma ( c - b );
    hf = gc * gcab / ( gca * gcb );
    return hf;
  }
  else if ( 1.0 + x <= eps && r8_abs ( c - a + b - 1.0 ) <= eps )
  {
    g0 = sqrt ( pi ) * pow ( 2.0, - a );
    g1 = r8_gamma ( c );
    g2 = r8_gamma ( 1.0 + a / 2.0 - b );
    g3 = r8_gamma ( 0.5 + 0.5 * a );
    hf = g0 * g1 / ( g2 * g3 );
    return hf;
  }
  else if ( l2 || l3 )
  {
    if ( l2 )
    {
      nm = ( int ) ( r8_abs ( a ) );
    }

    if ( l3 )
    {
      nm = ( int ) ( r8_abs ( b ) );
    }

    hf = 1.0;
    r = 1.0;

    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 )
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }

    return hf;
  }
  else if ( l4 || l5 )
  {
    if ( l4 )
    {
      nm = ( int ) ( r8_abs ( c - a ) );
    }

    if ( l5 )
    {
      nm = ( int ) ( r8_abs ( c - b ) );
    }

    hf = 1.0;
    r  = 1.0;
    for ( k = 1; k <= nm; k++ )
    {
      r = r * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 )
        / ( k * ( c + k - 1.0 ) ) * x;
      hf = hf + r;
    }
    hf = pow ( 1.0 - x, c - a - b ) * hf;
    return hf;
  }

  aa = a;
  bb = b;
  x1 = x;

  if ( x < 0.0 )
  {
    x = x / ( x - 1.0 );
    if ( a < c && b < a && 0.0 < b )
    {
      a = bb;
      b = aa;
    }
    b = c - b;
  }

  if ( 0.75 <= x )
  {
    gm = 0.0;

    if ( r8_abs ( c - a - b - ( int ) ( c - a - b ) ) < 1.0E-15 )
    {
      m = (int) ( c - a - b );
      ga = r8_gamma ( a );
      gb = r8_gamma ( b );
      gc = r8_gamma ( c );
      gam = r8_gamma ( a + m );
      gbm = r8_gamma ( b + m );

      pa = r8_psi ( a );
      pb = r8_psi ( b );

      if ( m != 0 )
      {
        gm = 1.0;
      }

      for ( j = 1; j <= abs ( m ) - 1; j++ )
      {
        gm = gm * j;
      }

      rm = 1.0;
      for ( j = 1; j <= abs ( m ); j++ )
      {
        rm = rm * j;
      }

      f0 = 1.0;
      r0 = 1.0;;
      r1 = 1.0;
      sp0 = 0.0;;
      sp = 0.0;

      if ( 0 <= m )
      {
        c0 = gm * gc / ( gam * gbm );
        c1 = - gc * pow ( x - 1.0, m ) / ( ga * gb * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 )
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( a + k - 1.0 ) + 1.0 / ( b + k - 1.0 )
          - 1.0 / ( double ) ( k );
        }

        f1 = pa + pb + sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a ) / ( k * ( a + k - 1.0 ) )
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + ( 1.0 - a )
              / ( ( j + k ) * ( a + j + k - 1.0 ) )
              + 1.0 / ( b + j + k - 1.0 );
          }

          rp = pa + pb + 2.0 * el + sp + sm + log ( 1.0 - x );

          r1 = r1 * ( a + m + k - 1.0 ) * ( b + m + k - 1.0 )
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }
          hw = f1;
        }
        hf = f0 * c0 + f1 * c1;
      }
      else if ( m < 0 )
      {
        m = - m;
        c0 = gm * gc / ( ga * gb * pow ( 1.0 - x, m ) );
        c1 = - pow ( - 1.0, m ) * gc / ( gam * gbm * rm );

        for ( k = 1; k <= m - 1; k++ )
        {
          r0 = r0 * ( a - m + k - 1.0 ) * ( b - m + k - 1.0 )
            / ( k * ( k - m ) ) * ( 1.0 - x );
          f0 = f0 + r0;
        }

        for ( k = 1; k <= m; k++ )
        {
          sp0 = sp0 + 1.0 / ( double ) ( k );
        }

        f1 = pa + pb - sp0 + 2.0 * el + log ( 1.0 - x );
        hw = f1;

        for ( k = 1; k <= 250; k++ )
        {
          sp = sp + ( 1.0 - a )
            / ( k * ( a + k - 1.0 ) )
            + ( 1.0 - b ) / ( k * ( b + k - 1.0 ) );

          sm = 0.0;
          for ( j = 1; j <= m; j++ )
          {
            sm = sm + 1.0 / ( double ) ( j + k );
          }

          rp = pa + pb + 2.0 * el + sp - sm + log ( 1.0 - x );

          r1 = r1 * ( a + k - 1.0 ) * ( b + k - 1.0 )
            / ( k * ( m + k ) ) * ( 1.0 - x );

          f1 = f1 + r1 * rp;

          if ( r8_abs ( f1 - hw ) < r8_abs ( f1 ) * eps )
          {
            break;
          }

          hw = f1;
        }

        hf = f0 * c0 + f1 * c1;
      }
    }
    else
    {
      ga = r8_gamma ( a );
      gb = r8_gamma ( b );
      gc = r8_gamma ( c );
      gca = r8_gamma ( c - a );
      gcb = r8_gamma ( c - b );
      gcab = r8_gamma ( c - a - b );
      gabc = r8_gamma ( a + b - c );
      c0 = gc * gcab / ( gca * gcb );
      c1 = gc * gabc / ( ga * gb ) * pow ( 1.0 - x, c - a - b );
      hf = 0.0;
      hw = hf;
      r0 = c0;
      r1 = c1;

      for ( k = 1; k <= 250; k++ )
      {
        r0 = r0 * ( a + k - 1.0 ) * ( b + k - 1.0 )
          / ( k * ( a + b - c + k ) ) * ( 1.0 - x );

        r1 = r1 * ( c - a + k - 1.0 ) * ( c - b + k - 1.0 )
          / ( k * ( c - a - b + k ) ) * ( 1.0 - x );

        hf = hf + r0 + r1;

        if ( r8_abs ( hf - hw ) < r8_abs ( hf ) * eps )
        {
          break;
        }
        hw = hf;
      }
      hf = hf + c0 + c1;
    }
  }
  else
  {
    a0 = 1.0;

    if ( a < c && c < 2.0 * a && b < c && c < 2.0 * b )
    {
      a0 = pow ( 1.0 - x, c - a - b );
      a = c - a;
      b = c - b;
    }

    hf = 1.0;
    hw = hf;
    r = 1.0;

    for ( k = 1; k <= 250; k++ )
    {
      r = r * ( a + k - 1.0 ) * ( b + k - 1.0 )
        / ( k * ( c + k - 1.0 ) ) * x;

      hf = hf + r;

      if ( r8_abs ( hf - hw ) <= r8_abs ( hf ) * eps )
      {
        break;
      }

      hw = hf;
    }
    hf = a0 * hf;
  }

  if ( x1 < 0.0 )
  {
    x = x1;
    c0 = 1.0 / pow ( 1.0 - x, aa );
    hf = c0 * hf;
  }

  a = aa;
  b = bb;

  if ( 120 < k )
  {
    fprintf (stderr, "\nR8_HYPER_2F1 - Warning!\n  A large number of iterations were needed.\n  The accuracy of the results should be checked.\n");
  }

  return hf;
}
/******************************************************************************/

double r8_max ( double x, double y )

/******************************************************************************/
/*
  Purpose:
  
    R8_MAX returns the maximum of two R8's.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    18 August 2004
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, double X, Y, the quantities to compare.
  
    Output, double R8_MAX, the maximum of X and Y.
*/
{
  double value;

  if ( y < x )
  {
    value = x;
  }
  else
  {
    value = y;
  }
  return value;
}
/******************************************************************************/

void r8vec_copy ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:
  
    R8VEC_COPY copies an R8VEC.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    03 July 2005
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the number of entries in the vectors.
  
    Input, double A1[N], the vector to be copied.
  
    Input, double A2[N], the copy of A1.
*/
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    a2[i] = a1[i];
  }
  return;
}
/******************************************************************************/

double r8vec_dot ( int n, double a1[], double a2[] )

/******************************************************************************/
/*
  Purpose:
  
    R8VEC_DOT computes the dot product of a pair of R8VEC's.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    03 July 2005
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the number of entries in the vectors.
  
    Input, double A1[N], A2[N], the two vectors to be considered.
  
    Output, double R8VEC_DOT, the dot product of the vectors.
*/
{
  int i;
  double value;

  value = 0.0;
  for ( i = 0; i < n; i++ )
  {
    value = value + a1[i] * a2[i];
  }
  return value;
}
/******************************************************************************/

void r8vec_reverse ( int n, double a[] )

/******************************************************************************/
/*
  Purpose:
  
    R8VEC_REVERSE reverses the elements of an R8VEC.
  
  Discussion:
  
    An R8VEC is a vector of double precision values.
  
    Input:
  
      N = 5,
      A = ( 11.0, 12.0, 13.0, 14.0, 15.0 )
  
    Output:
  
      A = ( 15.0, 14.0, 13.0, 12.0, 11.0 )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    30 April 2006
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, int N, the number of entries in the array.
  
    Input/output, double A[N], the array to be reversed.
*/
{
  int i;
  double temp;

  for ( i = 0; i < n/2; i++ )
  {
    temp     = a[i];
    a[i]     = a[n-1-i];
    a[n-1-i] = temp;
  }
  return;
}
/******************************************************************************/

void radau_compute ( int n, double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    RADAU_COMPUTE computes a Radau quadrature rule.
  
  Discussion:
  
    The Radau rule is distinguished by the fact that the left endpoint
    (-1) is always an abscissa.
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x) = 1.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= NORDER ) WEIGHT(I) * F ( XTAB(I) )
  
    The quadrature rule will integrate exactly all polynomials up to
    X**(2*NORDER-2).
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 August 2007
  
  Author:
  
    Original MATLAB version by Greg von Winckel.
    C version by John Burkardt.
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Claudio Canuto, Yousuff Hussaini, Alfio Quarteroni, Thomas Zang,
    Spectral Methods in Fluid Dynamics,
    Springer, 1993,
    ISNB13: 978-3540522058,
    LC: QA377.S676.
  
    Francis Hildebrand,
    Section 8.11,
    Introduction to Numerical Analysis,
    Dover, 1987,
    ISBN13: 978-0486653631,
    LC: QA300.H5.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3,
    LC: QA47.M315.
  
  Parameters:
  
    Input, int N, the order of the rule.  N must be at least 1.
  
    Output, double X[N], W[N], the abscissas and weights
    of the rule.
*/
{
  int i;
  int iterate;
  int iterate_max = 25;
  int j;
  double pi = 3.141592653589793;
  double temp;
  double test;
  double tolerance;

  if ( n < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "RADAU_COMPUTE - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", n );
    fprintf ( stderr, "  ORDER must be at least 1.\n" );
    exit ( 1 );
  }

  tolerance = 100.0 * r8_epsilon ( );
/*
  Initial estimate for the abscissas is the Chebyshev-Gauss-Radau nodes.
*/
  for ( i = 0; i < n; i++ )
  {
    x[i] = - cos ( 2.0 * pi * ( double ) (         i )
                            / ( double ) ( 2 * n - 1 ) );
  }
  double xold[n];
  double p[n*(n+1)];
  iterate = 0;

  do
  {
    for ( i = 0; i < n; i++ )
    {
      xold[i] = x[i];
    }

    temp = 1.0;
    for ( j = 0; j < n + 1; j++ )
    {
      p[0+j*n] = temp;
      temp = -temp;
    }

    for ( i = 1; i < n; i++ )
    {
      p[i+0*n] = 1.0;
    }
    for ( i = 1; i < n; i++ )
    {
      p[i+1*n] = x[i];
    }

    for ( j = 2; j <= n; j++ )
    {
      for ( i = 1; i < n; i++ )
      {
        p[i+j*n] = ( ( double ) ( 2 * j - 1 ) * x[i] * p[i+(j-1)*n]
                   + ( double ) (   - j + 1 ) *        p[i+(j-2)*n] )
                   / ( double ) (     j     );
      }
    }
    for ( i = 1; i < n; i++ )
    {
      x[i] = xold[i] - ( ( 1.0 - xold[i] ) / ( double ) ( n ) )
        * ( p[i+(n-1)*n] + p[i+n*n] ) / ( p[i+(n-1)*n] - p[i+n*n] );
    }
    test = 0.0;
    for ( i = 0; i < n; i++ )
    {
      test = r8_max ( test, r8_abs ( x[i] - xold[i] ) );
    }
    iterate = iterate + 1;
  } while ( tolerance < test && iterate < iterate_max );

  w[0] = 2.0 / ( double ) ( n * n );
  for ( i = 1; i < n; i++ )
  {
    w[i] = ( 1.0 - x[i] ) / pow ( ( double ) ( n ) * p[i+(n-1)*n], 2 );
  }

  return;
}
/******************************************************************************/

void radau_set ( int order, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    RADAU_SET sets abscissas and weights for Radau quadrature.
  
  Discussion:
  
    The Radau rule is distinguished by the fact that the left endpoint
    (-1) is always an abscissa.
  
    The integration interval is [ -1, 1 ].
  
    The weight function is w(x-1] = 1.0.
  
    The integral to approximate:
  
      Integral ( -1 <= X <= 1 ) F(X) dX
  
    The quadrature rule:
  
      Sum ( 1 <= I <= ORDER ) WEIGHt[I) * F ( XTAb[I) )
  
    The quadrature rule will integrate exactly all polynomials up to
    X**(2*ORDER-2).
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    01 May 2006
  
  Author:
  
    John Burkardt
  
  Reference:
  
    Milton Abramowitz, Irene Stegun,
    Handbook of Mathematical Functions,
    National Bureau of Standards, 1964,
    ISBN: 0-486-61272-4,
    LC: QA47.A34.
  
    Arthur Stroud, Don Secrest,
    Gaussian Quadrature Formulas,
    Prentice Hall, 1966,
    LC: QA299.4G3S7.
  
    Daniel Zwillinger, editor,
    CRC Standard Mathematical Tables and Formulae,
    30th Edition,
    CRC Press, 1996,
    ISBN: 0-8493-2479-3.
  
  Parameters:
  
    Input, int ORDER, the order of the rule.
    ORDER must be between 1 and 15.
  
    Output, double XTAB[ORDER], the abscissas of the rule.
  
    Output, double WEIGHT[ORDER], the weights of the rule.
    The weights are positive.  The weights are not symmetric.
    The weights should sum to 2.  WEIGHt[1) should equal 2 / ORDER**2.
*/
{
  if ( order == 1 )
  {
    xtab[1-1] =   - 1.0E+00;
    weight[1-1] =   2.0E+00;
  }
  else if ( order == 2 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =    1.0E+00 / 3.0E+00;

    weight[1-1] =  0.5E+00;
    weight[2-1] =  1.5E+00;
  }
  else if ( order == 3 )
  {
    xtab[1-1] =   - 1.0E+00;
    xtab[2-1] =   - 0.289897948556635619639456814941E+00;
    xtab[3-1] =     0.689897948556635619639456814941E+00;

    weight[1-1] =  0.222222222222222222222222222222E+00;
    weight[2-1] =  0.102497165237684322767762689304E+01;
    weight[3-1] =  0.752806125400934550100150884739E+00;
  }
  else if ( order == 4 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.575318923521694112050483779752E+00;
    xtab[3-1] =    0.181066271118530578270147495862E+00;
    xtab[4-1] =    0.822824080974592105208907712461E+00;

    weight[1-1] =  0.125E+00;
    weight[2-1] =  0.657688639960119487888578442146E+00;
    weight[3-1] =  0.776386937686343761560464613780E+00;
    weight[4-1] =  0.440924422353536750550956944074E+00;
  }
  else if ( order == 5 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.720480271312438895695825837750E+00;
    xtab[3-1] =  - 0.167180864737833640113395337326E+00;
    xtab[4-1] =    0.446313972723752344639908004629E+00;
    xtab[5-1] =    0.885791607770964635613757614892E+00;

    weight[1-1] =  0.08E+00;
    weight[2-1] =  0.446207802167141488805120436457E+00;
    weight[3-1] =  0.623653045951482508163709823153E+00;
    weight[4-1] =  0.562712030298924120384345300681E+00;
    weight[5-1] =  0.287427121582451882646824439708E+00;
  }
  else if ( order == 6 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.802929828402347147753002204224E+00;
    xtab[3-1] =  - 0.390928546707272189029229647442E+00;
    xtab[4-1] =    0.124050379505227711989974959990E+00;
    xtab[5-1] =    0.603973164252783654928415726409E+00;
    xtab[6-1] =    0.920380285897062515318386619813E+00;

    weight[1-1] =  0.555555555555555555555555555556E-01;
    weight[2-1] =  0.319640753220510966545779983796E+00;
    weight[3-1] =  0.485387188468969916159827915587E+00;
    weight[4-1] =  0.520926783189574982570229406570E+00;
    weight[5-1] =  0.416901334311907738959406382743E+00;
    weight[6-1] =  0.201588385253480840209200755749E+00;
  }
  else if ( order == 7 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.853891342639482229703747931639E+00;
    xtab[3-1] =  - 0.538467724060109001833766720231E+00;
    xtab[4-1] =  - 0.117343037543100264162786683611E+00;
    xtab[5-1] =    0.326030619437691401805894055838E+00;
    xtab[6-1] =    0.703842800663031416300046295008E+00;
    xtab[7-1] =    0.941367145680430216055899446174E+00;

    weight[1-1] =  0.408163265306122448979591836735E-01;
    weight[2-1] =  0.239227489225312405787077480770E+00;
    weight[3-1] =  0.380949873644231153805938347876E+00;
    weight[4-1] =  0.447109829014566469499348953642E+00;
    weight[5-1] =  0.424703779005955608398308039150E+00;
    weight[6-1] =  0.318204231467301481744870434470E+00;
    weight[7-1] =  0.148988471112020635866497560418E+00;
  }
  else if ( order == 8 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.887474878926155707068695617935E+00;
    xtab[3-1] =  - 0.639518616526215270024840114382E+00;
    xtab[4-1] =  - 0.294750565773660725252184459658E+00;
    xtab[5-1] =    0.943072526611107660028971153047E-01;
    xtab[6-1] =    0.468420354430821063046421216613E+00;
    xtab[7-1] =    0.770641893678191536180719525865E+00;
    xtab[8-1] =    0.955041227122575003782349000858E+00;

    weight[1-1] =  0.03125E+00;
    weight[2-1] =  0.185358154802979278540728972699E+00;
    weight[3-1] =  0.304130620646785128975743291400E+00;
    weight[4-1] =  0.376517545389118556572129261442E+00;
    weight[5-1] =  0.391572167452493593082499534004E+00;
    weight[6-1] =  0.347014795634501280228675918422E+00;
    weight[7-1] =  0.249647901329864963257869293513E+00;
    weight[8-1] =  0.114508814744257199342353728520E+00;
  }
  else if ( order == 9 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.910732089420060298533757956283E+00;
    xtab[3-1] =  - 0.711267485915708857029562959544E+00;
    xtab[4-1] =  - 0.426350485711138962102627520502E+00;
    xtab[5-1] =  - 0.903733696068532980645444599064E-01;
    xtab[6-1] =    0.256135670833455395138292079035E+00;
    xtab[7-1] =    0.571383041208738483284917464837E+00;
    xtab[8-1] =    0.817352784200412087992517083851E+00;
    xtab[9-1] =    0.964440169705273096373589797925E+00;

    weight[1-1] =  0.246913580246913580246913580247E-01;
    weight[2-1] =  0.147654019046315385819588499802E+00;
    weight[3-1] =  0.247189378204593052361239794969E+00;
    weight[4-1] =  0.316843775670437978338000849642E+00;
    weight[5-1] =  0.348273002772966594071991031186E+00;
    weight[6-1] =  0.337693966975929585803724239792E+00;
    weight[7-1] =  0.286386696357231171146705637752E+00;
    weight[8-1] =  0.200553298024551957421165090417E+00;
    weight[9-1] =  0.907145049232829170128934984159E-01;
  }
  else if ( order == 10 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.927484374233581078117671398464E+00;
    xtab[3-1] =  - 0.763842042420002599615429776011E+00;
    xtab[4-1] =  - 0.525646030370079229365386614293E+00;
    xtab[5-1] =  - 0.236234469390588049278459503207E+00;
    xtab[6-1] =    0.760591978379781302337137826389E-01;
    xtab[7-1] =    0.380664840144724365880759065541E+00;
    xtab[8-1] =    0.647766687674009436273648507855E+00;
    xtab[9-1] =    0.851225220581607910728163628088E+00;
    xtab[10-1] =   0.971175180702246902734346518378E+00;

    weight[1-1] =  0.02E+00;
    weight[2-1] =  0.120296670557481631517310522702E+00;
    weight[3-1] =  0.204270131879000675555788672223E+00;
    weight[4-1] =  0.268194837841178696058554475262E+00;
    weight[5-1] =  0.305859287724422621016275475401E+00;
    weight[6-1] =  0.313582457226938376695902847302E+00;
    weight[7-1] =  0.290610164832918311146863077963E+00;
    weight[8-1] =  0.239193431714379713376571966160E+00;
    weight[9-1] =  0.164376012736921475701681668908E+00;
    weight[10-1] = 0.736170054867584989310512940790E-01;
  }
  else if ( order == 11 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.939941935677027005913871284731E+00;
    xtab[3-1] =  - 0.803421975580293540697597956820E+00;
    xtab[4-1] =  - 0.601957842073797690275892603234E+00;
    xtab[5-1] =  - 0.351888923353330214714301017870E+00;
    xtab[6-1] =  - 0.734775314313212657461903554238E-01;
    xtab[7-1] =    0.210720306228426314076095789845E+00;
    xtab[8-1] =    0.477680647983087519467896683890E+00;
    xtab[9-1] =    0.705777100713859519144801128840E+00;
    xtab[10-1] =   0.876535856245703748954741265611E+00;
    xtab[11-1] =   0.976164773135168806180508826082E+00;

    weight[1-1] =  0.165289256198347107438016528926E-01;
    weight[2-1] =  0.998460819079680638957534695802E-01;
    weight[3-1] =  0.171317619206659836486712649042E+00;
    weight[4-1] =  0.228866123848976624401683231126E+00;
    weight[5-1] =  0.267867086189684177806638163355E+00;
    weight[6-1] =  0.285165563941007337460004408915E+00;
    weight[7-1] =  0.279361333103383045188962195720E+00;
    weight[8-1] =  0.250925377697128394649140267633E+00;
    weight[9-1] =  0.202163108540024418349931754266E+00;
    weight[10-1] = 0.137033682133202256310153880580E+00;
    weight[11-1] = 0.609250978121311347072183268883E-01;
  }
  else if ( order == 12 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.949452759204959300493337627077E+00;
    xtab[3-1] =  - 0.833916773105189706586269254036E+00;
    xtab[4-1] =  - 0.661649799245637148061133087811E+00;
    xtab[5-1] =  - 0.444406569781935851126642615609E+00;
    xtab[6-1] =  - 0.196994559534278366455441427346E+00;
    xtab[7-1] =    0.637247738208319158337792384845E-01;
    xtab[8-1] =    0.319983684170669623532789532206E+00;
    xtab[9-1] =    0.554318785912324288984337093085E+00;
    xtab[10-1] =   0.750761549711113852529400825472E+00;
    xtab[11-1] =   0.895929097745638894832914608454E+00;
    xtab[12-1] =   0.979963439076639188313950540264E+00;

    weight[1-1] =  0.138888888888888888888888888888E-01;
    weight[2-1] =  0.841721349386809762415796536813E-01;
    weight[3-1] =  0.145563668853995128522547654706E+00;
    weight[4-1] =  0.196998534826089634656049637969E+00;
    weight[5-1] =  0.235003115144985839348633985940E+00;
    weight[6-1] =  0.256991338152707776127974253598E+00;
    weight[7-1] =  0.261465660552133103438074715743E+00;
    weight[8-1] =  0.248121560804009959403073107079E+00;
    weight[9-1] =  0.217868879026192438848747482023E+00;
    weight[10-1] = 0.172770639313308564306065766966E+00;
    weight[11-1] = 0.115907480291738392750341908272E+00;
    weight[12-1] = 0.512480992072692974680229451351E-01;
  }
  else if ( order == 13 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.956875873668299278183813833834E+00;
    xtab[3-1] =  - 0.857884202528822035697620310269E+00;
    xtab[4-1] =  - 0.709105087529871761580423832811E+00;
    xtab[5-1] =  - 0.519197779050454107485205148087E+00;
    xtab[6-1] =  - 0.299201300554509985532583446686E+00;
    xtab[7-1] =  - 0.619016986256353412578604857936E-01;
    xtab[8-1] =    0.178909837597084635021931298881E+00;
    xtab[9-1] =    0.409238231474839556754166331248E+00;
    xtab[10-1] =   0.615697890940291918017885487543E+00;
    xtab[11-1] =   0.786291018233046684731786459135E+00;
    xtab[12-1] =   0.911107073689184553949066402429E+00;
    xtab[13-1] =   0.982921890023145161262671078244E+00;

    weight[1-1] =  0.118343195266272189349112426036E-01;
    weight[2-1] =  0.719024162924955289397537405641E-01;
    weight[3-1] =  0.125103834331152358133769287976E+00;
    weight[4-1] =  0.171003460470616642463758674512E+00;
    weight[5-1] =  0.206960611455877074631132560829E+00;
    weight[6-1] =  0.230888862886995434012203758668E+00;
    weight[7-1] =  0.241398342287691148630866924129E+00;
    weight[8-1] =  0.237878547660712031342685189180E+00;
    weight[9-1] =  0.220534229288451464691077164199E+00;
    weight[10-1] = 0.190373715559631732254759820746E+00;
    weight[11-1] = 0.149150950090000205151491864242E+00;
    weight[12-1] = 0.992678068818470859847363877478E-01;
    weight[13-1] = 0.437029032679020748288533846051E-01;
  }
  else if ( order == 14 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.962779269978024297120561244319E+00;
    xtab[3-1] =  - 0.877048918201462024795266773531E+00;
    xtab[4-1] =  - 0.747389642613378838735429134263E+00;
    xtab[5-1] =  - 0.580314056546874971105726664999E+00;
    xtab[6-1] =  - 0.384202003439203313794083903375E+00;
    xtab[7-1] =  - 0.168887928042680911008441695622E+00;
    xtab[8-1] =    0.548312279917645496498107146428E-01;
    xtab[9-1] =    0.275737205435522399182637403545E+00;
    xtab[10-1] =   0.482752918588474966820418534355E+00;
    xtab[11-1] =   0.665497977216884537008955042481E+00;
    xtab[12-1] =   0.814809550601994729434217249123E+00;
    xtab[13-1] =   0.923203722520643299246334950272E+00;
    xtab[14-1] =   0.985270697947821356698617003172E+00;

    weight[1-1] =  0.102040816326530612244897959184E-01;
    weight[2-1] =  0.621220169077714601661329164668E-01;
    weight[3-1] =  0.108607722744362826826720935229E+00;
    weight[4-1] =  0.149620539353121355950520836946E+00;
    weight[5-1] =  0.183127002125729654123867302103E+00;
    weight[6-1] =  0.207449763335175672668082886489E+00;
    weight[7-1] =  0.221369811499570948931671683021E+00;
    weight[8-1] =  0.224189348002707794238414632220E+00;
    weight[9-1] =  0.215767100604618851381187446115E+00;
    weight[10-1] = 0.196525518452982430324613091930E+00;
    weight[11-1] = 0.167429727891086278990102277038E+00;
    weight[12-1] = 0.129939668737342347807425737146E+00;
    weight[13-1] = 0.859405354429804030893077310866E-01;
    weight[14-1] = 0.377071632698969142774627282919E-01;
  }
  else if ( order == 15 )
  {
    xtab[1-1] =  - 1.0E+00;
    xtab[2-1] =  - 0.967550468197200476562456018282E+00;
    xtab[3-1] =  - 0.892605400120550767066811886849E+00;
    xtab[4-1] =  - 0.778685617639031079381743321893E+00;
    xtab[5-1] =  - 0.630779478886949283946148437224E+00;
    xtab[6-1] =  - 0.455352905778529370872053455981E+00;
    xtab[7-1] =  - 0.260073376740807915768961188263E+00;
    xtab[8-1] =  - 0.534757226797460641074538896258E-01;
    xtab[9-1] =    0.155410685384859484319182024964E+00;
    xtab[10-1] =   0.357456512022127651195319205174E+00;
    xtab[11-1] =   0.543831458701484016930711802760E+00;
    xtab[12-1] =   0.706390264637572540152679669478E+00;
    xtab[13-1] =   0.838029000636089631215097384520E+00;
    xtab[14-1] =   0.932997190935973719928072142859E+00;
    xtab[15-1] =   0.987166478414363086378359071811E+00;

    weight[1-1] =  0.888888888888888888888888888889E-02;
    weight[2-1] =  0.542027800486444943382142368018E-01;
    weight[3-1] =  0.951295994604808992038477266346E-01;
    weight[4-1] =  0.131875462504951632186262157944E+00;
    weight[5-1] =  0.162854477303832629448732245828E+00;
    weight[6-1] =  0.186715145839450908083795103799E+00;
    weight[7-1] =  0.202415187030618429872703310435E+00;
    weight[8-1] =  0.209268608147694581430889790306E+00;
    weight[9-1] =  0.206975960249553755479027321787E+00;
    weight[10-1] = 0.195637503045116116473556617575E+00;
    weight[11-1] = 0.175748872642447685670310440476E+00;
    weight[12-1] = 0.148179527003467253924682058743E+00;
    weight[13-1] = 0.114135203489752753013075582569E+00;
    weight[14-1] = 0.751083927605064397329716653914E-01;
    weight[15-1] = 0.328643915845935322530428528231E-01;
  }
  else
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "RADAU_SET - Fatal error!\n" );
    fprintf ( stderr, "  Illegal value of ORDER = %d\n", order );
    fprintf ( stderr, "  Legal values are 1 to 15.\n" );
    exit ( 1 );
  }

  return;
}
/******************************************************************************/

void rule_adjust ( double a, double b, double c, double d, int order,
  double x[], double w[] )

/******************************************************************************/
/*
  Purpose:
  
    RULE_ADJUST maps a quadrature rule from [A,B] to [C,D].
  
  Discussion:
  
    Most quadrature rules are defined on a special interval, like
    [-1,1] or [0,1].  To integrate over an interval, the abscissas
    and weights must be adjusted.  This can be done on the fly,
    or by calling this routine.
  
    If the weight function W(X) is not 1, then the weight vector W will
    require further adjustment by the user.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    30 April 2006
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, double A, B, the endpoints of the definition interval.
  
    Input, double C, D, the endpoints of the integration interval.
  
    Input, int ORDER, the number of abscissas and weights.
  
    Input/output, double X[ORDER], W[ORDER], the abscissas
    and weights.
*/
{
  int i;

  for ( i = 0; i < order; i++ )
  {
    x[i] = ( ( b - x[i]     ) * c
           + (     x[i] - a ) * d )
           / ( b        - a );
  }

  for ( i = 0; i < order; i++ )
  {
    w[i] = ( ( d - c ) / ( b - a ) ) * w[i];
  }

  return;
}
/******************************************************************************/

int s_eqi ( char *s1, char *s2 )

/******************************************************************************/
/*
  Purpose:
  
    S_EQI reports whether two strings are equal, ignoring case.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    05 May 2003
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, char *S1, char *S2, pointers to two strings.
  
    Output, bool S_EQI, is true if the strings are equal.
*/
{
  int i;
  int nchar;
  int nchar1;
  int nchar2;

  nchar1 = strlen ( s1 );
  nchar2 = strlen ( s2 );
  if ( nchar1 < nchar2 )
  {
    nchar = nchar1;
  }
  else
  {
    nchar = nchar2;
  }
/*
  The strings are not equal if they differ over their common length.
*/
  for ( i = 0; i < nchar; i++ )
  {

    if ( ch_cap ( s1[i] ) != ch_cap ( s2[i] ) )
    {
      return 0;
    }
  }
/*
  The strings are not equal if the longer one includes nonblanks
  in the tail.
*/
  if ( nchar < nchar1 )
  {
    for ( i = nchar; i < nchar1; i++ )
    {
      if ( s1[i] != ' ' )
      {
        return 0;
      }
    }
  }
  else if ( nchar < nchar2 )
  {
    for ( i = nchar; i < nchar2; i++ )
    {
      if ( s2[i] != ' ' )
      {
        return 0;
      }
    }
  }

  return 1;
}
/******************************************************************************/

double sum_sub ( double func ( double x ), double a, double b, int nsub,
  int order, double xlo, double xhi, double xtab[], double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    SUM_SUB carries out a composite quadrature rule.
  
  Discussion:
  
    SUM_SUB assumes the original rule was written for [XLO,XHI].
  
    The integration interval is [ A, B ].
  
    The integral to approximate:
  
      Integral ( A <= X <= B ) F(X) dX
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    02 May 2006
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, double FUNC ( double X ), the name of the function which
    evaluates the integrand.
  
    Input, double A, B, the lower and upper limits of integration.
  
    Input, int NSUB, the number of equal subintervals into
    which the finite interval (A,B) is to be subdivided for
    higher accuracy.  NSUB must be at least 1.
  
    Input, int ORDER, the order of the rule.
    ORDER must be at least 1.
  
    Input, double XLO, XHI, the left and right endpoints of
    the interval over which the quadrature rule was defined.
  
    Input, double XTAB[ORDER], the abscissas of a quadrature
    rule for the interval [XLO,XHI].
  
    Input, double WEIGHT[ORDER], the weights of the
    quadrature rule.
  
    Output, double SUM_SUB, the approximate value of the integral.
*/
{
  double a_sub;
  double b_sub;
  int i;
  int j;
  double quad_sub;
  double result;
  double result_sub;
  double x;
  double volume;
  double volume_sub;

  if ( a == b )
  {
    result = 0.0;
    return result;
  }

  if ( order < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SUM_SUB - Fatal error!\n" );
    fprintf ( stderr, "  Nonpositive value of ORDER = %d\n", order );
    exit ( 1 );
  }

  if ( nsub < 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SUM_SUB - Fatal error!\n" );
    fprintf ( stderr, "  Nonpositive value of NSUB = %d\n", nsub );
    exit ( 1 );
  }

  if ( xlo == xhi )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SUM_SUB - Fatal error!\n" );
    fprintf ( stderr, "  XLO = XHI.\n" );
    exit ( 1 );
  }

  volume = 0.0;
  result = 0.0;

  for ( j = 1; j <= nsub; j++ )
  {
    a_sub = ( ( double ) ( nsub - j + 1 ) * a
            + ( double ) (        j - 1 ) * b )
            / ( double ) ( nsub         );

    b_sub = ( ( double ) ( nsub - j )     * a
            + ( double ) (        j )     * b )
            / ( double ) ( nsub     );

    quad_sub = 0.0;
    for ( i = 0; i < order; i++ )
    {
      x = ( ( xhi - xtab[i]       ) * a_sub
          + (       xtab[i] - xlo ) * b_sub )
          / ( xhi           - xlo );
      quad_sub = quad_sub + weight[i] * func ( x );
    }

    volume_sub = ( b - a ) / ( ( xhi - xlo ) * ( double ) ( nsub ) );
    result_sub = quad_sub * volume_sub;

    volume = volume + volume_sub;
    result = result + result_sub;
  }

  return result;
}
/******************************************************************************/

double summer ( double func ( double x ), int order, double xtab[],
  double weight[] )

/******************************************************************************/
/*
  Purpose:
  
    SUMMER carries out a quadrature rule over a single interval.
  
  Formula:
  
    RESULT = sum ( 1 <= I <= ORDER ) WEIGHT(I) * FUNC ( XTAB(I) )
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    28 April 2006
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, double FUNC ( double X ), the name of the function which
    evaluates the integrand.
  
    Input, int ORDER, the order of the rule.
  
    Input, double XTAB[ORDER], the abscissas of the rule.
  
    Input, double WEIGHT[ORDER], the weights of the rule.
  
    Output, double SUMMER, the approximate value of the integral.
*/
{
  int i;
  double result;

  result = 0.0;
  for ( i = 0; i < order; i++ )
  {
    result = result + weight[i] * func ( xtab[i] );
  }

  return result;
}
/******************************************************************************/

void summer_gk ( double func ( double x ), int orderg, double weightg[],
  double *resultg, int orderk, double xtabk[], double weightk[],
  double *resultk )

/******************************************************************************/
/*
  Purpose:
  
    SUMMER_GK carries out Gauss-Kronrod quadrature over a single interval.
  
  Discussion:
  
    The abscissas for the Gauss-Legendre rule of order ORDERG are
    not required, since they are assumed to be the even-indexed
    entries of the corresponding Kronrod rule.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    02 May 2006
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, double FUNC ( double x ), the name of the function which
    evaluates the integrand.
  
    Input, int ORDERG, the order of the Gauss-Legendre rule.
  
    Input, double WEIGHTG[ORDERG], the weights of the
    Gauss-Legendre rule.
  
    Output, double *RESULTG, the approximate value of the
    integral, based on the Gauss-Legendre rule.
  
    Input, int ORDERK, the order of the Kronrod rule.  ORDERK
    must equal 2 * ORDERG + 1.
  
    Input, double XTABK[ORDERK], the abscissas of the Kronrod rule.
  
    Input, double WEIGHTK[ORDERK], the weights of the Kronrod rule.
  
    Output, double *RESULTK, the approximate value of the integral,
    based on the Kronrod rule.
*/
{
  double fk;
  int i;

  if ( orderk != 2 * orderg + 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SUMMER_GK - Fatal error!\n" );
    fprintf ( stderr, "  ORDERK must equal 2 * ORDERG + 1.\n" );
    fprintf ( stderr, "  The input value was ORDERG = %d\n", orderg );
    fprintf ( stderr, "  The input value was ORDERK = %d\n", orderk );
    exit ( 1 );
  }

  *resultg = 0.0;
  *resultk = 0.0;

  for ( i = 0; i < orderk; i++ )
  {
    fk = func ( xtabk[i] );

    *resultk = *resultk + weightk[i] * fk;

    if ( ( i % 2 ) == 1 )
    {
      *resultg = *resultg + weightg[(i-1)/2] * fk;
    }
  }
  return;
}
/******************************************************************************/

void sum_sub_gk ( double func ( double x ), double a, double b, int nsub,
  int orderg, double weightg[], double *resultg, int orderk, double xtabk[],
  double weightk[], double *resultk, double *error )

/******************************************************************************/
/*
  Purpose:
  
    SUM_SUB_GK carries out a composite Gauss-Kronrod rule.
  
  Discussion:
  
    The integration interval is [ A, B ].
  
    The integral to approximate:
  
      Integral ( A <= X <= B ) F(X) dX
  
    The quadrature rule:
  
      H = ( B - A ) / NSUB
      XMID(J) = A + 0.5 * H * ( 2 * J - 1 )
  
      Sum ( 1 <= J <= NSUB )
        Sum ( 1 <= I <= ORDERK )
          WEIGHTK(I) * F ( XMID(J) + 0.5 * H * XTABK(I) )
  
    The Gauss-Legendre weights should be computed by LEGCOM or LEGSET.
    The Kronrod abscissas and weights should be computed by KRONSET.
  
    The orders of the Gauss-Legendre and Kronrod rules must satisfy
    ORDERK = 2 * ORDERG + 1.
  
    The Kronrod rule uses the abscissas of the Gauss-Legendre rule,
    plus more points, resulting in an efficient and higher order estimate.
  
    The difference between the Gauss-Legendre and Kronrod estimates
    is taken as an estimate of the error in the approximation to the
    integral.
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    02 May 2006
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    Input, double FUNC ( double X ), the name of the function which
    evaluates the integrand.
  
    Input, double A, B, the lower and upper limits of integration.
  
    Input, int NSUB, the number of equal subintervals into
    which the finite interval (A,B) is to be subdivided for
    higher accuracy.  NSUB must be at least 1.
  
    Input, int ORDERG, the order of the Gauss-Legendre rule.
    ORDERG must be at least 1.
  
    Input, double WEIGHTG[ORDERG], the weights of the
    Gauss-Legendre rule.
  
    Output, double *RESULTG, the approximate value of the
    integral based on the Gauss-Legendre rule.
  
    Input, int ORDERK, the order of the Kronrod rule.
    ORDERK must be at least 1.
  
    Input, double XTABK[ORDERK], the abscissas of the
    Kronrod rule.
  
    Input, double WEIGHTK[ORDERK], the weights of the
    Kronrod rule.
  
    Output, double *RESULTK, the approximate value of the
    integral based on the Kronrod rule.
  
    Output, double *ERROR, an estimate of the approximation
    error.  This is computed by taking the sum of the absolute values of
    the differences between the Gauss-Legendre and Kronrod rules
    over each subinterval.  This is usually a good estimate of
    the error in the value RESULTG.  The error in the Kronrod
    estimate RESULTK is usually much smaller.
*/
{
  double fk;
  double h;
  int i;
  int j;
  double partg;
  double partk;
  double x;
  double xmid;

  *resultg = 0.0;
  *resultk = 0.0;
  *error = 0.0;

  if ( a == b )
  {
    return;
  }

  if ( orderk != 2 * orderg + 1 )
  {
    fprintf ( stderr, "\n" );
    fprintf ( stderr, "SUM_SUB_GK - Fatal error!\n" );
    fprintf ( stderr, "  ORDERK must equal 2 * ORDERG + 1.\n" );
    fprintf ( stderr, "  The input value was ORDERG = %d\n", orderg );
    fprintf ( stderr, "  The input value was ORDERK = %d\n", orderk );
    exit ( 1 );
  }

  h = ( b - a ) / ( double ) ( nsub );

  for ( j = 0; j < nsub; j++ )
  {
    xmid = a + 0.5 * h * ( double ) ( 2 * j + 1 );

    partg = 0.0;
    partk = 0.0;

    for ( i = 0; i < orderk; i++ )
    {
      x = xmid + 0.5 * h * xtabk[i];
      fk = func ( x );
      partk = partk + 0.5 * h * weightk[i] * fk;

      if ( ( i % 2 ) == 1 )
      {
        partg = partg + 0.5 * h * weightg[(i-1)/2] * fk;
      }
    }
    *resultg = *resultg + partg;
    *resultk = *resultk + partk;
    *error = *error + r8_abs ( partk - partg );
  }
  return;
}
/******************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:
  
    TIMESTAMP prints the current YMDHMS date as a time stamp.
  
  Example:
  
    31 May 2001 09:45:54 AM
  
  Licensing:
  
    This code is distributed under the GNU LGPL license.
  
  Modified:
  
    24 September 2003
  
  Author:
  
    John Burkardt
  
  Parameters:
  
    None
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
