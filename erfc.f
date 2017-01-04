	FUNCTION erfc ( x )

c    *******************************************************************
c    ** approximation to the complementary error function             **
c    **                                                               **
c    ** reference:                                                    **
c    **                                                               **
c    ** abramowitz and stegun, handbook of mathematical functions,    **
c    **    national bureau of standards, formula 7.1.26               **
c    *******************************************************************

        double precision        a1, a2, a3, a4, a5, p

        parameter ( a1 = 0.254829592, a2 = -0.284496736 )
        parameter ( a3 = 1.421413741, a4 = -1.453152027 )
        parameter ( a5 = 1.061405429, p  =  0.3275911   )

        double precision        t, x, xsq, tp, erfc

c    *******************************************************************

        t  = 1.0 / ( 1.0 + p * x )
        xsq = x * x

        tp = t * ( a1 + t * ( a2 + t * ( a3 + t * ( a4 + t * a5 ) ) ) )

        erfc = tp * exp ( -xsq )

        return
        end