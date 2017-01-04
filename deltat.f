	FUNCTION deltat( tempr )

c    *******************************************************************
c    ** calculates thickness of polymer brush as a function			**
c    ** of temp it gives delta in nm									**
c    *******************************************************************

        double precision deltat, tempr, tstar

        tstar = ( tempr - 5 ) / 80d0
        deltat = ( 4.5 + 87.36 * exp ( -exp ( 1.6 * tstar + 0.15 ) ) )

        end
c
c