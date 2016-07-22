!Numeric solution to non-linear equations using Newton and bisection methods
!
!(c) 2016 Matteo Paolieri
!LICENSE: The MIT License (MIT)

module functions 

    contains

    double precision function es1(x) !function to solve
        double precision :: x
        es1 = EXP(x)-(x+1)**(x-1) 
    end
    
    double precision function es1d(x) !derivative of es1
        double precision :: x
        es1d = EXP(x)-(x+1)**(x-2)*(x+(x+1)*log(x+1)-1)
    end
    
end module functions

module metodi

    contains

    subroutine bisection (a,b,f,tolx,x,steps) !bisection method
    
        implicit double precision (A-H, O-Z) 

        interface
            double precision function f(x)
                implicit double precision (A-H, O-Z)
            end function f
        end interface
    
        fa = f(a)
        fb = f(b)
        x = (a+b)/2.0
        fx = f(x)
        
        c = log(b-a)/log(2.0) !log_2(b-a)
        d = log(tolx)/log(2.0) !log_2(tolx)
        kmax = ceiling(real(c-d)) 

        steps = kmax
        
        do k = 2, kmax
        
!~      f1x = |f'(x)| (using incremental ratio) 
        
            f1x = abs((fb-fa)/(b-a))
            if (abs(fx) .le. tolx*f1x) then !arrest criteria
                steps = k-1
                exit 
            else if (fa*fx > 0) then 
                b = x 
                fb = fx
            else 
                a = x
                fa = fx
            end if
            x = (a+b)/2
            fx = f(x)
        end do
        
!~        arrest criteria: |f(x_i)| <= tolx*|f'(xi)|

    end subroutine bisection
    
    subroutine newton(x0, f, f1, tolx, kmax, x, steps) 

      implicit double precision (A-H, O-Z)
        interface
            double precision function f(x)
                implicit double precision  (A-H, O-Z)
            end function f
            
            double precision function f1(x)
                implicit double precision  (A-H, O-Z)
            end function f1
        end interface 
        
        fx = f(x0)
        f1x = f1(x0)
        x = x0 - fx/f1x ! first step of Newton's method
        k = 0
        
        do while ((k<kmax) .and. ( abs(x-x0)>tolx)) ! arrest criteria = |x(k+1)-x(k)|<tolx
            k = k+1
            x0 = x  
            fx = f(x0) 
            f1x = f1(x0)
            x = x-fx/f1x
        end do
    
        if (abs(x-x0)>tolx) then
            print *, "Newton: the method does not converge."
        end if
        
        steps = k
    
    end subroutine newton
    
end module metodi

program nonlin_eq
    
    use metodi
    use functions
    
    implicit double precision (A-H, O-Z)

        a = 2.0
        b = 4.5
        tolx = 1E-10
        
        call bisection(a, b, es1, tolx, xb, stepsb) !xb and stepsb are output
        
        print *, "______ BISECTION ______"
        print '(a, f10.5)', "x = ", xb
        print '(a, f10.0)', "steps = ", stepsb
        
        x0 = 3.0
        kmax = 100
        
        call newton(x0, es1, es1d, tolx, kmax, xn, stepsn) !xn and stepsn are output
        
        print *, "______ NEWTON ______"
        print '(a, f10.5)', "x = ", xn
        print '(a, f10.0)', "steps = ", stepsn
    
end program nonlin_eq
