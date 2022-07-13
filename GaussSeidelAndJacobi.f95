program GaussSeidelAndJacobi
    ! Lets call this matrix - A. I was resolving equation of Ax=y
    ! where y is in a shape of  {1 ..... N}^T
    !
    ! |3    1    0.2   0
    ! |1    3     1    0.2
    ! |0.2  1     3    1
    ! |0   0.2    1    3
    ! |0    0    0.2   1

    implicit none
    integer, parameter :: matrix_size = 100
    integer :: i,j,ii

    ! Variables for Jacobi method
    real(kind=selected_real_kind(16)), dimension (matrix_size,matrix_size) :: A_2, M
    real(kind=selected_real_kind(16)), dimension (matrix_size) :: N_2, x_tmp1, x_tmp2
    real(kind=selected_real_kind(16)) :: distance, prev_distance = 100, eps = 1E-8
    real(kind=selected_real_kind(16)) :: max_delta
    ! --------------------------------------------------------------------------

    real(kind=selected_real_kind(16)), dimension (matrix_size) :: AD, x, N, D,xi,xt
    real(kind=selected_real_kind(16)), dimension (matrix_size-1) :: Au1, Al1,u1,l1
    real(kind=selected_real_kind(16)), dimension (matrix_size-2) :: Au2, Al2,u2,l2
    logical :: thesame = .false.
    integer :: iter_counter = 0, iter_counter2 = 0

    ! SETTING UP ARRAYS WHICH REPRESENT DIAGONALS OF CONSIDERED MATRIX
    ! AND VECTORS WHICH ARE PART OF THE EQUATION
    ! OR WILL BE USED FOR COMPUTATION

    do i=1,matrix_size,1
        x(i)=real(i)
        AD(i)=3.0
        N(i)=1.0/AD(i)
        D(i)=x(i)*N(i)

        if(i > 1) then
            Al1(i-1)=1.0
            Au1(i-1)=1.0
            u1(i-1)=Au1(i-1)*N(i-1)
            l1(i-1)=Al1(i-1)*N(i-1)
        end if
        if(i > 2) then
            Al2(i-2)=0.2
            Au2(i-2)=0.2
            u2(i-2)=Au1(i-2)*N(i-2)
            l2(i-2)=Al1(i-2)*N(i-2)
        end if
        xi(i)=0.0
        xt(i)=0.0

    end do

    !--------------------------------------------------------------------
    ! Gauss-Seidel method
    
    ! Creating new files for results
    open(unit=1,file="gauss-seidel.dat", status="new", &
            action="write", position="append")
    open(unit=2,file="result-vector-gauss.dat",status="new", &
            action="write", position="append")
    print *, "Results showing covergence of Gauss-Seidel method in following iterations."
    print *
    print *, "First column is # of itereation, second - relative error:"
    print *

    do while(.true.)
        iter_counter = iter_counter + 1
        xi(1) = (D(1) - (u1(1) * xt(2)) - (u2(1) * xt(3)))
        xi(2) = (D(2) - (u1(2) * xt(3)) - (u2(2) * xt(4)) - (l1(1) * xi(1)))

        do i=1,96
            xi(i+2) = (D(i+2) - (u1(i+2) * xt(i+3)) - &
            (u2(i+2) * xt(i+4)) - (l1(i+1) * xi(i+1)) - (l2(i) * xi(i)))
        end do

        xi(99) = (D(99) - (l1(98) * xi(98)) - (l2(97) * xi(97)))  - (u1(99) * xt(100))
        xi(100) = (D(100) - (l1(99) * xi(99)) - (l2(98) * xi(98)))


        do ii=1,matrix_size
            max_delta = abs(xi(1) - xt(1))
            if(max_delta < abs(xi(ii)-xt(ii))) max_delta=abs(xi(ii)-xt(ii))
            xt(ii) = xi(ii)           
        end do


        ! Writing relative errors to file and printing to std output
        write(1,fmt="(I0,':',2X,E15.9)") iter_counter, max_delta
        write(*,fmt="(I0,':',2X,E15.9)") iter_counter, max_delta

        if(max_delta < eps) exit
    end do

    ! Writing result array to file
    do i=1,matrix_size
        write(2,fmt="(F16.7)") xt(i)
    end do

    close(1)
    close(2)

    !--------------------------------------------------------------
    ! Jacobi method

    ! Setting up matrix and other used variables
 
    iter_counter2 = 0
    do i=1,matrix_size
        do j=1,matrix_size
            if(i.eq.j) then
                A_2(i,j)=3.0
            else if(i-1.eq.j) then
                A_2(i,j)=1.0
            else if(i.eq.j-1) then
                A_2(i,j)=1.0
            else if(i-2.eq.j) then
                A_2(i,j)=0.2
            else if(i.eq.j-2) then
                A_2(i,j)=0.2
            end if
        end do
    end do

    do i=1,matrix_size
        N(i) = 1/A_2(i,i)
        x_tmp1(1) = 0.0
    end do

    do i=1,matrix_size
        do j=1,matrix_size
            if(i.eq.j) then
                M(i,j)=0.0
            else
                M(i,j) = M(i,j) - (A_2(i,j) * N(i))
            end if
        end do
    end do

    ! Creating new files for results
    open(unit=3,file="jacobi.dat", status="new", &
            action="write", position="append")
    open(unit=4,file="result-vector-jacobi.dat",status="new", &
            action="write", position="append")
    
    print *
    print *, "Results showing covergence of Jacobi method in following iterations."
    print *
    print *, "First column is # of itereation, second - relative error:"
    print *

    do while(thesame.neqv..true.)
        iter_counter2 = iter_counter2 + 1
        do i=1,matrix_size
            x_tmp2(i) = N(i) * x(i)
            do j=1,matrix_size
                x_tmp2(i) = x_tmp2(i) + M(i,j) *x_tmp1(j)
            end do
        end do

        do ii=1,matrix_size
            max_delta = abs(x_tmp2(1) - x_tmp1(1))
            if(max_delta < abs(x_tmp2(ii)-x_tmp1(ii))) max_delta=abs(x_tmp2(ii)-x_tmp1(ii))
            x_tmp1(ii) = x_tmp2(ii)
        end do

        ! Writing relative errors to file and printing to std output
        write(3,fmt="(I0,':',2X,E15.9)") iter_counter2, max_delta
        write(*,fmt="(I0,':',2X,E15.9)") iter_counter2, max_delta

        ! exit if desired accuracy is reached
        if(max_delta.le.eps) exit
    end do

    ! Writing result array to file
    do i=1,matrix_size
        write(4,fmt="(F16.7)") x_tmp2(i)
    end do

    close(3)
    close(4)

end program GaussSeidelAndJacobi
