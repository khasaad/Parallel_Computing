MODULE grad_conj2

  IMPLICIT NONE

CONTAINS
  !! Gradient Conjugue pour une iteration en temps
  SUBROUTINE gc(Nx,Ny,U,W,b,Nl,epsilon,coeff_a,coeff_bx,coeff_by, Me,coef_mixte,Np)
    INTEGER, Intent(IN)       :: Nx, Ny, Nl, Me,Np
    REAL(kind=8), Intent(IN)  :: epsilon, coeff_a, coeff_bx, coeff_by,coef_mixte
    REAL(kind=8), DIMENSION(:), Intent(IN)    :: b
    REAL(kind=8), DIMENSION(:), Intent(INOUT) :: U,W
    REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: kappa, r, d
    REAL(kind=8)  :: residu,drl,dwl,alpha,betal,beta
    INTEGER :: i,l

    ALLOCATE(kappa(Nx*Ny),r(Nx*Ny),d(Nx*Ny))

     !initialisation Gradient conjugue
     DO i=1,Nx*Ny
        kappa(i) = U(i)
        r(i)     = W(i) - b(i)
        d(i)     = W(i) - b(i)
     END DO

     residu = 0.d0
     DO i=1,Nx*Ny
        residu = residu + r(i)*r(i) !SUM(r*r)
     ENDDO

     ! boucle du Gradient conjugue
     l=1
     DO WHILE (l<=Nl .and. ( SQRT(residu) .ge. epsilon ))

        !Calcul de W=Ad
        Call prodMatVec_W(W, d, coeff_a, coeff_bx, coeff_by, Nx, Ny, Me,coef_mixte,Np)

        drl = 0.0d0
        dwl = 0.0d0
        DO i = 1, Nx*Ny
           drl = drl + d(i)*r(i)
           dwl = dwl + d(i)*w(i)
        ENDDO
        alpha = drl/dwl
        DO i=1,Nx*Ny
           kappa(i) = kappa(i) - alpha*d(i)
           r(i) = r(i) - alpha*W(i)
        END DO
        betal=0.d0
        DO i=1,Nx*Ny
           betal=betal+ r(i)*r(i)
        ENDDO
        ! beta = SUM(r*r)/residu
        beta = betal/residu
        DO i=1,Nx*Ny
           d(i) = r(i) + beta*d(i)
        ENDDO
        residu = 0.d0
        DO i=1,Nx*Ny
           residu    = residu+r(i)*r(i)!SUM(r*r)
        ENDDO
        l=l+1
        !Fin Gradient conjugue
     ENDDO
     print*, 'itérations gradient ', l

     DO i=1,Nx*Ny
        U(i)=kappa(i)
     ENDDO

     DEALLOCATE(kappa,r,d)
  END SUBROUTINE gc


  !! Calcul du produit A*d
  SUBROUTINE prodMatVec_W(W, vec, a, bx, by, Nx, Ny, Me,coef_mixte,Np)
    INTEGER, Intent(IN)      :: Nx, Ny, Me,Np
    REAL(kind=8), Intent(IN) :: a, bx, by,coef_mixte
    REAL(kind=8), DIMENSION(:), Intent(IN)    :: vec
    REAL(kind=8), DIMENSION(:), Intent(INOUT) :: W
    INTEGER :: q, p, i

    !!PRINT*, " "
    !!PRINT*, "Me  Nx  Ny  i  q  p"
    DO i= 1, Nx*Ny
       CALL indices(i, q ,p, Ny)

       !!PRINT*,Me,Nx,Ny,i,q,p

       IF (q==1 .AND. p==1) THEN
         if ( Me==0 ) then
           W(i) = a*vec(i) - bx*vec(i+Ny)- by*vec(i+1)
         ELSE
           W(i) = (a+coef_mixte)*vec(i) - bx*vec(i+Ny)- by*vec(i+1)
         end if
       ELSE IF (q==1 .AND. p==Ny) THEN
         if ( Me==0 ) then
           W(i) = a*vec(i) - bx*vec(i+Ny) - by*vec(i-1)
         else
           W(i) = (a+coef_mixte)*vec(i) - bx*vec(i+Ny) - by*vec(i-1)
         end if
       ELSE IF (q==Nx .AND. p==1) THEN
         if ( Me==Np-1 ) then
           W(i) = a*vec(i)- bx*vec(i-Ny)- by*vec(i+1)
         else
           W(i) = (a+coef_mixte)*vec(i)- bx*vec(i-Ny)- by*vec(i+1)
         end if
       ELSE IF (q==Nx .AND. p==Ny ) THEN
         if ( Me==Np-1 ) then
           W(i) = a*vec(i) - bx*vec(i-Ny) - by*vec(i-1)
         else
           W(i) = (a+coef_mixte)*vec(i) - bx*vec(i-Ny) - by*vec(i-1)
         end if
       ELSE IF (q==1) THEN
         if ( Me==0 ) then
           W(i) = a*vec(i) - bx*vec(i+Ny)- by*vec(i+1) - by*vec(i-1)
         else
           W(i) = (a+coef_mixte)*vec(i) - bx*vec(i+Ny)- by*vec(i+1) - by*vec(i-1)
         end if
       ELSE IF (q==Nx) THEN
         if ( Me==Np-1 ) then
           W(i) = a*vec(i) - bx*vec(i-Ny)- by*vec(i+1) - by*vec(i-1)
         else
           W(i) = (a+coef_mixte)*vec(i) - bx*vec(i-Ny)- by*vec(i+1) - by*vec(i-1)
         end if
       ELSE IF (p==1) THEN
          W(i) = a*vec(i) - bx*vec(i+Ny)- bx*vec(i-Ny)- by*vec(i+1)
       ELSE IF (p==Ny) THEN
          W(i) = a*vec(i) - bx*vec(i+Ny)- bx*vec(i-Ny) - by*vec(i-1)
       ELSE
          W(i) = a*vec(i) - bx*vec(i+Ny)- bx*vec(i-Ny)- by*vec(i+1) - by*vec(i-1)
       END IF
    END DO
  END SUBROUTINE prodMatVec_W


  !! Calcul de la position (q,p) du i-eme terme d'un vecteur
  SUBROUTINE indices(i, q, p, Ny)
    INTEGER, Intent(IN)    :: i, Ny
    INTEGER, Intent(INOUT) :: p, q
    p = MOD(i,Ny)
    q = i/Ny+1
    IF (p==0) THEN
       p=Ny
       q=q-1
    END IF
  END SUBROUTINE indices


  !! Calcul de la matrice
  subroutine matrice(A,A_ij,A_diag,cpt,N,Nx,Ny,Me,Npc,alpha,beta,dx,coeff_a,coeff_bx,coeff_by)
    integer, intent(in)                           :: N,Nx,Ny,Me,Npc
    real(kind=8), dimension(N*5), intent(out)     :: A
    integer, dimension(2,N*5), intent(out)        :: A_ij
    real(kind=8), dimension(N), intent(out)       :: A_diag
    integer, intent(out)                          :: cpt
    real(kind=8)                                  :: alpha, beta, dx, coeff_a, coeff_bx, coeff_by
    integer                                       :: i,i1,iN,j,k

    if ( Me==0 ) then
      i1=1
      iN=Nx*Ny
    ELSE
      i1=1+Ny
      iN=(Nx+1)*Ny
    end if
    cpt=0
    do i = 1, Ny
      if ( Me>0 ) then
        cpt=cpt+1
        A(cpt)=alpha/dx+beta
        A_ij(1,cpt)=i
        A_ij(2,cpt)=i
        A_diag(i)=A(cpt)
        cpt=cpt+1
        A(cpt)=-alpha/dx
        A_ij(1,cpt)=i
        A_ij(2,cpt)=i+Ny
      end if
      if ( Me<Npc-1 ) then
        cpt=cpt+1
        A(cpt)=alpha/dx+beta
        A_ij(1,cpt)=iN+i
        A_ij(2,cpt)=iN+i
        A_diag(iN+i)=A(cpt)
        cpt=cpt+1
        A(cpt)=-alpha/dx
        A_ij(1,cpt)=iN+i
        A_ij(2,cpt)=iN+i-Ny
      end if
    end do

    do i = 1,Nx
      do j = 1,Ny
        k=i1-1+j+(i-1)*Ny
        cpt=cpt+1
        A(cpt)=coeff_a
        A_ij(1,cpt)=k
        A_ij(2,cpt)=k
        A_diag(k)=A(cpt)

        if ( me>0 .or. i>1 ) then
          cpt=cpt+1
          A(cpt)=-coeff_bx
          A_ij(1,cpt)=k
          A_ij(2,cpt)=k-Ny
        end if

        if ( me<npc-1 .or. i<Nx ) then
          cpt=cpt+1
          A(cpt)=-coeff_bx
          A_ij(1,cpt)=k
          A_ij(2,cpt)=k+Ny
        end if

        if ( j>1 ) then
          cpt=cpt+1
          A(cpt)=-coeff_by
          A_ij(1,cpt)=k
          A_ij(2,cpt)=k-1
        end if

        if ( j<Ny ) then
          cpt=cpt+1
          A(cpt)=-coeff_by
          A_ij(1,cpt)=k
          A_ij(2,cpt)=k+1
        end if
      end do
    end do
  end subroutine matrice



  subroutine preconditionning(A,A_ij,Rhs,N,cpt,dFdiag)
    integer , intent(in) :: N ,cpt
    double precision, dimension(1:N) ,intent(inout):: Rhs ,dFdiag
    double precision,dimension(1:N*5),intent(inout) :: A
    integer,dimension(2,1:N*5),intent(in) :: A_ij

    integer :: k

    do k=1,N
      if (abs(dFdiag(k)).gt.1e-8) then
        rhs(k) = rhs(k)/  dFdiag(k)
      end if
    end do

    do k=1,cpt
      if (dFdiag(A_ij(1,k)).gt.1e-8) then
        A(k) = A(k)/  dFdiag(A_ij(1,k)) ! A(i,j) = A(i,j) / A(i,i)
      end if
    end do
  end subroutine preconditionning

  subroutine BiCGstab(A,A_ij,Rhs,X,N,cpt)
    integer , intent(in) :: N ,cpt
    double precision, dimension(1:N) ,intent(in):: Rhs
    double precision,dimension(1:N*5),intent(in) :: A
    integer,dimension(2,1:N*5),intent(in) :: A_ij
    double precision, dimension(1:N) ,intent(inout):: X

    integer :: itMax,i,it
    double precision, dimension(1:N) ::  r0,r0hat,v0,p0,p,s,ti,Prod

    double precision :: rho0,alpha,omega0,ts,tt,rho,rv,eps,conv_abs,conv_rel,beta
    double precision :: normB

    conv_rel=1.d-8
    conv_abs=1.d-10
    itMax=10000

    call ProdMatVec(A,A_ij,X,Prod,N,cpt)

    normB = 0.0d0

    do i=1,N
      r0(i)=Rhs(i)-Prod(i)
      r0hat(i)=r0(i)
      v0(i)=0.d0
      p0(i)=0.d0
      normB=normB+Rhs(i)**2
    end do

    normB = sqrt(normB)

    rho0=1.d0
    alpha=1.d0
    omega0=1.d0

    it=1
    eps=conv_rel*normB+conv_abs+10.d0
    !print*, conv_rel*normB+conv_abs
    rv=1.0
    do while( (it .lt. itMax) .and. (eps .gt.conv_rel*normB+conv_abs  ) .and. (abs(rv)> 1.0e-30))

      rho=0.d0
      do i=1,N
        rho=rho+r0hat(i)*r0(i)
      end do

      beta=(rho/rho0)*(alpha/omega0)

      do i=1,N
        p(i)=r0(i)+beta*(p0(i)-omega0*v0(i))
      end do
      do i=1,5*N
      end do

      call ProdMatVec(A,A_ij,p,v0,N,cpt)

      rv=0.d0
      do i=1,N
        rv=rv+r0hat(i)*v0(i)
      end do
      do i=1,N
      end do

      if (abs(rv)> 1.0e-33) then
        alpha=rho/rv

        do i=1,N
          s(i)=r0(i)-alpha*v0(i)
        end do

        call ProdMatVec(A,A_ij,s,ti,N,cpt)

        ts=0.d0
        tt=0.d0
        do i=1,N
          ts=ts+ti(i)*s(i)
          tt=tt+ti(i)*ti(i)
        end do

        omega0=ts/tt

        eps=0.d0
        do i=1,N
          X(i)=X(i)+alpha*p(i)+omega0*s(i)
          r0(i)=s(i)-omega0*ti(i)
          eps=eps+r0(i)**2
          p0(i)=p(i)
        end do
        rho0=rho

        eps=sqrt(eps)
      end if
      it=it+1

    end do

    if (eps .gt.conv_rel*normB+conv_abs) then
      write(*,'(A,i8,A,ES10.2)') '     Fin bicg , iter =', it,' , eps =',eps
    end if
  end subroutine BiCGstab



  subroutine  ProdMatVec(A,A_ij,X,Prod,N,cpt)
    integer , intent(in) :: N ,cpt
    double precision,dimension(1:N*5),intent(in) :: A
    integer        ,dimension(2,1:N*5),intent(in) :: A_ij
    double precision, dimension(1:N) ,intent(in):: X
    double precision, dimension(1:N) , intent(out) ::Prod

    integer :: i

    Prod(:) = 0.0

    do  i = 1,cpt
      Prod (A_ij(1,i)) = Prod (A_ij(1,i)) +  A(i)*X(A_ij(2,i))
    end do
  end subroutine ProdMatVec


  subroutine Jacobi(A,A_ij,A_diag,N,Me,cpt,X,Rhs,dx)
    integer, intent(in)                           :: N,Me,cpt
    real(kind=8), dimension(N*5), intent(in)     :: A
    integer, dimension(2,N*5), intent(in)        :: A_ij
    real(kind=8), dimension(N), intent(in)       :: A_diag
    double precision, dimension(1:N) ,intent(inout) :: X
    double precision, dimension(1:N) ,intent(in):: Rhs
    real(kind=8), intent(in) :: dx
    double precision, dimension(1:N) :: Z,residu
    integer                                       :: k, it_max, i
    real(kind=8) :: norm_res, epsilon, A_ii

    it_max=10000

    call ProdMatVec(A,A_ij,X,residu,N,cpt)

    residu = Rhs - residu

    norm_res = sqrt(sum(dx*abs(residu**2)))

    k = 0

    epsilon = 1.d-10

    do while ( (norm_res .gt. epsilon).and.(k .lt. it_max) )

      Z = Rhs

      do k = 1,cpt

        if ( A_ij(1,k) .eq. A_ij(2,k) ) then

          Z(A_ij(1,k)) = Z(A_ij(1,k))

        else

          Z(A_ij(1,k)) = Z(A_ij(1,k)) - A(k)*X(A_ij(2,k))

        endif

      ENDDO

      do k = 1,N

        Z(k) = Z(k) / A_diag(k)

      ENDDO

      X = Z

      call ProdMatVec(A,A_ij,X,residu,N,cpt)

      residu = Rhs - residu

      norm_res = sqrt(sum(dx*abs(residu**2)))

      k = k+1

    ENDDO

end subroutine Jacobi




END MODULE grad_conj2
