MODULE grad_conj1

  IMPLICIT NONE

CONTAINS
  !! Gradient Conjugue pour une iteration en temps
  SUBROUTINE gc(Nx,Ny,U,W,b,Nl,epsilon,coeff_a,coeff_bx,coeff_by, Me)
    INTEGER, Intent(IN)       :: Nx, Ny, Nl, Me
    REAL(kind=8), Intent(IN)  :: epsilon, coeff_a, coeff_bx, coeff_by
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
        Call prodMatVec_W(W, d, coeff_a, coeff_bx, coeff_by, Nx, Ny, Me)

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

     DO i=1,Nx*Ny
        U(i)=kappa(i)
     ENDDO

     DEALLOCATE(kappa,r,d)
  END SUBROUTINE gc


  !! Calcul du produit A*d
  SUBROUTINE prodMatVec_W(W, vec, a, bx, by, Nx, Ny, Me)
    INTEGER, Intent(IN)      :: Nx, Ny, Me
    REAL(kind=8), Intent(IN) :: a, bx, by
    REAL(kind=8), DIMENSION(:), Intent(IN)    :: vec
    REAL(kind=8), DIMENSION(:), Intent(INOUT) :: W
    INTEGER :: q, p, i

    !!PRINT*, " "
    !!PRINT*, "Me  Nx  Ny  i  q  p"
    DO i= 1, Nx*Ny
       CALL indices(i, q ,p, Ny)

       !!PRINT*,Me,Nx,Ny,i,q,p

       IF (q==1 .AND. p==1) THEN
          W(i) = a*vec(i) - bx*vec(i+Ny)- by*vec(i+1)
       ELSE IF (q==1 .AND. p==Ny) THEN
          W(i) = a*vec(i) - bx*vec(i+Ny) - by*vec(i-1)
       ELSE IF (q==Nx .AND. p==1) THEN
          W(i) = a*vec(i)- bx*vec(i-Ny)- by*vec(i+1)
       ELSE IF (q==Nx .AND. p==Ny ) THEN
          W(i) = a*vec(i) - bx*vec(i-Ny) - by*vec(i-1)
       ELSE IF (q==1) THEN
          W(i) = a*vec(i) - bx*vec(i+Ny)- by*vec(i+1) - by*vec(i-1)
       ELSE IF (q==Nx) THEN
          W(i) = a*vec(i) - bx*vec(i-Ny)- by*vec(i+1) - by*vec(i-1)
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

END MODULE grad_conj1
