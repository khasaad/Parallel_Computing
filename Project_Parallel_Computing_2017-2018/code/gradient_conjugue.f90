MODULE gradient_conjugue
  
  IMPLICIT NONE
  
CONTAINS
  !! Gradient Conjugue pour une iteration en temps
  SUBROUTINE gc(Nx,Ny,U,W,b,Nl,epsilon,coeff_a,coeff_bx,coeff_by,i1,iN,ideb,ifin,Statinfo,Me,Np)
    IMPLICIT NONE
    include 'mpif.h'
    
    INTEGER, Intent(IN)       :: Nx, Ny, Nl, i1, iN, ideb, ifin, Statinfo, Me, Np
    REAL(kind=8), Intent(IN)  :: epsilon, coeff_a, coeff_bx, coeff_by
    REAL(kind=8), DIMENSION(i1:iN), Intent(IN)    :: b
    REAL(kind=8), DIMENSION(ideb:ifin), Intent(INOUT) :: U
    REAL(kind=8), DIMENSION(i1:iN), Intent(INOUT) :: W
    REAL(kind=8), DIMENSION(:), ALLOCATABLE   :: kappa, r, d
    REAL(kind=8)  :: residu,drl,dwl,alpha,betal,beta,residuloc,dr,dw,residul
    INTEGER :: i,l
    
    ALLOCATE(kappa(Nx*Ny),r(Nx*Ny))  
    
    If (Me == 0) THEN
       ALLOCATE(d(i1:iN + Ny))
    ELSE IF (Me == Np-1) THEN
       ALLOCATE(d(i1-Ny:iN))
    Else
       ALLOCATE(d(i1-Ny:iN+ Ny))
    END IF
    
    ! Initialisation Gradient conjugue
    DO i=i1,iN
       kappa(i) = U(i) 
       r(i)     = W(i) - b(i)
       d(i)     = W(i) - b(i)
    END DO
    
    residuloc = SUM(r*r)
    CALL MPI_ALLREDUCE(residuloc,residu,1,MPI_REAL8,&
         & MPI_SUM, MPI_COMM_WORLD, Statinfo)   


    ! boucle du Gradient conjugue
    l=1
    DO WHILE (l<=Nl .and. ( SQRT(residu) .ge. epsilon ))

       ! Calcul de W=Ad
       Call prodMatVec_W(W, d, coeff_a, coeff_bx, coeff_by, Nx, Ny, i1, iN, ideb, ifin, Statinfo, Me, Np)
              
       drl = 0.0d0
       dwl = 0.0d0
       DO i = i1, iN
          drl = drl + d(i)*r(i)
          dwl = dwl + d(i)*W(i)
       ENDDO
        
       CALL MPI_ALLREDUCE(drl,dr,1,MPI_REAL8,&
            & MPI_SUM, MPI_COMM_WORLD, Statinfo)
       CALL MPI_ALLREDUCE(dwl,dw,1,MPI_REAL8,&
            & MPI_SUM, MPI_COMM_WORLD, Statinfo) 
       
       alpha = dr/dw
       
       DO i = i1, iN
          kappa(i) = kappa(i) - alpha*d(i)
          r(i) = r(i) - alpha*W(i)
       END DO

       residuloc = SUM(r*r)

       CALL MPI_ALLREDUCE(residuloc,residul,1,MPI_REAL8,&
         & MPI_SUM, MPI_COMM_WORLD, Statinfo)        
       beta = residul/residu 
       
       DO i = i1, iN
          d(i) = r(i) + beta*d(i)   
       ENDDO
       residu = residul

       l=l+1
       
       
    ENDDO

    DO i = i1, iN
       U(i)=kappa(i)  
    ENDDO

    DEALLOCATE(kappa,r,d) 
  END SUBROUTINE gc


  !! Calcul du produit A*d
  SUBROUTINE prodMatVec_W(W, vec, a, bx, by, Nx, Ny, i1, iN, ideb, ifin, Statinfo, Me, Np)
    
    IMPLICIT NONE
    include 'mpif.h'
    
    INTEGER, Intent(IN)      :: Nx, Ny, i1, iN, ideb, ifin, Me, Np
    REAL(kind=8), Intent(IN) :: a, bx, by
    REAL(kind=8), DIMENSION(ideb:ifin), Intent(INOUT)    :: vec
    REAL(kind=8), DIMENSION(i1:iN), Intent(INOUT) :: W
    INTEGER, Intent(IN) :: Statinfo
    INTEGER, DIMENSION (MPI_STATUS_SIZE)      :: Status
    INTEGER :: q, p, i, MsgTag
    MsgTag = 100
    
    !! Echange de message entre les processus

    If (Np>1) THEN
       If (Me == 0) THEN
          CALL MPI_SEND(vec(iN-Ny+1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Statinfo)
       ELSE IF (Me == Np-1) THEN
          CALL MPI_SEND(vec(i1),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Statinfo)   
       Else
          CALL MPI_SEND(vec(iN-Ny+1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Statinfo)   
          CALL MPI_SEND(vec(i1),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Statinfo)
       END IF

       If (Me == 0) THEN
          CALL MPI_RECV(vec(iN+1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
       ELSE IF (Me == Np-1) THEN
          CALL MPI_RECV(vec(i1-Ny),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)   
       Else
          CALL MPI_RECV(vec(iN+1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
          CALL MPI_RECV(vec(i1-Ny),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
       END IF
    END IF

    
    
    DO i= i1, iN
       CALL indices(i, q ,p, Ny)
       
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
          !!print*,"me ",Me," ",i,p,q
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
  
END MODULE gradient_conjugue
