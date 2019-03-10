PROGRAM Sequentiel2

  use functions2
  use grad_conj2

  IMPLICIT NONE
  Include 'mpif.h'

  !!Déclaration des variables
  INTEGER      :: unit,Ifreq
  INTEGER      :: Nx,Ny,Maxiter, Nl, kt,i,q,p,N
  INTEGER      :: userChoice, nrec, it, cpt
  INTEGER      :: i1, iN, l1, lN, lMax, Me, Statinfo, Np, ideb, ifin
  INTEGER      :: myrank, ok, MsgTag
  REAL(KIND=8) :: Tmax,dt,dx,dy,CFL,t,Lx,Ly,Diff,erreur,erreurloc,time,tl
  REAL(KIND=8) :: coeff_a, coeff_bx, coeff_by,coef_mixte
  REAL(KIND=8) :: alpha, beta,coef_mixte_2mb
  REAL(KIND=8) :: epsilon, SchwartzCD, CD_max, CDl,CDl2
  REAL(KIND=8) :: tstart,tend
  CHARACTER(LEN=14) :: nameSol, nameSolBis
  INTEGER, DIMENSION (MPI_STATUS_SIZE)    :: Status
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: U, g, test, vecTemp_g, vecTemp_d,K_d,K_g
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: b, W, msg_d, msg_g, A, A_diag
  Integer, DIMENSION(:,:), ALLOCATABLE    :: A_ij

  Call MPI_INIT(Statinfo)
  MsgTag = 100
  tstart= MPI_WTIME()
  Call MPI_COMM_RANK(MPI_COMM_WORLD, Me, Statinfo)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, Np, Statinfo)

  !! Lecture des donnees
  unit = 10
  OPEN(unit,file='data2.txt', form='formatted', status='old')
  READ(unit,*) userChoice
  READ(unit,*) Nx,Ny
  READ(unit,*) Lx,Ly,Diff
  READ(unit,*) Tmax
  READ(unit,*) CFL
  READ(unit,*) epsilon,Nl,Ifreq
  READ(unit,*) nrec
  READ(unit,*) alpha,beta
  CLOSE(unit)
  print*, 'alpha', alpha, 'beta', beta

  !Vérification de la taille de la matrice par rapport au nombre de procs
  IF (Np>Nx .OR. Np>Ny) THEN
     PRINT *,' '
     PRINT *,'----Taille du maillage incompatible avec le nombre de procs----'
     PRINT *,' '
     STOP
  END IF

  !! Verification du choix du probleme
  IF (userChoice <1 .OR. userChoice>3) THEN
     PRINT *,' '
     PRINT *,'----Choix du probleme invalide----'
     PRINT *,'Il n''y a que trois problemes resolvable par ce programme, numerotes de 1 à 3'
     PRINT *,' '
     STOP
  END IF

  CALL charge(Me, Np, Nx, i1, iN, l1, lN)
  IF (Me == 0) THEN
     lMax=lN-l1

     !Vérification de la taille de recouvrement
     IF (nrec>lMax .OR. nrec<1) THEN
        PRINT *,' '
        PRINT *,'----Taille de recouvrement incorrecte----'
        PRINT *,' '
        STOP
     END IF
  END IF

  if ( Me==0) then
    N=(lN-l1+2)*Ny
  else if ( Me==Np-1 ) then
    l1=l1-nrec
    N=(lN-l1+2)*Ny
  else
    l1=l1-nrec
    N=(lN-l1+3)*Ny
  end if

  !! Allocation dynamique de la memoire
  allocate(A(N*5), A_ij(2,N*5), A_diag(N), b(N), U(N))
  allocate(msg_d(Ny), msg_g(Ny), vecTemp_d(Ny), vecTemp_g(Ny), K_d(Ny), K_g(Ny))


  !! Determination des pas d'espaces et de temps
  dx = Lx/(Nx+1)
  dy = Ly/(Ny+1)
  If (dx>dy) Then
     dt = CFL*dx*dx/2.
  else
     dt = CFL*dy*dy/2.
  end if


  !! Initialisation des données
  U  = 0.0d0
  !g  = 0.0d0
  !W  = 0.0d0
  t  = 0.0d0
  SchwartzCD = 1.d-6
  cpt=0

  !! Mise en mémoire des trois coefficients de la
  !! matrice A du système linéiare A*U=B à résoudre
  coeff_a = 1.+2.*Diff*dt*(1./(dx*dx)+1./(dy*dy))
  coeff_bx = Diff*dt/(dx*dx)
  coeff_by = Diff*dt/(dy*dy)


  !! Determination du nombre d'itérations a effectuer
  !! pour calculer la solution au temps final Tmax
  Maxiter= Int(Tmax/dt)
  IF( (Tmax - Maxiter*dt) > 1.d-6 ) Maxiter=Maxiter+1

  !! Boucle en temps
  DO kt = 1, Maxiter
     !!dt = MIN (dt, Tmax-t)
     CD_max = 1
     it=0

     DO WHILE (CD_max > SchwartzCD .and. it < 10000)
        !! Calcul de W=AU
        !CALL prodMatVec_W(W, U, coeff_a, coeff_bx, coeff_by, lN-l1+1, Ny, Me,coef_mixte,Np)

        !! Echange de message entre proc
        IF (Np>1) THEN
           IF (Me == 0) THEN
              K_d(1:Ny)=(beta+alpha/dx)*U(N-(nrec+2)*Ny+1:N-(nrec+1)*Ny)-alpha/dx*U(N-(nrec+1)*Ny+1:N-nrec*Ny)
              CALL MPI_SEND(K_d(1:Ny),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Statinfo)
              vecTemp_d(1:Ny) = U(N-(nrec+2)*Ny+1:N-(nrec+1)*Ny)
              CALL MPI_RECV(Msg_d(1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
           ELSE IF (Me == Np-1) THEN
             CALL MPI_RECV(Msg_g(1),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
              K_g(1:Ny)=(beta+alpha/dx)*U(1+(nrec+1)*Ny:(nrec+2)*Ny)-alpha/dx*U(1+nrec*Ny:(nrec+1)*Ny)
              CALL MPI_SEND(K_g(1:Ny),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Statinfo)
              vecTemp_g(1:Ny) = U(1+(nrec+1)*Ny:(nrec+2)*Ny)

           ELSE
             CALL MPI_RECV(Msg_g(1),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)

              K_d(1:Ny)=(beta+alpha/dx)*U(N-(nrec+2)*Ny+1:N-(nrec+1)*Ny)-alpha/dx*U(N-(nrec+1)*Ny+1:N-nrec*Ny)
              K_g(1:Ny)=(beta+alpha/dx)*U(1+(nrec+1)*Ny:(nrec+2)*Ny)-alpha/dx*U(1+nrec*Ny:(nrec+1)*Ny)
              CALL MPI_SEND(K_d(1:Ny),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Statinfo)
              CALL MPI_RECV(Msg_d(1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
              CALL MPI_SEND(K_g(1:Ny),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Statinfo)
              vecTemp_d(1:Ny) = U(N-(nrec+2)*Ny+1:N-(nrec+1)*Ny)
              vecTemp_g(1:Ny) = U(1+(nrec+1)*Ny:(nrec+2)*Ny)
           END IF

        END IF
        !! Calcul du second membre b
        call SecondMembre(b,msg_d,msg_g,U,N,Ny,dt,t,dx,dy,Lx,Ly,coeff_bx,coeff_by,userChoice,me,Np,l1-1)
        call matrice(A,A_ij,A_diag,cpt,N,lN-l1+1,Ny,Me,Np,alpha,beta,dx,coeff_a,coeff_bx,coeff_by)
        call preconditionning(A,A_ij,b,N,cpt,A_diag)
        call BiCGstab(A,A_ij,b,U,N,cpt)
        !CALL gc(lN-l1+1,Ny,U,W,b,Nl,epsilon,coeff_a,coeff_bx,coeff_by, Me,coef_mixte,Np)

        !! Verification de la convergence des boucles de Schwartz
        IF (Np>1) THEN
           CDl  = 0
           IF (Me == 0) THEN
              CDl = 0.
              DO i=1,Ny
                 CDl = CDl + (vecTemp_d(i)-U(N-(nrec+2)*Ny+i))*(vecTemp_d(i)-U(iN-(nrec+2)*Ny+i))
              END DO
           ELSE IF (Me == Np-1) THEN
              DO i=1,Ny
                 CDl = CDl + (vecTemp_g(i)-U(i+(nrec+1)*Ny))*(vecTemp_g(i)-U(i+(nrec+1)*Ny))
              END DO
           ELSE
              CDl2 = 0
              DO i=1,Ny
                 CDl  = CDl + (vecTemp_d(i)-U(N-(nrec+2)*Ny+i))*(vecTemp_d(i)-U(N-(nrec+2)*Ny+i))
                 CDl2 = CDl2 + (vecTemp_g(i)-U(i+(nrec+1)*Ny))*(vecTemp_g(i)-U(i+(nrec+1)*Ny))
              END DO
              IF (CDl2 > CDl) THEN
                 CDl = CDl2
              END IF
           END IF

           CALL MPI_ALLREDUCE(CDl,CD_max,1,MPI_REAL8, MPI_MAX,MPI_COMM_WORLD, Statinfo)
           CD_max = sqrt(CD_max*dx)
        ELSE
           CD_max = 1.d-7
        END IF

        it = it+1

        ! IF ((Me==0)) THEN
        !   PRINT*,"Boucle", kt,"iterations :", it,"err",CD_max
        ! END IF
     END DO

     IF ((Me==0)) THEN
       PRINT*,"Boucle", kt,"; nombre d'iterations :", it,"err",CD_max
     END IF

     t = t + dt;
  ENDDO

  IF ((Me==0)) THEN
     PRINT*,"t finale", t
  END IF

  !!Calcul de l'erreur seulement pour les choix 1 et 2
  IF (userChoice == 1 .OR. userChoice == 2) THEN
     erreurloc=0.
     erreur=0.
     DO i=1,N
        CALL indices(i, q, p, Ny)
        erreurloc = erreurloc +(U(i)-solExacte((q+l1-1)*dx,p*dy,t,userChoice))**2.
     END DO
     CALL MPI_ALLREDUCE(erreurloc,erreur,1,MPI_REAL8,&
         & MPI_SUM, MPI_COMM_WORLD, Statinfo)
     erreur = sqrt(erreur*dx)
     IF (Me == 0) THEN
       PRINT *,"L'erreur est : ", erreur
     END IF
  END IF

   !! Ecriture des resultats
  unit = 11

  !!Appel de la subroutine rename pour renomer le fichier résultat selon le proc
  CALL Rename(Me,nameSol, nameSolBis)
  OPEN(unit,file=nameSol, form='formatted', status='unknown')
  DO i = 1,N
     CALL indices(i, q, p, Ny)
     WRITE(unit,*) (q+l1-1)*dx, p*dy, U(i)
  ENDDO
  CLOSE(unit)
  !unit = 12
  !OPEN(unit,file=nameSolBis, form='formatted', status='unknown')
  !DO i = ideb,ifin
  !   CALL indices(i, q, p, Ny)
  !   WRITE(unit,*) q*dx, p*dy, U(i)
  !ENDDO
  !CLOSE(unit)
  !print*, 'namesolbis'
  !!DEALLOCATE(U,g,b,W,test,msg_d,msg_g,vecTemp_d,vecTemp_g,K_d,K_g)

  deallocate(A, A_ij, A_diag, b, U)
  deallocate(msg_d, msg_g, vecTemp_d, vecTemp_g, K_d, K_g)
  tend=MPI_WTIME()
  CALL MPI_ALLREDUCE(tend-tstart,time,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,Statinfo)
  IF (Me==0) THEN
    PRINT*, "Le temps d'execution est :", time
  END IF

  !!Ecriture des differents temps de calcul
  !!Pour le calcul du speedup et de l'efficacité
  !!Sauvegarde du temps d'execution sur 1 processeur
  !!dans le fichier "proc"
  IF (Me==0) then
     OPEN(unit,file='proc', form='formatted', status='unknown')
     IF (Np==1) THEN
        WRITE(unit,*) time
        tl=time
     END IF
     CLOSE(unit)

     !!Programme executé avec plus d'un processeur :
     !!Recuperation du temps d'exectution sur 1 proc dans le fichier "proc"
     !!Puis ecriture dans le fichier "temps"
     !!du nombre de proc, du speedup et de l'efficacité
     !! ATTENTION : Il peut y avoir une erreur à l'execution
     !!             si le fichier proc a été supprimé et que
     !!             le programme est directement executé sur 1 processeur
     !!             mais les fichiers résultas sont quand même crées
     OPEN(unit,file='proc', form='formatted', status='unknown',iostat = ok)
     IF (Np/=0) THEN
        READ(unit,*,iostat=ok)tl
     END IF
     CLOSE(unit)

     OPEN(unit,file='temps', form='formatted', status='old', position='append')
     WRITE(unit,*) Np, tl/(Np*time),tl/time
     CLOSE(unit)
  END IF

  CALL MPI_FINALIZE(Statinfo)

CONTAINS
  Subroutine Charge(Me, Np, N, i1, iN, l1, lN)
    INTEGER, Intent (IN)  :: Me, Np, N
    INTEGER, Intent (OUT) :: i1, iN, l1, lN
    INTEGER               :: r, q

    q = N/Np
    r = N-q*Np

    IF (Me <r) THEN
       l1=Me*(q+1)+1
       lN=(Me+1)*(q+1)
    ELSE
       l1=1+r+Me*q
       lN=l1+q-1
    END IF

    i1 = (l1-1)*Ny+1
    iN = lN*Ny

  END Subroutine Charge

  Subroutine Rename(Me,name, name2)
    INTEGER, Intent(IN) :: Me
    CHARACTER(LEN=3)    :: tn
    CHARACTER(LEN=14), Intent(OUT) :: name, name2
    INTEGER :: i1,i2,i3
    i1 = Me/100
    i2 =( Me - 100*i1)/10
    i3 = Me - 100*i1 -10*i2
    tn = char(i1+48)//char(i2+48)//char(i3+48)
    name='sol'//tn//'.dat'
    name2='solbis'//tn//'.dat'
  END Subroutine Rename

  !!Calcul second membre
  subroutine SecondMembre(b,msg_d,msg_g,U,N,Ny,dt,t,dx,dy,Lx,Ly,coeff_bx,coeff_by,userChoice,me,Np,l_prec)

    integer, intent(in)                :: me, Np, N, Ny, userChoice,l_prec
    real(kind=8), dimension(:), intent(in)   :: msg_d, msg_g, U
    real(kind=8), dimension(:), intent(out)  :: b
    real(kind=8), intent(in)                 :: dt, dx, dy, Lx, Ly, t, coeff_bx, coeff_by
    integer                            :: i, p, q


    if ( me==0 ) then
      b(N-Ny+1:N)=msg_d(1:Ny)
      DO i=1,N-Ny
        CALL indices(i, q, p, Ny)
        b(i) = dt*f((q+l_prec)*dx,p*dy,t,Lx,Ly,userChoice)+U(i)
        IF (q==1) THEN
          b(i) = b(i) + coeff_bx*funct_h((q+l_prec-1)*dx,p*dy,t,userChoice)
        END IF
        IF (p==1) THEN
          b(i) = b(i) + coeff_by*funct_g((q+l_prec)*dx,(p-1)*dy,t,userChoice)
        END IF
        IF (p==Ny) THEN
          b(i) = b(i) + coeff_by*funct_g((q+l_prec)*dx,(p+1)*dy,t,userChoice)
        END IF
      END DO

    else if ( me==Np-1 ) then
      b(1:Ny)=msg_g(1:Ny)
      DO i=Ny+1,N
         CALL indices(i, q, p, Ny)
         b(i) = dt*f((q+l_prec)*dx,p*dy,t,Lx,Ly,userChoice)+U(i)
         IF (q==N/Ny) THEN
           b(i) = b(i) + coeff_bx*funct_h((q+l_prec+1)*dx,p*dy,t,userChoice)
         END IF
         IF (p==1) THEN
           b(i) = b(i) + coeff_by*funct_g((q+l_prec)*dx,(p-1)*dy,t,userChoice)
         END IF
         IF (p==Ny) THEN
           b(i) = b(i) + coeff_by*funct_g((q+l_prec)*dx,(p+1)*dy,t,userChoice)
         END IF
      END DO

    ELSE
      b(1:Ny)=msg_g(1:Ny)
      b(N-Ny+1:N)=msg_d(1:Ny)
      DO i=Ny+1,N-Ny
         CALL indices(i, q, p, Ny)
         b(i) = dt*f((q+l_prec)*dx,p*dy,t,Lx,Ly,userChoice)+U(i)
         IF (p==1) THEN
           b(i) = b(i) + coeff_by*funct_g((q+l_prec)*dx,(p-1)*dy,t,userChoice)
         END IF
         IF (p==Ny) THEN
           b(i) = b(i) + coeff_by*funct_g((q+l_prec)*dx,(p+1)*dy,t,userChoice)
         END IF

      END DO
    end if
  end subroutine SecondMembre



END PROGRAM Sequentiel2
