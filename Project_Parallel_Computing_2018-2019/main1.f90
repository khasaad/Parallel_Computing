PROGRAM Sequentiel1

  use functions1
  use grad_conj1

  IMPLICIT NONE
  Include 'mpif.h'

  !!Déclaration des variables
  INTEGER      :: unit,Ifreq
  INTEGER      :: Nx,Ny,Maxiter, Nl, kt,i,q,p
  INTEGER      :: userChoice, nrec
  INTEGER      :: i1, iN, l1, lN, lMax, Me, Statinfo, Np, ideb, ifin
  INTEGER      :: myrank, ok, MsgTag
  REAL(KIND=8) :: Tmax,dt,dx,dy,CFL,t,Lx,Ly,Diff,erreur,erreurloc,time,tl
  REAL(KIND=8) :: coeff_a, coeff_bx, coeff_by
  REAL(KIND=8) :: epsilon
  REAL(KIND=8) :: tstart,tend
  CHARACTER(LEN=14) :: nameSol, nameSolBis
  INTEGER, DIMENSION (MPI_STATUS_SIZE)    :: Status
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: U, g, test
  REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: b, W, msg_d, msg_g

  Call MPI_INIT(Statinfo)
  MsgTag = 100
  tstart= MPI_WTIME()
  Call MPI_COMM_RANK(MPI_COMM_WORLD, Me, Statinfo)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, Np, Statinfo)

  !! Lecture des donnees
  unit = 10
  OPEN(unit,file='data1.txt', form='formatted', status='old')
  READ(unit,*) userChoice
  READ(unit,*) Nx,Ny
  READ(unit,*) Lx,Ly,Diff
  READ(unit,*) Tmax
  READ(unit,*) CFL
  READ(unit,*) epsilon,Nl,Ifreq
  READ(unit,*) nrec
  CLOSE(unit)

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
     IF (nrec>lMax .OR. nrec<0) THEN
        PRINT *,' '
        PRINT *,'----Taille de recouvrement incorrect----'
        PRINT *,' '
        STOP
     END IF
  END IF

  !! Allocation dynamique de la memoire
  If (Me == 0) THEN
     ALLOCATE(U(i1:iN),test(i1:iN))
     ALLOCATE(g(i1:iN),b(i1:iN),W(i1:iN))
     ALLOCATE(msg_d(Ny))
     ideb=i1
     ifin=iN
  Else
     ALLOCATE(U(i1-nrec*Ny:iN),test(i1-nrec*Ny:iN))
     ALLOCATE(g(i1-nrec*Ny:iN),b(i1-nrec*Ny:iN),W(i1-nrec*Ny:iN))
     ALLOCATE(msg_g(Ny))
     IF (Me < Np-1) THEN
        ALLOCATE(msg_d(Ny))
     END IF
     l1=l1-nrec
     ideb=i1-nrec*Ny
     ifin=iN
  END IF

  do i=ideb,ifin
     test(i) = i
  end do

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
  g  = 0.0d0
  W  = 0.0d0
  t  = 0.0d0

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
     dt = MIN (dt, Tmax-t)

     !! Calcul de W=AU
     CALL prodMatVec_W(W, U, coeff_a, coeff_bx, coeff_by, lN-l1+1, Ny, Me)
     !! Echange de message entre proc
     IF (Np>1) THEN
        IF (Me == 0) THEN
           CALL MPI_SEND(U(iN-(nrec+1)*Ny+1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Statinfo)
        ELSE IF (Me == Np-1) THEN
           CALL MPI_SEND(U(i1),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Statinfo)
        ELSE
           CALL MPI_SEND(U(iN-(nrec+1)*Ny+1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Statinfo)
           CALL MPI_SEND(U(i1),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Statinfo)
        END IF

        If (Me == 0) THEN
           CALL MPI_RECV(Msg_d(1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
        ELSE IF (Me == Np-1) THEN
           CALL MPI_RECV(Msg_g(1),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
        ELSE
           CALL MPI_RECV(Msg_d(1),Ny,MPI_REAL8,Me+1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
           CALL MPI_RECV(Msg_g(1),Ny,MPI_REAL8,Me-1,MsgTag,MPI_COMM_WORLD,Status,Statinfo)
        END IF
     END IF

     !! Calcul du second membre b
     DO i=ideb,ifin
        CALL indices(i, q, p, Ny)
        b(i) = dt*f(q*dx,p*dy,t,Lx,Ly,userChoice)+U(i)


        IF (q==l1) THEN
           IF (Me==0) THEN
              b(i) = b(i) + coeff_bx*funct_h((q-1)*dx,p*dy,t,userChoice)
           ELSE
              b(i) = b(i) + coeff_bx*msg_g(p)
           END IF
        END IF
        IF (q==lN) THEN
           IF (Me==Np-1) THEN
              b(i) = b(i) + coeff_bx*funct_h((q+1)*dx,p*dy,t,userChoice)
           ELSE
              b(i) = b(i) + coeff_bx*msg_d(p)
           END IF
        END IF
        IF (p==1) THEN
          b(i) = b(i) + coeff_by*funct_g(q*dx,(p-1)*dy,t,userChoice)
        END IF
        IF (p==Ny) THEN
           b(i) = b(i) + coeff_by*funct_g(q*dx,(p+1)*dy,t,userChoice)
        END IF

     END DO

     CALL gc(lN-l1+1,Ny,U,W,b,Nl,epsilon,coeff_a,coeff_bx,coeff_by, Me)

     t = t + dt;
  ENDDO

  !!Calcul de l'erreur seulement pour les choix 1 et 2
  IF (userChoice == 1 .OR. userChoice == 2) THEN
     erreurloc=0.
     erreur=0.
     DO i=i1,iN
        CALL indices(i, q, p, Ny)
        erreurloc = erreurloc +(U(i)-solExacte(q*dx,p*dy,t,userChoice))**2.
     END DO
     CALL MPI_ALLREDUCE(erreurloc,erreur,1,MPI_REAL8,&
         & MPI_SUM, MPI_COMM_WORLD, Statinfo)
     erreur = sqrt(erreur)
    !  IF (Me == 0) THEN
    !    PRINT *,"L'erreur est : ", erreur
    !  END IF
  END IF

   !! Ecriture des resultats
  unit = 11

  !!Appel de la subroutine rename pour renomer le fichier résultat selon le proc
  CALL Rename(Me,nameSol, nameSolBis)

  OPEN(unit,file=nameSol, form='formatted', status='unknown')
  DO i = i1,iN
     CALL indices(i, q, p, Ny)
     WRITE(unit,*) q*dx, p*dy, U(i)
  ENDDO
  CLOSE(unit)

  unit = 12
  OPEN(unit,file=nameSolBis, form='formatted', status='unknown')
  DO i = ideb,ifin
     CALL indices(i, q, p, Ny)
     WRITE(unit,*) q*dx, p*dy, U(i)
  ENDDO
  CLOSE(unit)

  DEALLOCATE(U,g,b,W)

  tend=MPI_WTIME()
  CALL MPI_ALLREDUCE(tend-tstart,time,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,Statinfo)
  IF (Me==0) THEN
    PRINT*, time
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

END PROGRAM Sequentiel1
