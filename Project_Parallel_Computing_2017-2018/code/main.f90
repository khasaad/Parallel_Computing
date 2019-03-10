PROGRAM Equation_de_la_chaleur

  use functions
  use gradient_conjugue

  IMPLICIT NONE
  Include 'mpif.h'
  
  !!Déclaration des variables
  INTEGER           :: unit,Ifreq
  INTEGER           :: Nx,Ny,Maxiter, Nl, kt, l,i,j,q,p
  INTEGER           :: myrank,ok
  INTEGER           :: i1, iN, Me, Statinfo, Np,ideb,ifin
  INTEGER           :: userChoice
  REAL(KIND=8)      :: Tmax,dt,dx,dy,CFL,t,Lx,Ly,Diff,erreur,erreurloc,time,time1
  REAL(KIND=8)      :: coeff_a, coeff_bx, coeff_by
  REAL(KIND=8)      :: epsilon,residu
  REAL(KIND=8)      :: tstart,tend
  CHARACTER(LEN=12) :: nameSol
  REAL(KIND=8), DIMENSION(:), allocatable :: U, g, test
  REAL(KIND=8), DIMENSION(:), allocatable :: b,W

  Call MPI_INIT(Statinfo) 
  tstart= MPI_WTIME()
  Call MPI_COMM_RANK(MPI_COMM_WORLD, Me, Statinfo)
  Call MPI_COMM_SIZE(MPI_COMM_WORLD, Np, Statinfo)

  !! Lecture des donnees
  unit = 10
  OPEN(unit,file='data.txt', form='formatted', status='old')
  READ(unit,*) userChoice
  READ(unit,*) Nx,Ny
  READ(unit,*) Lx,Ly,Diff
  READ(unit,*) Tmax
  READ(unit,*) CFL
  READ(unit,*) epsilon,Nl,Ifreq
  CLOSE(unit)

  !Vérification de la taille de la matrice par rapport au nombre de procs
  IF (Np>Nx .OR. Np>Ny) THEN
     PRINT *,'----Taille du maillage incompatible avec le nombre de procs----'
     STOP
  END IF

  
  !! Verification du choix du probleme
  IF (userChoice <1 .OR. userChoice>3) THEN
     PRINT *,'----Choix du probleme invalide----'
     PRINT *,'Il n''y a que trois problemes resolvable par ce programme, numerotes de 1 à 3'
     PRINT *,'----------------------------------'
     STOP
  END IF

  CALL charge(Me, Np, Nx, i1, iN)
  
  !! Allocation dynamique de la memoire
  ALLOCATE(g(i1:iN),b(i1:iN),W(i1-1:iN+1))
  If (Me == 0) THEN
     ALLOCATE(U(i1:iN+Ny),test(i1:iN+Ny))
     ideb=i1
     ifin=iN+Ny
  ELSE IF (Me == Np-1) THEN
     ALLOCATE(U(i1-Ny:iN),test(i1-Ny:iN))
     ideb=i1-Ny
     ifin=iN
  Else
     ALLOCATE(U(i1-Ny:iN+Ny),test(i1-Ny:iN+Ny))
     ideb=i1-Ny
     ifin=iN+Ny
  END IF
   
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
  test = 0.d0
  g  = 0.0d0
  W  = 0.0d0 
  t  = 0.0d0

  do i=i1,iN
     test(i) = i
  end do
  
  !! Mise en mémoire des trois coefficients de la
  !! matrice A du système linéiare A*U=B à résoudre
  coeff_a = 1.+2.*Diff*dt*(1./(dx*dx)+1./(dy*dy))
  coeff_bx = Diff*dt/(dx*dx)
  coeff_by = Diff*dt/(dy*dy)

  !! Détermination du nombre d'itérations à effectuer
  !! pour calculer la solution au temps final Tmax
  Maxiter= Int(Tmax/dt)
  IF( (Tmax - Maxiter*dt) > 1.d-6 ) Maxiter=Maxiter+1 
  
  !! Boucle en temps
  DO kt = 1, Maxiter
     dt = MIN (dt, Tmax-t)
     
     
     !! Calcul de W=AU
     CALL prodMatVec_W(W, U, coeff_a, coeff_bx, coeff_by, Nx, Ny, i1, iN, ideb, ifin, Statinfo, Me, Np)     
     
     !! Calcul du second membre b
     DO i=i1,iN
        CALL indices(i, q, p, Ny)
        b(i) = dt*f(q*dx,p*dy,t,Lx,Ly,userChoice)+U(i)
                
        IF (q==1) THEN
           b(i) = b(i) + coeff_bx*funct_h((q-1)*dx,p*dy,t,userChoice)
        END IF
        IF (q==Nx) THEN
          b(i) = b(i) + coeff_bx*funct_h((q+1)*dx,p*dy,t,userChoice)
        END IF
        IF (p==1) THEN
          b(i) = b(i) + coeff_by*funct_g(q*dx,(p-1)*dy,t,userChoice)
        END IF
        IF (p==Ny) THEN
           b(i) = b(i) + coeff_by*funct_g(q*dx,(p+1)*dy,t,userChoice)
        END IF
       
     END DO
     
     !print*,b(i1:iN),"-----------------"

     CALL gc(Nx,Ny,U,W,b,Nl,epsilon,coeff_a,coeff_bx,coeff_by,i1,iN,ideb,ifin,Statinfo,Me,Np)
     
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
     PRINT *,"L'erreur est : ", erreur
  END IF
  

  !! Ecriture des resultats
  unit = 11

  !!Appel de la subroutine rename pour renomer le fichier résultat selon le proc
  call Rename(Me,nameSol)

  OPEN(unit,file=nameSol, form='formatted', status='unknown')
  DO i = i1,iN
     CALL indices(i, q, p, Ny)
     WRITE(unit,*) q*dx, p*dy, U(i)
  ENDDO
  CLOSE(unit)  

  tend=MPI_WTIME()
  CALL MPI_ALLREDUCE(tend-tstart,time,1,MPI_REAL8,MPI_MAX,MPI_COMM_WORLD,Statinfo)

  print*,time
  DEALLOCATE(U,g,b,W)


  !!Ecriture des différents temps de calcul
  !!Pour le calcul du speedup et de l'efficacité on a besoin de connaitre le temps
  !!d'execution sur 1 processeur, pour cela si on lance le programme sur 1 proc
  !!on sauvegarde ce temps dans le fichier "proc"
  if(Me==0) then
     OPEN(unit,file='proc', form='formatted', status='unknown')
     if(Np==1)then
        WRITE(unit,*) time
        time1=time
     end if
     CLOSE(unit)

     !!Si on à lancé le programme sur plus de 1 processeur, on va cherche le temps d'exectution sur 1 proc dans le fichier "proc"
     !!et on ecrit ensuite dans le fichier temps le nombre de proc le speedup et l'efficacité
     !!Il peut y avoir une erreur à l'execution si on a supprimé le fichier proc et qu'on lance le programme
     !!directement sur plus de 1 processeur mais les fichiers résultas sont quand même crées
     OPEN(unit,file='proc', form='formatted', status='unknown',iostat = ok)
     if(Np/=0) then
        READ(unit,*,iostat=ok)time1
     end if
     CLOSE(unit)

     OPEN(unit,file='temps', form='formatted', status='old', position='append')
     WRITE(unit,*) Np, time1/(Np*time),time1/time
     CLOSE(unit)  
  end if

  CALL MPI_FINALIZE(Statinfo)


CONTAINS
  Subroutine Charge(Me,Np,N, i1, iN)
    INTEGER, Intent (IN)  :: Me, Np, N
    INTEGER, Intent (OUT) :: i1, iN
    INTEGER               :: r, q
    
    q = N/Np
    r = N-q*Np
    
    IF (Me <r) THEN
       i1=Me*(q+1)+1
       iN=(Me+1)*(q+1)
    ELSE
       i1=1+r+Me*q
       iN=i1+q-1
    END IF
    
    i1 = (i1-1)*Ny+1
    iN = iN*Ny
    
  END Subroutine Charge

  subroutine Rename(Me,name)
    integer           :: Me
    character(LEN=12) :: name
    character(LEN=3)  :: tn
    integer :: i1,i2,i3
    i1 = Me/100
    i2 =( Me - 100*i1)/10
    i3 = Me - 100*i1 -10*i2
    tn = char(i1+48)//char(i2+48)//char(i3+48)
    name='sol'//tn//'.dat'
  end subroutine Rename

END PROGRAM Equation_de_la_chaleur
