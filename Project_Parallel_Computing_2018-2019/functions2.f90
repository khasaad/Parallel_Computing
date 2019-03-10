Module functions2

  IMPLICIT NONE

CONTAINS
  !! Fonction du second membre
  FUNCTION f(x, y, t, Lx, Ly, num_pb)
    INTEGER, Intent(IN)      :: num_pb
    REAL(kind=8), Intent(IN) :: x, y, t, Lx, Ly
    REAL(kind=8)             :: f

    IF (num_pb == 1) THEN
       f = 2.*(y-y*y+x-x*x)
    ELSE IF (num_pb == 2) THEN
       f = sin(x)+cos(y)
    ELSE
       f = exp(-(x-0.5*Lx)**2.)*exp(-(y-0.5*Ly)**2.)*cos(3.141592653589*0.5*t)
    END IF

  END FUNCTION f


  !! Fonctions determinant les conditions aux bords
  FUNCTION funct_g(x, y, t, num_pb)
    INTEGER, Intent(IN)      :: num_pb
    REAL(kind=8), Intent(IN) :: x, y, t
    REAL(kind=8)             :: funct_g

    IF (num_pb == 2) THEN
       funct_g = sin(x) + cos(y)
    ELSE
       funct_g = 0.
    END IF

  END FUNCTION funct_g

  FUNCTION funct_h(x, y, t, num_pb)
    INTEGER, Intent(IN)      :: num_pb
    REAL(kind=8), Intent(IN) :: x, y, t
    REAL(kind=8)             :: funct_h

    IF (num_pb == 1) THEN
       funct_h = 0.
    ELSE IF (num_pb == 2) THEN
       funct_h = sin(x)+cos(y)
    ELSE
       funct_h = 1
    END IF

  END FUNCTION funct_h


  !! Solutions excactes des deux premiers problèmes
  FUNCTION solExacte(x, y, t, num_pb)
    INTEGER, Intent(IN)      :: num_pb
    REAL(kind=8), Intent(IN) :: x, y, t
    REAL(kind=8)             :: solExacte

    IF (num_pb == 1) THEN
       solExacte = x*(1.-x)*y*(1.-y)
    ELSE
       solExacte = sin(x)+cos(y)
    END IF

  END FUNCTION solExacte

END MODULE functions2
