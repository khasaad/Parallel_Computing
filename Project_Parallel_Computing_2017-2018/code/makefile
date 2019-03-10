FC = mpif90

OBJS = gradient_conjugue.o main.o functions.o

main.exe : $(OBJS)
	  $(FC) $(OBJS) -o $@

%.o: %.f90
	$(FC) -fcheck=all -fbounds-check -O0 -ffpe-trap=invalid -pedantic -g -c -O $<

main.o : functions.o gradient_conjugue.o

clean:
	rm -f *~ *.o *.exe *.mod fort.10 sol*
