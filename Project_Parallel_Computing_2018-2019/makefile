FC = mpif90

OBJS1 = gradient_conjugue1.o main1.o functions1.o

OBJS2 = gradient_conjugue2.o main2.o functions2.o

main1.exe : $(OBJS1)
	  $(FC) $(OBJS1) -o $@

main2.exe : $(OBJS2)
		$(FC) $(OBJS2) -o $@

%.o: %.f90
	$(FC) -fcheck=all -fbounds-check -O0 -ffpe-trap=invalid -pedantic -g -c -O $<

main1.o : functions1.o gradient_conjugue1.o

main2.o : functions2.o gradient_conjugue2.o

clean:
	rm -f *~ *.o *.exe *.mod fort.10 sol*
