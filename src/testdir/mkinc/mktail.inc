echo "pde.o: pde.f90" 1>>$f
echo "	"'$(F90) $(FFLAGS) '"$MODFLAG"'$(PHAML_MODDIR)'" $ZOLTANMOD -c pde.f90" 1>>$f
echo "" 1>>$f

echo "alittleclean:" 1>>$f
echo "	rm -f test*.exe phaml_slave" 1>>$f

echo "clean:" 1>>$f
echo "	rm -f *.o *.mod test*.exe phaml_slave phaml_graphics test*.out test*.err test*.diff" 1>>$f
