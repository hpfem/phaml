if [ ! $PHAML_PARALLEL = "sequential" -a ! $PHAML_PARALLEL = "messpass_nospawn" ]
then
echo "phaml_slave: "'\' 1>>$f
echo "	"'slave.o \' 1>>$f
echo "	pde.o" 1>>$f
echo "	"'$(LINKER) $(LINKFLAGS) -o phaml_slave \' 1>>$f
echo "	"'slave.o \' 1>>$f
echo "	pde.o "'\' 1>>$f
echo "	"'-L$(PHAML_LIBDIR) -lphaml \' 1>>$f
if [ -n "$HYPRELIBS" ]
then
echo "	$HYPRELIBS "'\' 1>>$f
fi
if [ -n "$ZOLTANLIBS" ]
then
echo "	$ZOLTANLIBS "'\' 1>>$f
fi
if [ -n "$ARPACKLIBS" ]
then
echo "	$ARPACKLIBS "'\' 1>>$f
fi
if [ -n "$PETSCLIBS" ]
then
echo "	$PETSCLIBS "'\' 1>>$f
fi
if [ -n "$BLOPEXLIBS" ]
then
echo "	$BLOPEXLIBS "'\' 1>>$f
fi
if [ -n "$MUMPSLIBS" ]
then
echo "	$MUMPSLIBS "'\' 1>>$f
fi
if [ -n "$SUPERLULIBS" ]
then
echo "	$SUPERLULIBS "'\' 1>>$f
fi
if [ -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$OTHERLIBS" ]
then
   echo "	$MESSPASSLIBS "'\' 1>>$f
else
   echo "	$MESSPASSLIBS" 1>>$f
fi
if [ -n "$LAPACKLIBS" ]
then
   if [ -n "$BLASLIBS" -o -n "$OTHERLIBS" ]
   then
      echo "	$LAPACKLIBS "'\' 1>>$f
   else
      echo "	$LAPACKLIBS " 1>>$f
   fi
fi
if [ -n "$BLASLIBS" ]
then
   if [ -n "$OTHERLIBS" ]
   then
      echo "	$BLASLIBS "'\' 1>>$f
   else
      echo "	$BLASLIBS " 1>>$f
   fi
fi
if [ -n "$OTHERLIBS" ]
then
echo "	$OTHERLIBS" 1>>$f
fi
echo "" 1>>$f
fi

if [ $PHAML_PARALLEL = "messpass_spawn" ]
then
echo 'slave.o: $(PHAML_SRCDIR)/slave.f90' 1>>$f
echo "	"'$(F90) $(FFLAGS) '"$MODFLAG"'$(PHAML_MODDIR)'" $ZOLTANMOD "'-o slave.o -c $(PHAML_SRCDIR)/slave.f90' 1>>$f
echo "" 1>>$f

else

echo 'phaml_slave:' 1>>$f

fi
