echo "phaml_graphics: "'\' 1>>$f
echo "	"'graphmain.o \' 1>>$f
echo "	pde.o" 1>>$f
echo "	"'$(LINKER) $(LINKFLAGS) -o phaml_graphics \' 1>>$f
echo "	"'graphmain.o \' 1>>$f
echo "	pde.o "'\' 1>>$f
echo "	"'-L$(PHAML_LIBDIR) -lphaml \' 1>>$f
if [ -n "$HYPRELIBS" ]
then
echo "	$HYPRELIBS "'\' 1>>$f
fi
if [ -n "$SUPERLULIBS" ]
then
echo "	$SUPERLULIBS "'\' 1>>$f
fi
if [ -n "$ARPACKLIBS" ]
then
echo "	$ARPACKLIBS "'\' 1>>$f
fi
if [ -n "$ZOLTANLIBS" ]
then
echo "	$ZOLTANLIBS "'\' 1>>$f
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
echo "	$MESSPASSLIBS "'\' 1>>$f
echo "	$OGLLIBS "'\' 1>>$f
if [ -n "$LAPACKLIBS" ]
then
echo "	$LAPACKLIBS "'\' 1>>$f
fi
if [ -n "$BLASLIBS" ]
then
echo "	$BLASLIBS "'\' 1>>$f
fi
if [ -n "$OTHERLIBS" ]
then
echo "	$OTHERLIBS "'\' 1>>$f
fi
echo "	$XLIBS" 1>>$f
echo "" 1>>$f

echo 'graphmain.o: $(PHAML_SRCDIR)/graphmain.f90' 1>>$f
echo "	"'$(F90) $(FFLAGS) '"$MODFLAG"'$(PHAML_MODDIR) -o graphmain.o -c $(PHAML_SRCDIR)/graphmain.f90' 1>>$f
echo "" 1>>$f
