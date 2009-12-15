echo "$TESTN.exe: "'\' 1>>$f
echo "	$TESTN.o pde.o" 1>>$f
echo "	"'$(LINKER) $(LINKFLAGS)'" -o $TESTN.exe "'\' 1>>$f
echo "	$TESTN.o pde.o "'\' 1>>$f
if [ -n "$SUPERLULIBS" -o -n "$HYPRELIBS" -o -n "$ZOLTANLIBS" -o -n "$ARPACKLIBS" -o -n "$PETSCLIBS" -o -n "BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
then
echo "	"'-L$(PHAML_LIBDIR) -lphaml \' 1>>$f
else
echo "	"'-L$(PHAML_LIBDIR) -lphaml' 1>>$f
fi
if [ -n "$SUPERLULIBS" ]
then
   if [ -n "$HYPRELIBS" -o -n "$ZOLTANLIBS" -o -n "$ARPACKLIBS" -o -n "$PETSCLIBS" -o -n "BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "  $SUPERLULIBS "'\' 1>>$f
   else
echo "  $SUPERLULIBS " 1>>$f
   fi
fi
if [ -n "$HYPRELIBS" ]
then
   if [ -n "$ZOLTANLIBS" -o -n "$ARPACKLIBS" -o -n "$PETSCLIBS" -o -n "BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$HYPRELIBS "'\' 1>>$f
   else
echo "	$HYPRELIBS " 1>>$f
   fi
fi
if [ -n "$ZOLTANLIBS" ]
then
   if [ -n "$ARPACKLIBS" -o -n "$PETSCLIBS" -o -n "BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$ZOLTANLIBS "'\' 1>>$f
   else
echo "	$ZOLTANLIBS " 1>>$f
   fi
fi
if [ -n "$ARPACKLIBS" ]
then
   if [ -n "$PETSCLIBS" -o -n "BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$ARPACKLIBS "'\' 1>>$f
   else
echo "	$ARPACKLIBS " 1>>$f
   fi
fi
if [ -n "$PETSCLIBS" ]
then
   if [ -n "BLOPEXLIBS" -o -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$PETSCLIBS "'\' 1>>$f
   else
echo "	$PETSCLIBS " 1>>$f
   fi
fi
if [ -n "$BLOPEXLIBS" ]
then
   if [ -n "$MUMPSLIBS" -o -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$BLOPEXLIBS "'\' 1>>$f
   else
echo "	$BLOPEXLIBS " 1>>$f
   fi
fi
if [ -n "$MUMPSLIBS" ]
then
   if [ -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$MESSPASSLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$MUMPSLIBS "'\' 1>>$f
   else
echo "	$MUMPSLIBS " 1>>$f
   fi
fi
if [ -n "$MESSPASSLIBS" ]
then
   if [ -n "$LAPACKLIBS" -o -n "$BLASLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$MESSPASSLIBS "'\' 1>>$f
   else
echo "	$MESSPASSLIBS " 1>>$f
   fi
fi
# OGLLIBS must come before anything that gives a directory containing glut or gl
if [ -n "$PHAML_GETS_GRAPHICSLIBS" ]
then
echo "  $OGLLIBS "'\' 1>>$f
fi
if [ -n "$LAPACKLIBS" ]
then
   if [ -n "$BLASLIBS" -o -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$LAPACKLIBS "'\' 1>>$f
   else
echo "	$LAPACKLIBS " 1>>$f
   fi
fi
if [ -n "$BLASLIBS" ]
then
   if [ -n "$PHAML_GETS_GRAPHICSLIBS" -o -n "$OTHERLIBS" ]
   then
echo "	$BLASLIBS "'\' 1>>$f
   else
echo "	$BLASLIBS " 1>>$f
   fi
fi
if [ -n "$OTHERLIBS" ]
then
   if [ -n "$PHAML_GETS_GRAPHICSLIBS" ]
   then
echo "	$OTHERLIBS "'\' 1>>$f
   else
echo "	$OTHERLIBS " 1>>$f
   fi
fi
if [ -n "$PHAML_GETS_GRAPHICSLIBS" ]
then
echo "	$XLIBS" 1>>$f
fi
echo "" 1>>$f

# write the rules for compiling files

echo "$TESTN.o: $TESTN.f90 pde.o" 1>>$f
echo "	"'$(F90) $(FFLAGS) '"$MODFLAG"'$(PHAML_MODDIR)'" $ZOLTANMOD -c $TESTN.f90" 1>>$f
echo "" 1>>$f
