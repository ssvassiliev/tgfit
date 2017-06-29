# Makefile for TgFit on LINUX # 
CC=gcc
# Intel: CFLAGS = -mkl -xHost
CFLAGS = -O3 -Wno-unused-result --fast-math -I/usr/local/include/plplot -I/opt/intel/mkl/include/ 
LIB_DIR =  -L/opt/intel/mkl/lib/intel64 -L/usr/local/lib/plplot
LIBS = -lgsl -lblas -llapack -lplplot -lmkl_rt -lm
BINDIR = /usr/local/bin
#CONV=conv.c
CONV=conv-mkl.o

OBJ1 = tgfit-lg.o utils.o kinetic.o $(CONV)  derivs.o mrqmin.o mem-alloc.o matr.o spline.o gra.o load.o save.o
OBJ2 = tgfit-auto.o utils.o kinetic.o $(CONV) derivs.o mrqmin.o mem-alloc.o matr.o spline.o load.o save.o
OBJ3 = mod_ans.o utils.o kinetic.o $(CONV)  derivs.o mrqmin.o mem-alloc.o matr.o spline.o load.o save.o
OBJ4 = tgfit-random-init.c utils.o kinetic.o $(CONV)  derivs.o mrqmin.o mem-alloc.o matr.o spline.o load.o save.o
OBJ5 = tgfit-grow-fit.c utils.o kinetic.o $(CONV) derivs.o mrqmin.o mem-alloc.o matr.o spline.o load.o save.o

tgfit : $(OBJ1)
	$(CC)  $(CFLAGS) $(OBJ1) -o bin/tgfit $(LIB_DIR) $(LIBS) 
auto : $(OBJ2)
	$(CC) $(CFLAGS) $(OBJ2) -o bin/tgfit-auto $(LIB_DIR) $(LIBS)
mod : $(OBJ3)
	$(CC) $(CFLAGS) $(OBJ3) -o bin/modans $(LIB_DIR) $(LIBS)
rand : $(OBJ4)
	$(CC) $(CFLAGS) $(OBJ4) -o bin/rand $(LIB_DIR) $(LIBS)
grow : $(OBJ5)
	$(CC) $(CFLAGS) $(OBJ5) -o bin/grow $(LIB_DIR) $(LIBS)
all : tgfit auto mod rand grow
install : 
	cd bin; \ 
	cp * $(BINDIR); cp -r .tgfit ~
tar :
	tar -cvf tgfit.tar *.c *.h Makefile .tgfitrc data ; \
	gzip tgfit.tar
clean :
	rm -f bin/* *.o; cd data; rm -f DAS.dat residuals.dat autocorr.dat \
	model.dat derivs.dat fit.dat trace.dat
