CPP=g++ 
INCLUDE=-I./include -I/usr/include/mysql
LIBS=
CPPFLAGS=-g -O3 -w $(INCLUDE) $(LIBS)
PROGRAM=SFMULT

all: $(PROGRAM)

OBJ = \
    sfmult_anxnyn.o \
    sfmult_anxnyt.o \
    sfmult_anxnyt_k.o \
    sfmult_anxtyn.o \
    sfmult_anxtyt.o \
    sfmult_anxtyt_k.o \
    sfmult_atxnyn.o \
    sfmult_atxnyt.o \
    sfmult_atxtyn.o \
    sfmult_atxtyn_k.o \
    sfmult_atxtyt.o \
    sfmult_atxtyt_k.o \
    sfmult.o \
    sfmult_vector_1.o \
    sfmult_vector_k.o \
    sfmult_xA.o \
    ssmult_transpose.o

$(OBJ): sfmult.h

.c.o:
	g++ -O -c $<

$(PROGRAM): $(OBJ)

clean:
	rm -f *.o
