# Makefile for 3D FDTD Program

CC = /usr/bin/gcc
OBJS = control.o fdtd.o normal.o mur1st.o
WFLAGS = -Wall -W -Wformat=2 -Wcast-qual -Wcast-align -Wwrite-strings -Wconversion -Wfloat-equal -Wpointer-arith
CFLAGS = $(WFLAGS) -O2 -fopenmp -I/svtool/local/include
#CFLAGS = -g -Wall -O0
DEBUGFLAG = -pg
INCLUDES = -I. -I/opt/local/include
LIBS = -L. -L/opt/local/lib -L/svtool/local/lib
TARG = control
LIBOPT = -lm -lglut -lGLU

all: $(TARG)
$(TARG): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) $(LIBOPT) -o $@ $(OBJS)

.c.o:
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -o $@ -c $<

.PHONY: clean
clean:
	rm $(OBJS) $(TARG)

test:
	./$(TARG) 1.0
