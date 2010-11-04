# Makefile for 3D FDTD Program

CC = /usr/bin/gcc
OBJS = fdtd.o normal.o mur1st.o 
CFLAGS = -Wall -O2 #-fopenmp
#CFLAGS = -g -Wall -O0 
DEBUGFLAG = -pg
INCLUDES = -I. -I/opt/local/include
LIBS = -L. -L/opt/local/lib 
TARG = fdtd
LIBOPT = -lm

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
