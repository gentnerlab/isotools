# $Id: Makefile,v 1.3 2011/01/06 22:21:47 samn Exp $ 


# from http://www.metalshell.com/view/tutorial/120/
# Automated Makefile



CC = g++

CFLAGS = -Wall -O2

COMPILE = $(CC) $(CFLAGS) -c

OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))



all: isoi



isoi: $(OBJFILES)

	$(CC) -o isoi $(OBJFILES)



%.o: %.cpp

	$(COMPILE) -o $@ $<

clean:
	rm -rf *.o

install:
	mv isoi ../../