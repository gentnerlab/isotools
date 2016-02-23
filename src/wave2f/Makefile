# $Id: Makefile,v 1.2 2011/01/15 20:24:52 samn Exp $ 


# from http://www.metalshell.com/view/tutorial/120/
# Automated Makefile



CC = g++

CFLAGS = -Wall -O2

COMPILE = $(CC) $(CFLAGS) -c

OBJFILES := $(patsubst %.cpp,%.o,$(wildcard *.cpp))



all: wave2f



wave2f: $(OBJFILES)

	$(CC) -o wave2f $(OBJFILES)



%.o: %.cpp

	$(COMPILE) -o $@ $<

clean:
	rm -rf *.o

install:
	mv wave2f ../../