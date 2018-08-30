CC = gcc
CFLAGS = -O3 -Wall
#CFLAGS = -g -Wall
LFLAGS = -lm
PROGRAMS = OpenGPTv1
OBJECTS  = gptMain.o multiMatch.o singleMatch.o gpt.o utility.o eval.o

.SUFFIXES: .c .o

all:    $(PROGRAMS)

$(PROGRAMS): $(OBJECTS)
	$(CC) -o $(PROGRAMS) $(OBJECTS) $(LFLAGS) 

.c.o:
	$(CC) $(CFLAGS) -c $<

.PHONY: clean
clean: 
	rm -f *.o $(PROGRAMS)

eval.o: parameters.h eval.h gpt.h utility.h
gpt.o:  parameters.h eval.h gpt.h utility.h
gptMain.o:      parameters.h eval.h gpt.h multiMatch.h singleMatch.h utility.h
multiMatch.o:   parameters.h eval.h gpt.h multiMatch.h utility.h
singleMatch.o:  parameters.h eval.h gpt.h singleMatch.h multiMatch.h utility.h
utility.o:      parameters.h eval.h gpt.h utility.h

