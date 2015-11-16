CC=g++   
CFLAGS=-O2 -g  
IFLAGS=-I /usr/include/opencv 
LDFLAGS=-lm -lcv -lhighgui -lcvaux -lcxcore
THREAD=-pthread


main: imageOperations.o multiThread.o activeContour.o main.o
	${CC} ${CFLAGS} -o main main.o activeContour.o imageOperations.o multiThread.o ${IFLAGS} ${LDFLAGS} ${THREAD}

activeContour.o : activeContour.cpp activeContour.h imageOperations.h multiThread.h
	${CC}  -c ${CFLAGS} activeContour.cpp imageOperations.o multiThread.o ${IFLAGS} ${LDFLAGS} ${THREAD}

imageOperations.o : imageOperations.cpp imageOperations.h
	${CC}  -c ${CFLAGS} imageOperations.cpp ${IFLAGS} ${LDFLAGS} ${THREAD}

multiThread.o : multiThread.cpp multiThread.h imageOperations.h	
	${CC}  -c ${CFLAGS} multiThread.cpp imageOperations.o ${IFLAGS} ${LDFLAGS} ${THREAD}

main.o : main.cpp
	${CC}  -c ${CFLAGS} main.cpp ${IFLAGS} ${LDFLAGS} ${THREAD}

clean :
	rm -rf *.o 
	rm main

all : 
	make main
