CC = mpicxx
CFLAGS = -Wall -O3 --std=c++11

esercizio10.2.exe : esercizio10.2.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio10.2.o -o esercizio10.2.exe
esercizio10.2.o : esercizio10.2.cpp
	$(CC) -c esercizio10.2.cpp -o esercizio10.2.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o esercizio10.2.exe
