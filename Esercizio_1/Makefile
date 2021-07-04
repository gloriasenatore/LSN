CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio1.exe : esercizio1.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio1.o -o esercizio1.exe
esercizio1.o : esercizio1.cpp
	$(CC) -c esercizio1.cpp -o esercizio1.o $(CFLAGS)

esercizio1.2.exe : esercizio1.2.o random.o funzioni.o
	$(CC) random.o esercizio1.2.o funzioni.o -o esercizio1.2.exe
esercizio1.2.o : esercizio1.2.cpp
	$(CC) -c esercizio1.2.cpp -o esercizio1.2.o $(CFLAGS)

esercizio1.3.exe : esercizio1.3.o random.o funzioni.o
	$(CC) random.o esercizio1.3.o funzioni.o -o esercizio1.3.exe
esercizio1.3.o : esercizio1.3.cpp
	$(CC) -c esercizio1.3.cpp -o esercizio1.3.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o esercizio1.exe
	rm *.o esercizio1.2.exe
	rm *.o esercizio1.3.exe

