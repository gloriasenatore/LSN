CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio8.1.exe : esercizio8.1.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio8.1.o -o esercizio8.1.exe
esercizio8.1.o : esercizio8.1.cpp
	$(CC) -c esercizio8.1.cpp -o esercizio8.1.o $(CFLAGS)

esercizio8.2.exe : esercizio8.2.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio8.2.o -o esercizio8.2.exe
esercizio8.2.o : esercizio8.2.cpp
	$(CC) -c esercizio8.2.cpp -o esercizio8.2.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o esercizio8.1.exe
	rm *.o esercizio8.2.exe
