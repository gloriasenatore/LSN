CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio2.1.exe : esercizio2.1.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio2.1.o -o esercizio2.1.exe
esercizio2.1.o : esercizio2.1.cpp
	$(CC) -c esercizio2.1.cpp -o esercizio2.1.o $(CFLAGS)

esercizio2.2.exe : esercizio2.2.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio2.2.o -o esercizio2.2.exe
esercizio2.2.o : esercizio2.2.cpp
	$(CC) -c esercizio2.2.cpp -o esercizio2.2.o $(CFLAGS)

main.exe : main.o funzioni.o random.o
	$(CC) random.o funzioni.o main.o -o main.exe
main.o : main.cpp
	$(CC) -c main.cpp -o main.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)

clean :
	rm *.o main.exe seed.out
	rm *.o esercizio2.1.exe
	rm *.o esercizio2.2.exe
