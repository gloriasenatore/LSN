CC = g++
CFLAGS = -Wall -O3 --std=c++11

esercizio5_unif.exe : esercizio5_unif.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio5_unif.o -o esercizio5_unif.exe
esercizio5_unif.o : esercizio5_unif.cpp
	$(CC) -c esercizio5_unif.cpp -o esercizio5_unif.o $(CFLAGS)

esercizio5_gauss.exe : esercizio5_gauss.o random.o funzioni.o
	$(CC) random.o funzioni.o esercizio5_gauss.o -o esercizio5_gauss.exe
esercizio5_gauss.o : esercizio5_gauss.cpp
	$(CC) -c esercizio5_gauss.cpp -o esercizio5_gauss.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
funzioni.o : funzioni.cpp funzioni.h
	$(CC) -c funzioni.cpp -o funzioni.o $(CFLAGS)
clean :
	rm *.o main.exe seed.out
	rm *.o esercizio5_unif.exe
	rm *.o esercizio5_gauss.exe
