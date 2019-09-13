CC=CC
CFLAGS=-Wall --pedantic -std=c++11 -O3 -DLOG -g

.PHONY: all clean
all: testeventlogger

testeventlogger: testeventlogger.cpp eventlogger.hpp
	$(CC) $(CFLAGS) -o $@ $<

clean:
	rm testeventlogger *.log
