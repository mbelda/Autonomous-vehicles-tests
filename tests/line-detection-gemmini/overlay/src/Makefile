CC=riscv64-unknown-linux-gnu-gcc
CFLAGS=-O2 -static
LDFLAGS=-static -lgcc -lm -lgcc

OBJS =  main.o


%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

all: line-detection-gemmini

line-detection-gemmini: $(OBJS)
	$(CC) -o line-detection-gemmini $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o line-detection-gemmini
