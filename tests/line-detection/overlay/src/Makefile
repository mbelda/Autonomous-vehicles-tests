CC=riscv64-unknown-linux-gnu-gcc
CFLAGS=-O2 -static
LDFLAGS=-static -lgcc -lm -lgcc

OBJS =  main.o routinesCPU.o


%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

all: line-detection

line-detection: $(OBJS)
	$(CC) -o line-detection $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o line-detection
