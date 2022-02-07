CC=riscv64-unknown-linux-gnu-gcc
CFLAGS=-O2 -static
LDFLAGS=-static -lgcc -lm -lgcc

OBJS = img0.o main.o routinesCPU.o


%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

all: line-detection-bare-br

line-detection-bare-br: $(OBJS)
	$(CC) -o line-detection-bare-br $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o line-detection-bare-br
