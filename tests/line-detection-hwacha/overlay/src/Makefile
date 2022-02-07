CC=~/firesim/target-design/chipyard/esp-tools-install/bin/riscv64-unknown-linux-gnu-gcc
#CC=riscv64-unknown-linux-gnu-gcc
CFLAGS=-O2 -static -march=rv64gcxhwacha -DUSE_HWACHA 
LDFLAGS=-static -lgcc -lm -lgcc -march=rv64gcxhwacha

OBJS =  main.o routinesCPU.o vec_add_asm.o vec_mul_asm.o


%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

all: line-detection-hwacha

line-detection-hwacha: $(OBJS)
	$(CC) -o line-detection-hwacha $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o line-detection-hwacha
