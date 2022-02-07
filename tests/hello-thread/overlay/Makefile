CC=riscv64-unknown-linux-gnu-gcc
CFLAGS=-O2 -static -march=rv64gc
LDFLAGS=-static -lgcc -lm -lpthread -lgcc -march=rv64gc

OBJS = main.c \


%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

all: hello-thread

hello-thread: $(OBJS)
	$(CC) -o hello-thread $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o hello-thread
