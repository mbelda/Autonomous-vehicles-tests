CC=riscv64-unknown-linux-gnu-gcc
CFLAGS=-O2 -static
LDFLAGS=-static -lgcc -lm -lpthread -lgcc

OBJS = main.c \


%.o: %.c %.h
	$(CC) $(CFLAGS) -c $< -o $@

all: hello-world

hello-world: $(OBJS)
	$(CC) -o hello-world $(OBJS) $(LDFLAGS)

clean:
	rm -f *.o hello-world
