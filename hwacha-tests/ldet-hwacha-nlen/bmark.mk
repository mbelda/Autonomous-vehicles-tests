

#=======================================================================
# UCB CS250 Makefile fragment for benchmarks
#-----------------------------------------------------------------------
#
# Each benchmark directory should have its own fragment which
# essentially lists what the source files are and how to link them
# into an riscv and/or host executable. All variables should include
# the benchmark name as a prefix so that they are unique.
#

ldet-hwacha-nlen_c_src = \
        routines.c \
        main.c \
        syscalls.c \

ldet-hwacha-nlen_riscv_src = \
        crt.S \
        vec_mul_asm.S

ldet-hwacha-nlen_c_objs     = $(patsubst %.c, %.o, $(ldet-hwacha-nlen_c_src))
ldet-hwacha-nlen_riscv_objs = $(patsubst %.S, %.o, $(ldet-hwacha-nlen_riscv_src))

ldet-hwacha-nlen_host_bin = ldet-hwacha-nlen.host
$(ldet-hwacha-nlen_host_bin): $(ldet-hwacha-nlen_c_src)
        $(HOST_COMP) $^ -o $(ldet-hwacha-nlen_host_bin)

ldet-hwacha-nlen_riscv_bin = ldet-hwacha-nlen.riscv
$(ldet-hwacha-nlen_riscv_bin): $(ldet-hwacha-nlen_c_objs) $(ldet-hwacha-nlen_riscv_objs)
        $(RISCV_LINK) $(ldet-hwacha-nlen_c_objs) $(ldet-hwacha-nlen_riscv_objs) -o $(ldet-hwacha-nlen_riscv_bin) $(RISCV_LINK_OPTS)

junk += $(ldet-hwacha-nlen_c_objs) $(ldet-hwacha-nlen_riscv_objs) \
        $(ldet-hwacha-nlen_host_bin) $(ldet-hwacha-nlen_riscv_bin)
