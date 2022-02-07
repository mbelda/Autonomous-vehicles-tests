#=======================================================================
# UCB CS250 Makefile fragment for benchmarks
#-----------------------------------------------------------------------
#
# Each benchmark directory should have its own fragment which
# essentially lists what the source files are and how to link them
# into an riscv and/or host executable. All variables should include
# the benchmark name as a prefix so that they are unique.
#

line_detection_c_src = \
	main.c \
	routinesCPU.c \
	img0.c \
	syscalls.c \

line_detection_riscv_src = \
	crt.S \
	vec_add_asm.S \

line_detection_c_objs     = $(patsubst %.c, %.o, $(line_detection_c_src))
line_detection_riscv_objs = $(patsubst %.S, %.o, $(line_detection_riscv_src))

line_detection_host_bin = line-detection.host
$(line_detection_host_bin): $(line_detection_c_src)
	$(HOST_COMP) $^ -o $(line_detection_host_bin)

line_detection_riscv_bin = line-detection.riscv
$(line_detection_riscv_bin): $(line_detection_c_objs) $(line_detection_riscv_objs)
	$(RISCV_LINK) $(line_detection_c_objs) $(line_detection_riscv_objs) -o $(line_detection_riscv_bin) $(RISCV_LINK_OPTS)

junk += $(line_detection_c_objs) $(line_detection_riscv_objs) \
        $(line_detection_host_bin) $(line_detection_riscv_bin)
