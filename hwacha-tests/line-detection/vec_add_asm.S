#include "vec-util.h"

.text
.align 2

.globl vec_add_asm
.type  vec_add_asm,@function

# assumes calling convention:
# a0 has int n
# a1 has float* result  <---
# a2 has float* x
# a3 has float* y
vec_add_asm:
    li a4, VCFG(0, 0, 2, 1)
    vsetcfg a4
stripmine:
    vsetvl t0, a0 #a0 is requested vec len, actual is placed in t0
    vmca va0, a2
    vmca va1, a3
    la t5, vvadd_v
    vf 0(t5)
    slli t1, t0, 1
    add a2, a2, t1
    add a3, a3, t1
    sub a0, a0, t0
    bnez a0, stripmine
    fence
    ret

# vector thread asm
.align 3
vvadd_v:
    vpset vp0
    vlbu vv0, va0
    vlbu vv1, va1
    vadd vv1, vv0, vv1
    vsb vv1, va1
    vstop
