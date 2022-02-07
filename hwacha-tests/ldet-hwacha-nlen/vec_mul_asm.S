## Hwacha v4 VVADD ASM code

#include "vec-util.h"

.text
.align 2

.globl vec_mul_asm
.type  vec_mul_asm,@function

# assumes calling convention:
# a0 has int n
# a1 has int* x
# a2 has int* y
vec_mul_asm:
    li a4, VCFG(0, 2, 0, 1) #(64, 32, 16, pred)
    vsetcfg a4
stripmine:
    vsetvl t0, a0 #a0 is requested vec len, actual is placed in t0
    vmca va0, a1
    vmca va1, a2
    la t5, vmul_v
    vf 0(t5)
    slli t1, t0, 1
    add a1, a1, t1
    add a2, a2, t1
    sub a0, a0, t0
    bnez a0, stripmine
    fence
    ret

# vector thread asm
.align 3
vmul_v:
    vpset vp0
    vlw vv0, va0
    vlw vv1, va1
    vmul vv1, vv0, vv1
    vsw vv1, va1
    vstop
