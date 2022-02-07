## Hwacha v4 SAXPY ASM code

#include "vec-util.h"

.text
.align 2

.globl vec_mul_asm
.type  vec_mul_asm,@function

# assumes calling convention:
# a0 has int n
# a1 has int* x
# a2 has int* y
vec_saxpy_asm:
    li t0, VCFG(0, 2, 0, 1) #Solo necesito 2 vectores de enteros y 1 predicado (64, 32, 16, pred)
    vsetcfg t0
    vsetvl t0, a0 #a0 is requested vec len, actual is placed in t0
    vmca va0, a1
    vmca va1, a2
    la t5, mul_v
    vf 0(t5)
    fence
    ret

# vector thread asm
.align 3
mul_v:
    #vpset vp0
    vlw vv0, va0
    vlw vv1, va1
    # appears to be faster if we write into vv1 instead of using a new reg?
    #vfmadd.s vv1, vs1, vv0, vv1
    vmul vv1, vv0, vv1 #Multiplico a*b
    vsw vv1, va1
    vstop
    