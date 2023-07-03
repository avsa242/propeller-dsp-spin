{
    --------------------------------------------
    Filename: dsp.fft.spin
    Author: Jesse Burt
    Description: Fast Fourier Transform
    Started Jan 25, 2011
    Updated Jul 03, 2023
    See end of file for terms of use.
    --------------------------------------------

    NOTE: This is based on heater_fft.spin,
    originally by Michael Rychlik, modified for formatting
    and unused code removed for simplicity and clarity

    Specifications:
        In-place radix-2 decimation in time FFT
        1024-point FFT: 34ms (* on original demo test data)

    Usage:
        1) Application should reserve two 1024 arrays of LONGs for input and output data, x and y
        2) Place input signal into the x array (values 13-bit signed; -4096..4095)
        3) Clear the y array.
        4) Call butterflies() giving a command and addresses of the arrays as parameters.
        5) The command is one or more of the following bits set:
               CMD_DECIMATE   - Perform bit-reversal reordering on the data, results in x and y
               CMD_BUTTERFLY  - Perform the actual FFT butterfly calculations, results in x and y
               CMD_MAGNITUDE  - Convert resulting x and y values to frequency magnitudes in x

        6) For some applications it may be as well to write the input data directly into x
            LONG by LONG in the correct bit-reversed order. Then remove the CMD_DECIMATE in the
            butterfly() call.
            This could move 10 percent or so of the processing time to the input COG.

    References:
        https://forums.parallax.com/discussion/128292/heaters-fast-fourier-transform/p1
        https://forums.parallax.com/discussion/127306/fourier-for-dummies-under-construction/p1
}

#define USE_FASTER_MULT
#define USE_FASTER_SQRT

CON
    'Specify size of FFT buffer here with length and log base 2 of the length.
    'N.B. Changing this will require changing the "twiddle factor" tables.
    '     and may also require changing the fixed point format (if going bigger)
    FFT_SIZE        = 1024
    FFT_RANGE       = 4096                      ' amplitude
    LOG2_FFT_SIZE   = 10

    CMD_DECIMATE    = %0001
    CMD_BUTTERFLY   = %0010
    CMD_MAGNITUDE   = %0100
    CMD_TEST        = %1000

VAR

    long mailboxp

PUB start(mailp)

    mailboxp := mailp
    long[mailboxp] := 0
    cognew(@bfly, mailp)     'Check error?

PUB butterflies(cmd, bxp, byp)

    long[mailboxp + 4] := bxp                       ' Address of x buffer
    long[mailboxp + 8] := byp                       ' Address of y buffer
    long[mailboxp{+0}] := cmd                       ' Do butterflies and/or decimation
    repeat while long[mailboxp]

DAT
                org     0
bfly            mov     mb_ptr, par
                rdlong  command, mb_ptr wz          ' Wait for run command in mailbox
        if_z    jmp     #bfly

                add     mb_ptr, #4                  ' Fetch x array address from mbox
                rdlong  bx_ptr, mb_ptr

                add     mb_ptr, #4                  ' Fetch y array address from mbox
                rdlong  by_ptr, mb_ptr
                sub     mb_ptr, #8

                test    command, #CMD_DECIMATE wz   ' Bit reversal required on data?
        if_z    jmp     #:no_decimate

' Radix-2 decimation in time. (The bit reversal stage)
' Moves every sample of bx to a postion given by reversing the bits of its original array index.
' This is a direct translation of the Spin decimate above, original Spin code used as comments.
' N.B. Only the x array is bit-reversed it is up to the app to clear y.

                mov     c, fft_size_                ' repeat i from 0 to FFT_SIZE - 1
                mov     b, #0

:dloop          mov     a, b                        ' revi := i >< LOG2_FFT_SIZE
                mov     rev_a, a
                rev     rev_a, #32 - LOG2_FFT_SIZE

                cmp     a, rev_a wc                 ' if i < revi
        if_nc   jmp     #:skip_rev

                shl     a, #2                       ' Times 4 as we are reading longs
                shl     rev_a, #2

                mov     hub_ptr, bx_ptr             ' tx1 := long[bxp + i * 4]
                add     hub_ptr, a
                rdlong  tx, hub_ptr

                mov     hub_rev_ptr, bx_ptr         ' long[bxp + i * 4] := long[bxp + revi * 4]
                add     hub_rev_ptr, rev_a
                rdlong  ty, hub_rev_ptr
                wrlong  ty, hub_ptr

                wrlong  tx, hub_rev_ptr             ' long[bxp + revi * 4] := tx1

:skip_rev       add     b, #1
                djnz    c, #:dloop

:no_decimate
                test    command, #CMD_BUTTERFLY wz  ' Perform butterflies?
        if_z    jmp     #:no_butterfly

' Apply FFT butterflies to N complex samples in buffers bx and by, in time decimated order!
' Resulting FFT is produced in bx and by in the correct order.
' This is a direct translation from the Spin code above, original Spin code in comments.

                mov     flight_max, fft_size_       ' flight_max := FFT_SIZE / 2
                sar     flight_max, #1
                mov     wangleSkip, fft_size_       ' wangleSkip := FFT_SIZE * 4
                shl     wangleSkip, #2

                mov     butterflySpan, #4           ' butterflySpan := 4
                mov     butterfly_max, #1           ' butterfly_max := 1
                mov     flightSkip, #4              ' flightSkip := 4

                ' Loop through all the decimation levels
                mov     level, #LOG2_FFT_SIZE       ' level := LOG2_FFT_SIZE
:lloop                                              ' repeat
                mov     b0x_ptr, bx_ptr             ' b0x_ptr := @bx
                mov     b0y_ptr, by_ptr             ' b0y_ptr := @by

                mov     b1x_ptr, b0x_ptr            ' b1x_ptr := b0x_ptr + butterflySpan
                add     b1x_ptr, butterflySpan

                mov     b1y_ptr, b0y_ptr            ' b1y_ptr := b0y_ptr + butterflySpan
                add     b1y_ptr, butterflySpan

                ' Loop though all the flights in a level
                mov     flight, flight_max          ' flight := flight_max
:floop                                              ' repeat
{new}           mov     wangle, #0

                ' Loop through all the butterflies in a flight
                mov     butterfly, butterfly_max    ' butterfly := butterfly_max

                ' Do the initial pass optimization, when W = [1,0] we don't need to multiply.
                ' c = 1 (well, 4096/4096), d = 0
                mov     k2, #0                      ' k2 := (d * (a + b)) / 4096
                rdlong  a, b1x_ptr                  ' a := long[b1x_ptr]
                mov     k1, a                       ' k1 := (a * (c + d)) / 4096
                neg     k3, a                       ' k3 := (c * (b - a)) / 4096
                rdlong  b, b1y_ptr                  ' b := long[b1y_ptr]
                add     k3, b                       ' k3 := (c * (b - a)) / 4096 (cont.)
                jmp     #:continue_bloop
              
:bloop          ' repeat                            ' At last...the butterfly.
                rdlong  a, b1x_ptr                  ' a := long[b1x_ptr]

                ' Precompute the optimization for c=0, d=-1
                neg     k1, a                       ' k1 := (a * (c + d)) / 4096
                neg     k2, a                       ' k2 := (d * (a + b)) / 4096

                rdlong  b, b1y_ptr                  ' b := long[b1y_ptr]

                ' Precompute the optimization for c=0, d=-1
                sub     k2, b                       ' k2 := (d * (a + b)) / 4096 (cont.)
                mov     k3, #0                      ' k3 := (c * (b - a)) / 4096

                mov     c, wangle
{getcos}        add     c, sin_90                   ' For cosine, add 90°
                test    c, sin_90      wc           ' Get quadrant 2|4 into c
                test    c, sin_180     wz           ' Get quadrant 3|4 into nz
                negc    c, c                        ' If quadrant 2|4, negate offset
                or      c, sin_table                ' OR in sin table address >> 1
                shl     c, #1                       ' Shift left to get final word address
                rdword  c, c                        ' Read word sample from $E000 to $F000
                negnz   c, c                        ' If quadrant 3|4, negate sample

                sar     c, #4 wz                    ' Scale to +/- 4095

        if_z    jmp     #:continue_bloop            ' if c==0, we already kave k1, k2, k3 calculated

                mov     d, wangle
{getsin}        test    d, sin_90      wc           ' Get quadrant 2|4 into c
                test    d, sin_180     wz           ' Get quadrant 3|4 into nz
                negc    d, d                        ' If quadrant 2|4, negate offset
                or      d, sin_table                ' OR in sin table address >> 1
                shl     d, #1                       ' Shift left to get final word address
                rdword  d, d                        ' Read word sample from $E000 to $F000
                negnz   d, d                        ' If quadrant 3|4, negate sample

                sar     d, #4                       ' Scale to +/- 4095
                neg     d, d                        ' We want -cos

                mov     m1, c                       ' k1 := (a * (c + d)) / 4096
                add     m1, d
                mov     m2, a
                call    #mul
                mov     k1, m1
                sar     k1, #15 - 3

                mov     m1, a                       ' k2 := (d * (a + b)) / 4096
                add     m1, b
                mov     m2, d
                call    #mul
                mov     k2, m1
                sar     k2, #15 - 3

                mov     m1, b                       ' k3 := (c * (b - a)) / 4096
                sub     m1, a
                mov     m2, c
                call    #mul
                mov     k3, m1
                sar     k3, #15 - 3

:continue_bloop

                mov     tx, k1                      ' tx := k1 - k2 (part I)
                mov     ty, k1                      ' ty := k1 + k3 (part I)

                rdlong  k1, b0x_ptr                 ' k1 := long[b0x_ptr]

                sub     tx, k2                      ' (part II) moved from above to take advantage of the hub wait times
                add     ty, k3                      ' ditto

                rdlong  k2, b0y_ptr                 ' k2 := long[b0y_ptr]

                mov     a, k1                       ' long[b1x_ptr] := k1 - tx
                sub     a, tx
                wrlong  a, b1x_ptr

                mov     a, k2                       ' long[b1y_ptr] := k2 - ty
                sub     a, ty
                wrlong  a, b1y_ptr

                mov     a, k1                       ' long[b0x_ptr] := k1 + tx
                add     a, tx
                wrlong  a, b0x_ptr

                mov     a, k2                       ' long[b0y_ptr] := k2 + ty
                add     a, ty
                wrlong  a, b0y_ptr

                add     b0x_ptr, #4                 ' b0x_ptr += 4
                add     b0y_ptr, #4                 ' b0y_ptr += 4

                add     b1x_ptr, #4                 ' b1x_ptr += 4
                add     b1y_ptr, #4                 ' b1y_ptr += 4

                add     wangle, wangleSkip          ' wangle += wangleSkip

                djnz    butterfly, #:bloop          ' while --butterfly <> 0

                add     b0x_ptr, flightSkip         ' b0x_ptr += flightSkip
                add     b0y_ptr, flightSkip         ' b0y_ptr += flightSkip
                add     b1x_ptr, flightSkip         ' b1x_ptr += flightSkip
                add     b1y_ptr, flightSkip         ' b1y_ptr += flightSkip
                djnz    flight, #:floop             ' while --flight <> 0

                shl     butterflySpan, #1           ' butterflySpan <<= 1
                shl     flightSkip, #1              ' flightSkip <<= 1

                shr     flight_max, #1              ' flight_max >>= 1

                shr     wangleSkip, #1
                shr     wSkip, #1                   ' wSkip >>= 1
                shl     butterfly_max, #1           ' butterfly_max <<= 1
                djnz    level, #:lloop              ' while --level <> 0
:no_butterfly
                test    command, #CMD_MAGNITUDE wz  ' Calculate magnitudes?
        if_z    jmp     #:no_magnitude

' Calculate magnitudes from the complex results in x and y. Results placed into x

                mov     c, fft_size_                ' repeat i from 0 to FFT_SIZE
                add     c, #1                       ' That is one more than half FFT_SIZE so as
                                                    ' to include the Nyquist frequency
                mov     b0x_ptr, bx_ptr
                mov     b0y_ptr, by_ptr

:mloop          rdlong  m1, b0x_ptr
                sar     m1, #LOG2_FFT_SIZE - 1
                mov     m2, m1
                call    #mul
                mov     input, m1

                rdlong  m1, b0y_ptr
                sar     m1, #LOG2_FFT_SIZE - 1
                mov     m2, m1
                call    #mul
                add     input, m1

                call    #sqrt

                wrlong  root, b0x_ptr               ' Write result to x array

                add     b0x_ptr, #4                 ' Next x and y element and loop
                add     b0y_ptr, #4
                djnz    c, #:mloop

:no_magnitude
                mov     command, #0
                wrlong  command, mb_ptr
                jmp     #bfly

mul             ' Account for sign
#ifdef USE_FASTER_MULT
                abs     m1, m1 wc
                negc    m2, m2
                abs     m2, m2 wc
                ' Make t2 the smaller of the 2 unsigned parameters
                mov     m3, m1
                max     m3, m2
                min     m2, m1
                ' Correct the sign of the adder
                negc    m2, m2
#else
                abs     m3, m1 wc
                negc    m2, m2
#endif
                ' My accumulator
                mov     m1, #0
                ' Do the work
:mul_loop       shr     m3, #1 wc,wz                ' Get the low bit of t2
        if_c    add     m1, m2                      ' If it was a 1, add adder to accumulator
                shl     m2, #1                      ' Shift the adder left by 1 bit
        if_nz   jmp     #:mul_loop                  ' Continue as long as there are no more 1's
mul_ret         ret

m1              long    0
m2              long    0
m3              long    0

#ifdef USE_FASTER_SQRT
' Faster code for square root (Chip Gracey after discussion with lonesock on Propeller Forums):
sqrt            mov     root, h40000000
                cmpsub  input, root  wc
                sumnc   root, h40000000
                shr     root, #1

                or      root, h10000000
                cmpsub  input, root  wc
                sumnc   root, h10000000
                shr     root, #1

                or      root, h04000000
                cmpsub  input, root  wc
                sumnc   root, h04000000
                shr     root, #1

                or      root, h01000000
                cmpsub  input, root  wc
                sumnc   root, h01000000
                shr     root, #1

                or      root, h00400000
                cmpsub  input, root  wc
                sumnc   root, h00400000
                shr     root, #1

                or      root, h00100000
                cmpsub  input, root  wc
                sumnc   root, h00100000
                shr     root, #1

                or      root, h00040000
                cmpsub  input, root  wc
                sumnc   root, h00040000
                shr     root, #1

                or      root, h00010000
                cmpsub  input, root  wc
                sumnc   root, h00010000
                shr     root, #1

                or      root, h00004000
                cmpsub  input, root  wc
                sumnc   root, h00004000
                shr     root, #1

                or      root, h00001000
                cmpsub  input, root  wc
                sumnc   root, h00001000
                shr     root, #1

                or      root, h00000400
                cmpsub  input, root  wc
                sumnc   root, h00000400
                shr     root, #1

                or      root, #$100
                cmpsub  input, root  wc
                sumnc   root, #$100
                shr     root, #1

                or      root, #$40
                cmpsub  input, root  wc
                sumnc   root, #$40
                shr     root, #1

                or      root, #$10
                cmpsub  input, root  wc
                sumnc   root, #$10
                shr     root, #1

                or      root, #$4
                cmpsub  input, root  wc
                sumnc   root, #$4
                shr     root, #1

                or      root, #$1
                cmpsub  input, root  wc
                sumnc   root, #$1
                shr     root, #1
sqrt_ret        ret

h10000000       long    $10000000
h04000000       long    $04000000
h01000000       long    $01000000
h00400000       long    $00400000
h00100000       long    $00100000
h00040000       long    $00040000
h00010000       long    $00010000
h00004000       long    $00004000
h00001000       long    $00001000
h00000400       long    $00000400

#else

sqrt            mov     root, #0                    ' Reset root
                mov     mask, h40000000             ' Reset mask (constant in register)
:sqloop         or      root, mask                  ' Set trial bit
                cmpsub  input, root wc              ' Subtract root from input if fits
                sumnc   root, mask                  ' Cancel trial bit, set root bit if fit
                shr     root, #1                    ' Shift root down
                shr     mask, #2                    ' Shift mask down
                tjnz    mask, #:sqloop              ' Loop until mask empty
sqrt_ret        ret
#endif

h40000000       long    $40000000

' Large constants
fft_size_       long    FFT_SIZE
sin_90          long    $0800
sin_180         long    $1000
sin_table       long    $E000 >> 1                  ' ROM sin table base shifted right

{ COG variables }
level           long    0
flight          long    0
butterfly       long    0
flight_max      long    0
wSkip           long    0
butterflySpan   long    0
butterfly_max   long    0
flightSkip      long    0
k1              long    0
k2              long    0
k3              long    0
a               long    0
b               long    0
c               long    0
d               long    0
tx              long    0
ty              long    0
b0x_ptr         long    0
b0y_ptr         long    0
b1x_ptr         long    0
b1y_ptr         long    0
mb_ptr          long    0
bx_ptr          long    0
by_ptr          long    0
wangle          long    0
wangleSkip      long    0

rev_a           long    0
hub_ptr         long    0
hub_rev_ptr     long    0
command         long    0
root            long    0
mask            long    0
input           long    0

                fit     496

DAT
{
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
}

