# mini_pvdr.py

mini_pvdr.py is a minimal (< 100 SLOC) implementation of "Phase Vocoder Done Right", primarily for (my personal) research/learning purposes.

Note that several processes are simplified in this implementation.  DO NOT USE IT TO EVALUATE THE ORIGINAL PAPER.

## Usage

```
./src/mini_pvdr.py stretch_factor src.wav dst.wav
```

It supports 16-bit mono wave files only.

## Implementation notes

- It is important to align the origin of fourier transform with the center of the window, because the positions of wave packets are scaled in phase integration.

- mini_pvdr.py uses a simple nearest-neighbor difference scheme rather than the four-point difference scheme used in the paper. While simple scheme yields better results in my tests, it may be caused by the incorrect implementation.

## References

- [Z. Průša and N. Holighaus, "Phase Vocoder Done Right", 2017](https://doi.org/10.23919/EUSIPCO.2017.8081353).
