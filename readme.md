# mini_pvdr.py

mini_pvdr.py is a minimal (< 100 SLOC) implementation of "Phase Vocoder Done Right", primarily for (my personal) research/learning purposes.

Note that several processes are simplified in this implementation.  The result must be worse than the original.  Do not use it to evaluate the method.

## Usage

```
./src/mini_pvdr.py stretch_factor src.wav dst.wav
```

It supports 16-bit mono wave files only.

## References

- [Z. Průša and N. Holighaus, "Phase Vocoder Done Right", 2017](https://doi.org/10.23919/EUSIPCO.2017.8081353).
