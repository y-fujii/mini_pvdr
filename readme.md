# mini_pvdr.py

mini_pvdr.py is a polyphonic audio time-stretcher and pitch-shifter using Phase Gradient Heap Integration method.  It is a small (&simeq; 100 SLOC) implementation, primarily intended for (my personal) learning/experiment/research purposes.

Note that several processes in this implementation are simplified from those described in the paper.  **Do not use it to evaluate the original paper.**

## Usage

```
pip install numpy
./src/mini_pvdr.py time_factor pitch_factor src.wav dst.wav
```

It supports 16-bit mono wave files only.

## Implementation notes

- It is important to align the origin of Fourier transform with the center of the window function, because the positions of wave packets are scaled in phase integration.

- mini_pvdr.py uses a simple nearest-neighbor difference scheme rather than the four-point difference scheme used in the paper.  While simple scheme yields better results in my tests, it may be caused by the incorrect implementation.

- mini_pvdr.py processes time-stretching and pitch-shifting simultaneously by using asymmetric FFT/inverse-FFT pairs, rather than resampling with time-stretching.  It requires non-power-of-two FFT.

## References

- [Z. Průša and N. Holighaus, "Phase Vocoder Done Right", 2017](https://doi.org/10.23919/EUSIPCO.2017.8081353) ([arXiv:2202.07382](https://arxiv.org/abs/2202.07382)).
