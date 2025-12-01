#!/usr/bin/env python3
import sys
import cmath
import heapq
import wave
import numpy as np


def calcPhase(pa0, cs0, cs1, Aa, As, M2):
    assert len(pa0) == len(cs0) and len(cs0) == len(cs1)

    ratio = As / Aa
    empty = [True for _ in cs1]
    pa1 = np.full(cs1.shape, np.nan)
    heap = [(-abs(v), i, True) for i, v in enumerate(cs0)]
    heapq.heapify(heap)
    count = 0

    while count < len(cs1):
        _, i, isPrev = heapq.heappop(heap)
        if isPrev:
            if empty[i]:
                dp = cmath.phase((cs1[i] / cs0[i]) * cmath.exp((-1.0j * cmath.pi * Aa / M2) * i))
                pa1[i] = pa0[i] + ratio * dp + (cmath.pi * As / M2) * i
                empty[i] = False
                heapq.heappush(heap, (-abs(cs1[i]), i, False))
                count += 1
        else:
            if i >= 1 and empty[i - 1]:
                dp = cmath.phase(cs1[i - 1] / cs1[i])
                pa1[i - 1] = pa1[i] + ratio * dp
                empty[i - 1] = False
                heapq.heappush(heap, (-abs(cs1[i - 1]), i - 1, False))
                count += 1
            if i < len(cs1) - 1 and empty[i + 1]:
                dp = cmath.phase(cs1[i + 1] / cs1[i])
                pa1[i + 1] = pa1[i] + ratio * dp
                empty[i + 1] = False
                heapq.heappush(heap, (-abs(cs1[i + 1]), i + 1, False))
                count += 1

    assert not any(empty) and not np.any(np.isnan(pa1))
    return pa1

def stretchTime(fa, Aa, As, L2, M2):
    N = (len(fa) - 2 * L2) // Aa
    window = np.sqrt(2.0 / 3.0) * np.sin((np.pi / (2 * L2)) * np.arange(2 * L2)) ** 2

    ca = np.empty([M2 + 1, N], np.complex128)
    for i in range(N):
        windowed = window * fa[Aa * i : Aa * i + 2 * L2]
        ca[:, i] = np.fft.rfft(np.concatenate([windowed[L2 :], np.zeros(2 * M2 - 2 * L2), windowed[: L2]]))

    cs = np.empty_like(ca)
    cs[:, 0] = ca[:, 0]
    ps = np.angle(ca[:, 0])
    for i in range(1, N):
        ps = calcPhase(ps, ca[:, i - 1], ca[:, i], Aa, As, M2)
        cs[:, i] = np.abs(ca[:, i]) * np.exp(1.0j * ps)

    fs = np.zeros(As * N + 2 * L2)
    for i in range(N):
        inversed = np.fft.irfft(cs[:, i])
        fs[As * i : As * i + 2 * L2] += window * np.concatenate([inversed[2 * M2 - L2 :], inversed[: L2]]).real

    return fs

def main():
    if len(sys.argv) != 4:
        print(f"Usage: {sys.argv[0]} stretch_factor src.wav dst.wav")
        return

    with wave.open(sys.argv[2]) as f:
        assert f.getnchannels() == 1 and f.getsampwidth() == 2
        rate = f.getframerate()
        fa = np.frombuffer(f.readframes(f.getnframes()), "<h")
    fa = fa.astype(np.float64) / 32768.0

    fs = stretchTime(fa, int(np.round(1024.0 / float(sys.argv[1]))), 1024, 2048, 4096)

    fs = np.clip(np.round(32768.0 * fs), -32768.0, 32767.0).astype("<h")
    with wave.open(sys.argv[3], "wb") as f:
        f.setnchannels(1)
        f.setsampwidth(2)
        f.setframerate(rate)
        f.writeframes(fs)

main()
