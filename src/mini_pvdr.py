#!/usr/bin/env python3
import sys
import heapq
import wave
import numpy as np


def wrap2Pi(x):
    return x - (2.0 * np.pi) * np.round(x / (2.0 * np.pi))

def calcPhase(pa0, cs0, cs1, Aa, As, M):
    assert len(pa0) == len(cs0) and len(cs0) == len(cs1)

    ratio = As / Aa
    empty = [True for _ in cs1]
    pa1 = np.full(cs1.shape, np.nan)
    heap = [(-np.abs(v), i, True) for i, v in enumerate(cs0)]
    heapq.heapify(heap)
    count = 0

    while count < len(cs1):
        _, i, isPrev = heapq.heappop(heap)
        if isPrev:
            if empty[i]:
                dp = wrap2Pi(np.angle(cs1[i]) - np.angle(cs0[i]) - (2.0 * np.pi * Aa / M) * i)
                pa1[i] = pa0[i] + ratio * dp + (2.0 * np.pi * As / M) * i
                empty[i] = False
                heapq.heappush(heap, (-np.abs(cs1[i]), i, False))
                count += 1
        else:
            if i >= 1 and empty[i - 1]:
                dp = wrap2Pi(np.angle(cs1[i - 1]) - np.angle(cs1[i]))
                pa1[i - 1] = pa1[i] + ratio * dp
                empty[i - 1] = False
                heapq.heappush(heap, (-np.abs(cs1[i - 1]), i - 1, False))
                count += 1
            if i < len(cs1) - 1 and empty[i + 1]:
                dp = wrap2Pi(np.angle(cs1[i + 1]) - np.angle(cs1[i]))
                pa1[i + 1] = pa1[i] + ratio * dp
                empty[i + 1] = False
                heapq.heappush(heap, (-np.abs(cs1[i + 1]), i + 1, False))
                count += 1

    assert not any(empty) and not np.any(np.isnan(pa1))
    return pa1

def stretchTime(fa, Aa, As, L, M):
    N = (len(fa) - L) // Aa
    window = np.sqrt(2.0 / 3.0) * np.sin((np.pi / L) * (1.0 / 2.0 + np.arange(L))) ** 2

    nfm = (1.0 / np.sqrt(M)) * np.exp((-2.0j * np.pi / M) * np.outer(np.arange(M // 2 + 1), np.arange(L) - (L - 1) / 2.0))
    ifm = (1.0 / np.sqrt(M)) * np.exp((+2.0j * np.pi / M) * np.outer(np.arange(L) - (L - 1) / 2.0, np.arange(M // 2 + 1)))
    ifm[:, 1:-1] *= 2.0
    #np.testing.assert_allclose((ifm @ nfm).real, np.eye(L), rtol = 1.0, atol = 1e-8)

    ca = np.empty([M // 2 + 1, N], np.complex128)
    for i in range(N):
        ca[:, i] = nfm @ (window * fa[Aa * i : Aa * i + L])

    cs = np.empty_like(ca)
    cs[:, 0] = ca[:, 0]
    ps = np.angle(ca[:, 0])
    for i in range(1, N):
        ps = calcPhase(ps, ca[:, i - 1], ca[:, i], Aa, As, M)
        cs[:, i] = np.abs(ca[:, i]) * np.exp(1.0j * ps)

    fs = np.zeros(As * N + L)
    for i in range(N):
        fs[As * i : As * i + L] += window * (ifm @ cs[:, i]).real

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

    fs = stretchTime(fa, int(np.round(1024.0 / float(sys.argv[1]))), 1024, 4096, 8192)

    fs = np.clip(np.round(32768.0 * fs), -32768.0, 32767.0).astype("<h")
    with wave.open(sys.argv[3], "wb") as f:
        f.setnchannels(1)
        f.setsampwidth(2)
        f.setframerate(rate)
        f.writeframes(fs)

main()
