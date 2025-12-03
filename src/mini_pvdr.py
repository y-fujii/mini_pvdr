#!/usr/bin/env python3
import sys
import cmath
import heapq
import wave
import numpy as np


def calcPhase(pa0, cs0, cs1, ra, rs):
    assert len(pa0) == len(cs0) and len(cs0) == len(cs1)

    ratio = rs / ra
    empty = [True for _ in cs1]
    pa1 = np.full(cs1.shape, np.nan)
    heap = [(-abs(v), i, True) for i, v in enumerate(cs0)]
    heapq.heapify(heap)
    count = 0

    while count < len(cs1):
        _, i, isPrev = heapq.heappop(heap)
        if isPrev:
            if empty[i]:
                dp = cmath.phase((cs1[i] / cs0[i]) * cmath.exp((-2.0j * cmath.pi * ra) * i))
                pa1[i] = pa0[i] + ratio * dp + (2.0 * cmath.pi * rs) * i
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

def process(fa, Aa, La, As, Ls, Mf):
    assert La % 2 == 0 and Ls % 2 == 0
    N = (len(fa) - La) // Aa

    window = np.sqrt(2.0 / 3.0) * np.sin((np.pi / La) * np.arange(La)) ** 2
    ca = np.empty([La * Mf // 2 + 1, N], np.complex128)
    for i in range(N):
        windowed = window * fa[Aa * i : Aa * i + La]
        ca[:, i] = np.fft.rfft(np.r_[windowed[La // 2 :], np.zeros(La * Mf - La), windowed[: La // 2]], norm = "forward")

    cs = np.empty_like(ca)
    cs[:, 0] = ca[:, 0]
    ps = np.angle(ca[:, 0])
    for i in range(1, N):
        ps = calcPhase(ps, ca[:, i - 1], ca[:, i], Aa / (La * Mf), As / (Ls * Mf))
        cs[:, i] = np.abs(ca[:, i]) * np.exp(1.0j * ps)

    window = np.sqrt(2.0 / 3.0) * np.sin((np.pi / Ls) * np.arange(Ls)) ** 2
    fs = np.zeros(As * N + Ls)
    for i in range(N):
        inversed = np.fft.irfft(np.r_[cs[: Ls * Mf // 2 + 1, i], np.zeros(max((Ls - La) * Mf // 2, 0))], norm = "forward")
        fs[As * i : As * i + Ls] += window * np.r_[inversed[-Ls // 2 :], inversed[: Ls // 2]]

    return fs

def main():
    if len(sys.argv) != 5:
        print(f"Usage: {sys.argv[0]} time_factor pitch_factor src.wav dst.wav")
        return

    with wave.open(sys.argv[3]) as f:
        assert f.getnchannels() == 1 and f.getsampwidth() == 2
        rate = f.getframerate()
        fa = np.frombuffer(f.readframes(f.getnframes()), "<h")
    fa = fa.astype(np.float64) / 32768.0

    kt, kp = float(sys.argv[1]), float(sys.argv[2])
    fs = process(fa, round(1024.0 / (kt * kp)), 4096, round(1024.0 / kp), 4 * round(1024.0 / kp), 2)

    fs = np.clip(np.round(32768.0 * fs), -32768.0, 32767.0).astype("<h")
    with wave.open(sys.argv[4], "wb") as f:
        f.setnchannels(1)
        f.setsampwidth(2)
        f.setframerate(rate)
        f.writeframes(fs)

main()
