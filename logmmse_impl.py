from scikits.audiolab import Sndfile, Format
import numpy as np
from scipy.special import *

np.seterr('raise')


def remove_absolute_silent(x):
    exist_zero_slice = False
    for i in range(len(x) - 60):
        if np.all(np.array(x[i:i+60])==0):
            x = np.delete(x,[i for i in range(i,i+60)])
            exist_zero_slice = True
            print "22"
            break
    return exist_zero_slice, x


def logmmse():
    file_name = "/home/sunya/Downloads/cs_jimi_wav/33670_22YS4VRIZ.wav"
    # file_name = "/home/sunya/Downloads/cs_jimi_wav/36418_CQDU35CLI.wav"
    out_name = "/home/sunya/out.wav"
    file_in = Sndfile(file_name, "r")

    sr = file_in.samplerate
    num_frames = file_in.nframes
    x = file_in.read_frames(num_frames)
    while True:
        exist_zero_slice, x = remove_absolute_silent(x)
        if not exist_zero_slice:
            break

    output_file = Sndfile(out_name, 'w',
                          Format(type=file_in.file_format, encoding='pcm16', endianness=file_in.endianness),
                          file_in.channels, sr)
    length = np.floor(20 * sr / 1000)  # frame size in samples
    if length % 2 == 1:
        length = length + 1
    PERC = 50  # window  overlap in percent of frame size
    len1 = np.floor(length * PERC / 100)
    len2 = length - len1
    w = np.hanning(length)  # defile window
    # Noise magnitude calculations - assuming that the first 6 frames is noise / silence
    nFFT = 2 * length
    noise_mean = np.zeros(int(nFFT), dtype=np.double)
    j = 1
    for m in range(0, 6):
        noise_mean = noise_mean + np.absolute(np.fft.fft(w * x[j:j + int(length)], int(nFFT), axis=0))
        j = j + int(length)
    noise_mu = noise_mean / 6
    noise_mu2 = np.power(noise_mu, 2)
    x_old = np.zeros((int(len1)), dtype=np.double)
    Xk_prev = np.zeros((int(len1), 1), dtype=np.double)
    Nframes = np.floor(len(x) / len2) - np.floor(length / len2)
    xfinal = np.zeros(int(Nframes) * int(len2), dtype=np.double)
    # start processing
    k = 1
    aa = 0.98
    mu = 0.98
    eta = 0.15

    ksi_min = 10 ^ (-25 / 10)

    for n in range(1, int(Nframes)):
        insign = w * x[k:k + int(length)]
        # print x[k:k+int(length)]
        spec = np.fft.fft(insign, int(nFFT))
        print spec
        sig = np.absolute(spec)  # compute the  magnitude
        sig2 = sig ** 2
        gammak = np.minimum(sig2 / noise_mu2, 40)  # limit post SNR to avoid overflows
        if n == 1:
            ksi = aa + (1 - aa) * np.maximum(gammak - 1, 0)
        else:
            ksi = aa * Xk_prev / noise_mu2 + (1 - aa) * np.maximum(gammak - 1, 0)  # a priori SNR
            ksi = np.maximum(ksi_min, ksi)  # limit ksi to - 25 dB

        log_sigma_k = gammak * ksi / (1 + ksi) - np.log(1 + ksi)
        vad_decision = np.sum(log_sigma_k) / length
        if vad_decision < eta:
            # noise only frame found
            noise_mu2 = mu * noise_mu2 + (1 - mu) * sig2
        # == =end of  vad == =
        A = ksi / (1 + ksi)  # % Log - MMSE  estimator
        vk = A * gammak
        ei_vk = 0.5 * expn(1, vk)
        hw = A * np.exp(ei_vk)
        sig = sig * hw
        Xk_prev = sig ** 2

        xi_w = np.fft.ifft(hw * spec, nFFT, axis=0)
        xi_w = np.real(xi_w)
        xfinal[k:k + int(len2)] = x_old + xi_w[0:int(len1)]
        x_old = xi_w[int(len1): int(length)]
        k = k + int(len2)
        # print xfinal
    output = np.array(xfinal * np.iinfo(np.int16).max, dtype=np.int16)
    output_file.write_frames(output)

logmmse()
# file_name = "/home/sunya/Downloads/cs_jimi_wav/33670_22YS4VRIZ.wav"
# file_in = Sndfile(file_name, "r")
# num_frames = file_in.nframes
# x = file_in.read_frames(num_frames)
# print len(x)
# exist_zero_slice,x=remove_absolute_silent(x)
# print len(x)
