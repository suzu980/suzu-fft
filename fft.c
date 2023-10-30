#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

int sample_rate = 64;
int signal_length = 64;
int N = 32;
float pi;

// FAST FOURIER TRANSFORM
void fft(float *x, int n, int s, complex float *out) {
  if (n <= 1) {
    out[0] = x[0];
		return;
  }
  // even
	float *even = x;
  fft(x, n / 2, s * 2, out);
  // odd
	float *odd = x;
  fft(x + s, n / 2, s * 2, out + n / 2);
  for (size_t k = 0; k < n / 2; ++k) {
    float a = -2 * (pi)*k / n;
    complex float t = out[k + n / 2] * (cos(a) + I * sin(a));
    complex float e = out[k];
    out[k] = e + t;
    out[(n / 2) + k] = e - t;
  }
}

// DISCRETE FOURIER TRANSFORM
void dft(float *x, int n, complex float *out, float *hann) {
  for (size_t f = 0; f < (int)(n / 2); ++f) {
    complex float new_dat = 0.0;
    for (size_t k = 0; k < n; ++k) {
      float y = x[k] *  hann[k];
      float a = -2 * (pi)*f * k / n;
      new_dat += y * (cos(a) + I * sin(a));
    }
    out[f] = new_dat;
  }
}

int main() {
  pi = atan(1) * 4;

  float hannCoeffs[N];
  for (size_t k = 0; k < N; ++k) {
    hannCoeffs[k] = 0.5 - 0.5 * cos(2 * pi * k / N - 1);
  }

  float xdata[signal_length];
  for (size_t i = 0; i < signal_length; ++i) {
    float t = (float)i / sample_rate;
    float f1 = 1;
    float f2 = 3;
    xdata[i] = sin(2 * pi * t * f1) + sin(2 * pi * t * f2);
  }
  // dft
  complex float dft_result[N / 2];
  dft(xdata, N, dft_result, hannCoeffs);
  for (size_t i = 0; i < (int)(N / 2); ++i) {
    float max_freq = (float)sample_rate / 2;
    float freq_per_bin = max_freq / ((float)N / 2.0);
    float q = i * freq_per_bin;
    float power_spec = cabsf(dft_result[i]);
    printf("DFT Frequency: %.2f\t%.2f\n", q, power_spec); 
  }
  printf("\n");
  // fft
  complex float fft_result[N / 2];
  float xdata2[N];
  for (size_t k = 0; k < N; ++k) {
    xdata2[k] = xdata[k] * hannCoeffs[k];
  }

  fft(xdata2, N, 1, fft_result);
  for (size_t i = 0; i < (int)(N / 2); ++i) {
    float max_freq = (float)sample_rate / 2;
    float freq_per_bin = max_freq / ((float)N / 2.0);
    float q = i * freq_per_bin;
    float power_spec = cabsf(fft_result[i]);
    printf("FFT Frequency: %.2f\t%.2f\n", q, power_spec); 
  }

  return 0;
}
