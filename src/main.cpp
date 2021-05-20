
#include "bitstream_autocorrelation.hpp"
#include "daisy_pod.h"
#define _USE_MATH_DEFINES
#include <cmath>

using namespace daisy;

static DaisyPod pod;
static BitstreamCorrelation correlation;

float testFrequency = 440;
unsigned int sampleCount = 0;

void AudioCallback(float* in, float* out, size_t size)
{
	auto buff_size = size / 2;

	float outl, outr, inl, inr;
	for (size_t i = 0; i < size; i += 2)
	{
		float value = sin(sampleCount++ / 44100.0f * testFrequency * 2 * pi());
		inl = value; // in[i];
		inr = value; // in[i + 1];

		// left out
		out[i] = inl;

		// right out
		out[i + 1] = inr;

		correlation.addSample(inl);
		float guess = correlation.guessFrequency();
		if (guess != 0.0f)
		{
			pod.seed.PrintLine("guessing %fhz\n", guess);
		}
	}
}

int main(void)
{
	// initialize pod hardware and oscillator daisysp module
	float sample_rate;

	// Inits and sample rate
	pod.Init();
	pod.seed.StartLog(true);

	// start callback
	pod.StartAdc();
	pod.StartAudio(AudioCallback);

	pod.seed.PrintLine("starting audio");

	while (1)
	{
	}
}
