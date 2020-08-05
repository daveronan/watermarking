/*
  ==============================================================================

    BufferTools.h
    Created: 16 Mar 2016 5:00:14pm
    Author:  cr

  ==============================================================================
*/

#ifndef BUFFERTOOLS_H_INCLUDED
#define BUFFERTOOLS_H_INCLUDED

#include "../JuceLibraryCode/JuceHeader.h"


class BufferTools
{
public:
	BufferTools();
	~BufferTools();
	void cleanAccumulatedBuffer();
	void prepareAccumulatedBuffer(int samplerate, int samplesperblock, int numChannels);

	void bufferAccumulator(AudioSampleBuffer buffer);

	void bufferAccumulatorOut(AudioSampleBuffer& buffer);

	
	bool isBuffer1Full;
	bool isBuffer2Full;


	bool isBufferEmpty;

	int windowSizeInBlocks;

	int blockCallCount;
	int blockCallCountOut;

	int buffer1Index;
	int buffer2Index;

	int bufferIndexOut;

	int bufferInUseIn;
	int bufferInUseOut;

	AudioSampleBuffer accumulatedBuffer1;
	AudioSampleBuffer accumulatedBuffer2;

private:
	
	int sampleRate;
	int samplesPerBlock;
	float windowSizeInMS;
	

};



#endif  // BUFFERTOOLS_H_INCLUDED
