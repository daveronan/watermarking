/*
  ==============================================================================

    BufferTools.cpp
    Created: 16 Mar 2016 5:00:14pm
    Author:  cr

  ==============================================================================
*/

#include "BufferTools.h"
#include "Watermark.h"


BufferTools::BufferTools(): windowSizeInBlocks(0), blockCallCount(0), blockCallCountOut(0), sampleRate(0), samplesPerBlock(0), windowSizeInMS(0)
{
	buffer1Index = 0;	// used to set the start sample for copying from the default buffer to the larger one
	isBuffer1Full = false;
	buffer2Index = 0;	// used to set the start sample for copying from the default buffer to the larger one
	isBuffer2Full = false;

	isBufferEmpty = true;

	bufferIndexOut = 0;
	bufferInUseIn = 1;
	bufferInUseOut = 1;
}

BufferTools::~BufferTools()
{
}

void BufferTools::cleanAccumulatedBuffer()
{
	for(int i = 0; i < accumulatedBuffer1.getNumChannels(); i++)
	{
		float* bufferData = accumulatedBuffer1.getWritePointer(i);

		for(int j =0; j < windowSizeInBlocks; j++)
		{
			bufferData[j] = 0;
		}
	}

	for(int i = 0; i < accumulatedBuffer2.getNumChannels(); i++)
	{
		float* bufferData = accumulatedBuffer2.getWritePointer(i);

		for(int j =0; j < windowSizeInBlocks; j++)
		{
			bufferData[j] = 0;
		}
	}
}

void BufferTools::prepareAccumulatedBuffer(int samplerate, int samplesperblock, int numChannels)
{	
	sampleRate = samplerate;
	samplesPerBlock = samplesperblock;

	//Same as algorithm window size.
	windowSizeInBlocks = NFREQ / samplesPerBlock;	
	accumulatedBuffer1.setSize(numChannels, NFREQ, 0, 1, 1);
	accumulatedBuffer2.setSize(numChannels, NFREQ, 0, 1, 1);

	cleanAccumulatedBuffer();
}

/*
** buffer comes in from processblock
** make sure to set size of acccumulatedBuffer generously
** buffer is copied to accumulatedBuffer
** check to see if accumulatedBuffer has fulfilled windowsizeinMS
** if it has: set a flag to notify that it is full
** if not: continue filling
** create function for processing once accumulatedBuffer is full
*/

void BufferTools::bufferAccumulator(AudioSampleBuffer buffer)
{
	if(bufferInUseIn == 1)
	{
		if(buffer1Index + buffer.getNumSamples() >= accumulatedBuffer1.getNumSamples())
		{
			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				accumulatedBuffer1.copyFrom(i, buffer1Index, buffer, i, 0, accumulatedBuffer1.getNumSamples() - buffer1Index);
			}

			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				accumulatedBuffer2.copyFrom(i, 0, buffer, i, accumulatedBuffer1.getNumSamples() - buffer1Index,  buffer1Index + buffer.getNumSamples() - accumulatedBuffer1.getNumSamples());
			}

			bufferInUseIn = 2;
			buffer2Index = buffer1Index + buffer.getNumSamples() - accumulatedBuffer1.getNumSamples();
			isBuffer2Full = false;

			buffer1Index = 0;
			isBuffer1Full = true;
			
		}
		else
		{
			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				accumulatedBuffer1.copyFrom(i, buffer1Index, buffer, i, 0, buffer.getNumSamples());
			}

			buffer1Index += buffer.getNumSamples();
		}
			
	}
	else if(bufferInUseIn ==2)
	{
		if(buffer2Index + buffer.getNumSamples() >= accumulatedBuffer2.getNumSamples())
		{
			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				accumulatedBuffer2.copyFrom(i, buffer2Index, buffer, i, 0, accumulatedBuffer2.getNumSamples() - buffer2Index);
			}

			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				accumulatedBuffer1.copyFrom(i, 0, buffer, i, accumulatedBuffer2.getNumSamples() - buffer2Index,  buffer2Index + buffer.getNumSamples() - accumulatedBuffer2.getNumSamples());
			}

			bufferInUseIn = 1;
			buffer1Index = buffer2Index + buffer.getNumSamples() - accumulatedBuffer2.getNumSamples();
			isBuffer1Full = false;

			buffer2Index = 0;
			isBuffer2Full = true;
		}
		else
		{
			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				accumulatedBuffer2.copyFrom(i, buffer2Index, buffer, i, 0, buffer.getNumSamples());
			}

			buffer2Index += buffer.getNumSamples();
		}
	}
}

void BufferTools::bufferAccumulatorOut(AudioSampleBuffer& buffer)
{
	if(bufferInUseOut == 1)
	{
		if (bufferIndexOut + buffer.getNumSamples() >= accumulatedBuffer1.getNumSamples())
		{
			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				buffer.copyFrom(i, 0, accumulatedBuffer1, i, bufferIndexOut, accumulatedBuffer1.getNumSamples() - bufferIndexOut);
			}

			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				buffer.copyFrom(i, 0, accumulatedBuffer2, i, 0, bufferIndexOut + buffer.getNumSamples() - accumulatedBuffer1.getNumSamples());
			}

			bufferIndexOut = bufferIndexOut + buffer.getNumSamples() - accumulatedBuffer1.getNumSamples();
			bufferInUseOut = 2;
		}
		else
		{
			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				buffer.copyFrom(i, 0, accumulatedBuffer1, i, bufferIndexOut, buffer.getNumSamples());
			}

			bufferIndexOut += buffer.getNumSamples();

			//if (bufferIndexOut + buffer.getNumSamples() == accumulatedBuffer1.getNumSamples())
			//{
			//	bufferIndexOut = 0;
			//	bufferInUseOut = 2;
			//}
		}	
	}
	else if(bufferInUseOut == 2)
	{
		if (bufferIndexOut + buffer.getNumSamples() >= accumulatedBuffer2.getNumSamples())
		{
			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				buffer.copyFrom(i, 0, accumulatedBuffer2, i, bufferIndexOut, accumulatedBuffer2.getNumSamples() - bufferIndexOut);
			}

			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				buffer.copyFrom(i, 0, accumulatedBuffer1, i, 0, bufferIndexOut + buffer.getNumSamples() - accumulatedBuffer2.getNumSamples());
			}

			bufferIndexOut = bufferIndexOut + buffer.getNumSamples() - accumulatedBuffer2.getNumSamples();
			bufferInUseOut = 1;
		}
		else
		{
			for (int i = 0; i < buffer.getNumChannels(); i++)
			{
				buffer.copyFrom(i, 0, accumulatedBuffer2, i, bufferIndexOut, buffer.getNumSamples());
			}

			bufferIndexOut += buffer.getNumSamples();

			//if (bufferIndexOut + buffer.getNumSamples() == accumulatedBuffer2.getNumSamples())
			//{
			//	bufferIndexOut = 0;
			//	bufferInUseOut = 1;
			//}
		}		
	}
	
	
}
