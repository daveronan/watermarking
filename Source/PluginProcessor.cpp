/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#include "PluginProcessor.h"
#include "PluginEditor.h"



//==============================================================================
AudioWatermarkingAudioProcessor::AudioWatermarkingAudioProcessor()
{	
	_watermarkToApply = "";
}

AudioWatermarkingAudioProcessor::~AudioWatermarkingAudioProcessor()
{
}

//==============================================================================
const String AudioWatermarkingAudioProcessor::getName() const
{
    return JucePlugin_Name;
}

bool AudioWatermarkingAudioProcessor::acceptsMidi() const
{
   #if JucePlugin_WantsMidiInput
    return true;
   #else
    return false;
   #endif
}

bool AudioWatermarkingAudioProcessor::producesMidi() const
{
   #if JucePlugin_ProducesMidiOutput
    return true;
   #else
    return false;
   #endif
}

double AudioWatermarkingAudioProcessor::getTailLengthSeconds() const
{
    return 0.0;
}

int AudioWatermarkingAudioProcessor::getNumPrograms()
{
    return 1;   // NB: some hosts don't cope very well if you tell them there are 0 programs,
                // so this should be at least 1, even if you're not really implementing programs.
}

int AudioWatermarkingAudioProcessor::getCurrentProgram()
{
    return 0;
}

void AudioWatermarkingAudioProcessor::setCurrentProgram (int index)
{
}

const String AudioWatermarkingAudioProcessor::getProgramName (int index)
{
    return String();
}

void AudioWatermarkingAudioProcessor::changeProgramName (int index, const String& newName)
{
}

//==============================================================================
void AudioWatermarkingAudioProcessor::prepareToPlay (double sampleRate, int samplesPerBlock)
{
    // Use this method as the place to do any pre-playback
    // initialisation that you need..
	
	_bufferTools.prepareAccumulatedBuffer(sampleRate, samplesPerBlock, this->getNumInputChannels());
	_watermark.initialize(this->getSampleRate(), this->getNumInputChannels(), _bufferTools.accumulatedBuffer1.getNumSamples());
	_watermark.readWMBits(_watermarkToApply);
	this->setNonRealtime(true);
	_samplesPerBlock = samplesPerBlock;
	this->suspendProcessing(false);
}

void AudioWatermarkingAudioProcessor::releaseResources()
{
    // When playback stops, you can use this as an opportunity to free up any
    // spare memory, etc.
	this->suspendProcessing(true);
	_watermark.cleanUp();
	_bufferTools.cleanAccumulatedBuffer();
	
}

#ifndef JucePlugin_PreferredChannelConfigurations
bool AudioWatermarkingAudioProcessor::setPreferredBusArrangement (bool isInput, int bus, const AudioChannelSet& preferredSet)
{
    // Reject any bus arrangements that are not compatible with your plugin

    const int numChannels = preferredSet.size();

   #if JucePlugin_IsMidiEffect
    if (numChannels != 0)
        return false;
   #elif JucePlugin_IsSynth
    if (isInput || (numChannels != 1 && numChannels != 2))
        return false;
   #else
    if (numChannels != 1 && numChannels != 2)
        return false;

    if (! AudioProcessor::setPreferredBusArrangement (! isInput, bus, preferredSet))
        return false;
   #endif

    return AudioProcessor::setPreferredBusArrangement (isInput, bus, preferredSet);
}
#endif

void AudioWatermarkingAudioProcessor::silenceBuffer(AudioSampleBuffer& buffer)
{
	for(int i =0; i < buffer.getNumChannels(); i++)
	{
		float* bufferData = buffer.getWritePointer(i);

		for(int j =0; j < buffer.getNumSamples(); j++)
		{
			bufferData[j] = 0;
		}
	}
}

void AudioWatermarkingAudioProcessor::performWMaction(AudioSampleBuffer& buffer)
{
	AudioPlayHead::CurrentPositionInfo info;
	this->getPlayHead()->getCurrentPosition(info);

	if(info.isPlaying || this->isNonRealtime())
	{
		//Only fill our buffer if there we're actually doing something.
		_bufferTools.bufferAccumulator(buffer);
	}

	if (_bufferTools.isBuffer1Full)
	{
		if(this->isNonRealtime())
		{
			//Offline, so add a watermark.
			_watermark.addWaterMark(_bufferTools.accumulatedBuffer1);
		}
		else if(info.isPlaying)
		{
			//Online, so we're only doing detection.
			_watermark.detectWaterMark(_bufferTools.accumulatedBuffer1);
		}

		//Buffer has been filled and processed, so reset it.
		_bufferTools.isBuffer1Full = false;
		_bufferTools.isBufferEmpty = false;
		_bufferTools.buffer1Index = 0;		
	}
	else if(_bufferTools.isBuffer2Full)
	{
		if(this->isNonRealtime())
		{
			//Offline, so add a watermark.
			_watermark.addWaterMark(_bufferTools.accumulatedBuffer2);
		}
		else if(info.isPlaying)
		{
			//Online, so we're only doing detection.
			_watermark.detectWaterMark(_bufferTools.accumulatedBuffer2);
		}

		//Buffer has been filled and processed, so reset it.
		_bufferTools.isBuffer2Full = false;
		_bufferTools.isBufferEmpty = false;
		_bufferTools.buffer2Index = 0;	
	}


	if(!_bufferTools.isBufferEmpty)
	{
		//Buffer is full, so now we can start writing to the output.
		_bufferTools.bufferAccumulatorOut(buffer);				
	}	
	else
	{
		//Silence until buffer is full.
		silenceBuffer(buffer);
	}

	if(info.isPlaying && !this->isNonRealtime())
	{
		//Silence audio if we're detecting watermarks
		silenceBuffer(buffer);
	}
	else if(!info.isPlaying)
	{
		silenceBuffer(buffer);
	}
}

void AudioWatermarkingAudioProcessor::processBlock (AudioSampleBuffer& buffer, MidiBuffer& midiMessages)
{
    const int totalNumInputChannels  = getTotalNumInputChannels();
    const int totalNumOutputChannels = getTotalNumOutputChannels();
	

	if(buffer.getNumSamples() > NFREQ)
	{
		int count = 0;
		while(count != -1)
		{
			AudioSampleBuffer tempBuffer;
			tempBuffer.setSize(buffer.getNumChannels(), NFREQ, 0, 1, 1);
			
			if((count * NFREQ) + NFREQ < buffer.getNumSamples())			
			{
				for(int ch =0; ch < buffer.getNumChannels(); ch++)
				{
					tempBuffer.copyFrom(ch, 0, buffer, ch, count * NFREQ, NFREQ);
				}

				performWMaction(tempBuffer);

				for(int ch =0; ch < buffer.getNumChannels(); ch++)
				{
					buffer.copyFrom(ch, count * NFREQ, tempBuffer, ch, 0, NFREQ);
				}

				count++;

			}
			else if((count * NFREQ) + NFREQ == buffer.getNumSamples())			
			{
				for(int ch =0; ch < buffer.getNumChannels(); ch++)
				{
					tempBuffer.copyFrom(ch, 0, buffer, ch, count * NFREQ, NFREQ);
				}

				performWMaction(tempBuffer);

				for(int ch =0; ch < buffer.getNumChannels(); ch++)
				{
					buffer.copyFrom(ch, count * NFREQ, tempBuffer, ch, 0, NFREQ);
				}

				count = -1;
			}
			else
			{
				tempBuffer.setSize(buffer.getNumChannels(), ((count * NFREQ) + NFREQ - buffer.getNumSamples()), 0, 1, 1);

				for(int ch =0; ch < buffer.getNumChannels(); ch++)
				{
					tempBuffer.copyFrom(ch, 0, buffer, ch, count * NFREQ, ((count * NFREQ) + NFREQ - buffer.getNumSamples()));
				}

				performWMaction(tempBuffer);

				for(int ch =0; ch < buffer.getNumChannels(); ch++)
				{
					buffer.copyFrom(ch, count * NFREQ, tempBuffer, ch, 0, ((count * NFREQ) + NFREQ - buffer.getNumSamples()));
				}

				count = -1;
			}

		}
	}
	else
	{
		performWMaction(buffer);
	}
	
}

//==============================================================================
bool AudioWatermarkingAudioProcessor::hasEditor() const
{
    return true; // (change this to false if you choose to not supply an editor)
}

AudioProcessorEditor* AudioWatermarkingAudioProcessor::createEditor()
{
    return new AudioWatermarkingAudioProcessorEditor (*this);
}

//==============================================================================
void AudioWatermarkingAudioProcessor::getStateInformation (MemoryBlock& destData)
{
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
}

void AudioWatermarkingAudioProcessor::setStateInformation (const void* data, int sizeInBytes)
{
    // You should use this method to restore your parameters from this memory block,
    // whose contents will have been created by the getStateInformation() call.
}

void AudioWatermarkingAudioProcessor::setWatermark(std::string watermark)
{
	_watermarkToApply = _watermark.stringToHex(watermark);
	_watermark.readWMBits(_watermarkToApply);
}

std::string AudioWatermarkingAudioProcessor::getWatermark()
{
	return _watermark.getUIUpdate();
}

//==============================================================================
// This creates new instances of the plugin..
AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new AudioWatermarkingAudioProcessor();
}
