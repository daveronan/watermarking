/*
  ==============================================================================

    This file was auto-generated!

    It contains the basic framework code for a JUCE plugin processor.

  ==============================================================================
*/

#ifndef PLUGINPROCESSOR_H_INCLUDED
#define PLUGINPROCESSOR_H_INCLUDED

#include <JuceHeader.h>
#include "watermark.h"
#include "BufferTools.h"



//==============================================================================
/**
*/
class AudioWatermarkingAudioProcessor  : public AudioProcessor
{
public:
    //==============================================================================
    AudioWatermarkingAudioProcessor();
    ~AudioWatermarkingAudioProcessor();

    //==============================================================================
    void prepareToPlay (double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
	void silenceBuffer(AudioSampleBuffer& buffer);
	void performWMaction(AudioSampleBuffer& buffer);
#ifndef JucePlugin_PreferredChannelConfigurations
    bool setPreferredBusArrangement (bool isInput, int bus, const AudioChannelSet& preferredSet) override;
   #endif

    void processBlock (AudioSampleBuffer&, MidiBuffer&) override;

    //==============================================================================
    AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram (int index) override;
    const String getProgramName (int index) override;
    void changeProgramName (int index, const String& newName) override;

    //==============================================================================
    void getStateInformation (MemoryBlock& destData) override;
    void setStateInformation (const void* data, int sizeInBytes) override;
	
	void setWatermark(std::string watermark);
	std::string getWatermark();

private:
    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (AudioWatermarkingAudioProcessor)
		
		//std::string UserParams[totalNumParam];		
		//bool UIUpdateFlag;
		int _samplesPerBlock;	
		Watermark _watermark;
		BufferTools _bufferTools;
		std::string _watermarkToApply;

};


#endif  // PLUGINPROCESSOR_H_INCLUDED
