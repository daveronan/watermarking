/*
  ==============================================================================

  This is an automatically generated GUI class created by the Projucer!

  Be careful when adding custom code to these files, as only the code within
  the "//[xyz]" and "//[/xyz]" sections will be retained when the file is loaded
  and re-saved.

  Created with Projucer version: 4.2.1

  ------------------------------------------------------------------------------

  The Projucer is part of the JUCE library - "Jules' Utility Class Extensions"
  Copyright (c) 2015 - ROLI Ltd.

  ==============================================================================
*/

#ifndef __JUCE_HEADER_5391EE8143A66EAD__
#define __JUCE_HEADER_5391EE8143A66EAD__

//[Headers]     -- You can add your own extra header files here --
#include "../JuceLibraryCode/JuceHeader.h"
//[/Headers]



//==============================================================================
/**
                                                                    //[Comments]
    An auto-generated component, created by the Projucer.

    Describe your class and how it works here!
                                                                    //[/Comments]
*/
class AudioWatermarkingAudioProcessorEditor  : public AudioProcessorEditor,
                                               public Timer,
                                               public ButtonListener
{
public:
    //==============================================================================
    AudioWatermarkingAudioProcessorEditor (AudioProcessor& p);
    ~AudioWatermarkingAudioProcessorEditor();

    //==============================================================================
    //[UserMethods]     -- You can add your own custom methods in this section.
	void timerCallback();
    //[/UserMethods]

    void paint (Graphics& g) override;
    void resized() override;
    void buttonClicked (Button* buttonThatWasClicked) override;



private:
    //[UserVariables]   -- You can add your own custom variables in this section.

    //[/UserVariables]

    //==============================================================================
    ScopedPointer<Label> _watermarkingLabel;
    ScopedPointer<TextEditor> _hexWatermarkCode;
    ScopedPointer<TextButton> _applyWatermarkButton;
    ScopedPointer<Label> labelWaterMarkFound;
    ScopedPointer<TextEditor> textEditorWaterMarkFound;


    //==============================================================================
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (AudioWatermarkingAudioProcessorEditor)
};

//[EndFile] You can add extra defines here...
//[/EndFile]

#endif   // __JUCE_HEADER_5391EE8143A66EAD__
