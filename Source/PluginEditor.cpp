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

//[Headers] You can add your own extra header files here...
#include "PluginProcessor.h"
//[/Headers]

#include "PluginEditor.h"



//[MiscUserDefs] You can add your own user definitions and misc code here...
//[/MiscUserDefs]

//==============================================================================
AudioWatermarkingAudioProcessorEditor::AudioWatermarkingAudioProcessorEditor (AudioProcessor& p)
    : AudioProcessorEditor(&p)
{
    //[Constructor_pre] You can add your own custom stuff here..
    //[/Constructor_pre]

    addAndMakeVisible (_watermarkingLabel = new Label ("new label",
                                                       TRANS("Audio Watermarking Tool\n"
                                                       "\n")));
    _watermarkingLabel->setFont (Font (41.10f, Font::plain));
    _watermarkingLabel->setJustificationType (Justification::centredLeft);
    _watermarkingLabel->setEditable (false, false, false);
    _watermarkingLabel->setColour (Label::textColourId, Colour (0xffdd8282));
    _watermarkingLabel->setColour (TextEditor::textColourId, Colours::black);
    _watermarkingLabel->setColour (TextEditor::backgroundColourId, Colour (0x00000000));

    addAndMakeVisible (_hexWatermarkCode = new TextEditor ("new text editor"));
    _hexWatermarkCode->setMultiLine (false);
    _hexWatermarkCode->setReturnKeyStartsNewLine (false);
    _hexWatermarkCode->setReadOnly (false);
    _hexWatermarkCode->setScrollbarsShown (true);
    _hexWatermarkCode->setCaretVisible (true);
    _hexWatermarkCode->setPopupMenuEnabled (true);
    _hexWatermarkCode->setText (String());

    addAndMakeVisible (_applyWatermarkButton = new TextButton ("new button"));
    _applyWatermarkButton->setButtonText (TRANS("Apply"));
    _applyWatermarkButton->addListener (this);

    addAndMakeVisible (labelWaterMarkFound = new Label ("new label",
                                                        TRANS("Watermark in current audio file:")));
    labelWaterMarkFound->setFont (Font (18.30f, Font::plain));
    labelWaterMarkFound->setJustificationType (Justification::centredLeft);
    labelWaterMarkFound->setEditable (false, false, false);
    labelWaterMarkFound->setColour (Label::textColourId, Colour (0xffa33131));
    labelWaterMarkFound->setColour (TextEditor::textColourId, Colours::black);
    labelWaterMarkFound->setColour (TextEditor::backgroundColourId, Colour (0x00000000));

    addAndMakeVisible (textEditorWaterMarkFound = new TextEditor ("new text editor"));
    textEditorWaterMarkFound->setMultiLine (false);
    textEditorWaterMarkFound->setReturnKeyStartsNewLine (false);
    textEditorWaterMarkFound->setReadOnly (true);
    textEditorWaterMarkFound->setScrollbarsShown (true);
    textEditorWaterMarkFound->setCaretVisible (false);
    textEditorWaterMarkFound->setPopupMenuEnabled (true);
    textEditorWaterMarkFound->setText (String());


    //[UserPreSize]
    //[/UserPreSize]

    setSize (460, 400);


    //[Constructor] You can add your own custom stuff here..
	startTimer(200);//starts timer with interval of 200mS
    //[/Constructor]
}

AudioWatermarkingAudioProcessorEditor::~AudioWatermarkingAudioProcessorEditor()
{
    //[Destructor_pre]. You can add your own custom destruction code here..
    //[/Destructor_pre]

    _watermarkingLabel = nullptr;
    _hexWatermarkCode = nullptr;
    _applyWatermarkButton = nullptr;
    labelWaterMarkFound = nullptr;
    textEditorWaterMarkFound = nullptr;


    //[Destructor]. You can add your own custom destruction code here..
    //[/Destructor]
}

//==============================================================================
void AudioWatermarkingAudioProcessorEditor::paint (Graphics& g)
{
    //[UserPrePaint] Add your own custom painting code here..
    //[/UserPrePaint]

    g.fillAll (Colour (0xff3c3333));

    //[UserPaint] Add your own custom painting code here..
    //[/UserPaint]
}

void AudioWatermarkingAudioProcessorEditor::resized()
{
    //[UserPreResize] Add your own custom resize code here..
    //[/UserPreResize]

    _watermarkingLabel->setBounds (8, 8, 456, 112);
    _hexWatermarkCode->setBounds (24, 144, 408, 24);
    _applyWatermarkButton->setBounds (152, 176, 150, 24);
    labelWaterMarkFound->setBounds (24, 224, 248, 24);
    textEditorWaterMarkFound->setBounds (280, 224, 150, 24);
    //[UserResized] Add your own custom resize handling here..
    //[/UserResized]
}

void AudioWatermarkingAudioProcessorEditor::buttonClicked (Button* buttonThatWasClicked)
{
    //[UserbuttonClicked_Pre]
    //[/UserbuttonClicked_Pre]

    if (buttonThatWasClicked == _applyWatermarkButton)
    {
        //[UserButtonCode__applyWatermarkButton] -- add your button handler code here..
		AudioWatermarkingAudioProcessor* ap = static_cast<AudioWatermarkingAudioProcessor*>(getAudioProcessor());
		ap->setWatermark(_hexWatermarkCode->getText().toStdString());

        //[/UserButtonCode__applyWatermarkButton]
    }

    //[UserbuttonClicked_Post]
    //[/UserbuttonClicked_Post]
}



//[MiscUserCode] You can add your own definitions of your custom methods or any other code here...
void AudioWatermarkingAudioProcessorEditor::timerCallback()
{
	AudioProcessor* ap = getAudioProcessor();

	AudioWatermarkingAudioProcessor* apWM = static_cast<AudioWatermarkingAudioProcessor*>(getAudioProcessor());

	AudioPlayHead* playHead = ap->getPlayHead();

	AudioPlayHead::CurrentPositionInfo pos;
	playHead->getCurrentPosition(pos);

	if (pos.isRecording || pos.isPlaying || pos.isLooping)
	{
		_applyWatermarkButton->setEnabled(false);
		_hexWatermarkCode->setEnabled(false);
		textEditorWaterMarkFound->clear();
		textEditorWaterMarkFound->setText(apWM->getWatermark(), dontSendNotification);
	}
	else
	{
		_applyWatermarkButton->setEnabled(true);
		_hexWatermarkCode->setEnabled(true);
	}

}

//[/MiscUserCode]


//==============================================================================
#if 0
/*  -- Projucer information section --

    This is where the Projucer stores the metadata that describe this GUI layout, so
    make changes in here at your peril!

BEGIN_JUCER_METADATA

<JUCER_COMPONENT documentType="Component" className="AudioWatermarkingAudioProcessorEditor"
                 componentName="" parentClasses="public AudioProcessorEditor, public Timer"
                 constructorParams="AudioProcessor&amp; p" variableInitialisers="AudioProcessorEditor(&amp;p)"
                 snapPixels="8" snapActive="1" snapShown="1" overlayOpacity="0.330"
                 fixedSize="0" initialWidth="460" initialHeight="400">
  <BACKGROUND backgroundColour="ff3c3333"/>
  <LABEL name="new label" id="98334ec62c7189ea" memberName="_watermarkingLabel"
         virtualName="" explicitFocusOrder="0" pos="8 8 456 112" textCol="ffdd8282"
         edTextCol="ff000000" edBkgCol="0" labelText="Audio Watermarking Tool&#10;&#10;"
         editableSingleClick="0" editableDoubleClick="0" focusDiscardsChanges="0"
         fontname="Default font" fontsize="41.100000000000001421" bold="0"
         italic="0" justification="33"/>
  <TEXTEDITOR name="new text editor" id="787f7a8f1f45a064" memberName="_hexWatermarkCode"
              virtualName="" explicitFocusOrder="0" pos="24 144 408 24" initialText=""
              multiline="0" retKeyStartsLine="0" readonly="0" scrollbars="1"
              caret="1" popupmenu="1"/>
  <TEXTBUTTON name="new button" id="f22159baeccef462" memberName="_applyWatermarkButton"
              virtualName="" explicitFocusOrder="0" pos="152 176 150 24" buttonText="Apply"
              connectedEdges="0" needsCallback="1" radioGroupId="0"/>
  <LABEL name="new label" id="a88033defa612e90" memberName="labelWaterMarkFound"
         virtualName="" explicitFocusOrder="0" pos="24 224 248 24" textCol="ffa33131"
         edTextCol="ff000000" edBkgCol="0" labelText="Watermark in current audio file:"
         editableSingleClick="0" editableDoubleClick="0" focusDiscardsChanges="0"
         fontname="Default font" fontsize="18.300000000000000711" bold="0"
         italic="0" justification="33"/>
  <TEXTEDITOR name="new text editor" id="7d6b6e53eae43ace" memberName="textEditorWaterMarkFound"
              virtualName="" explicitFocusOrder="0" pos="280 224 150 24" initialText=""
              multiline="0" retKeyStartsLine="0" readonly="1" scrollbars="1"
              caret="0" popupmenu="1"/>
</JUCER_COMPONENT>

END_JUCER_METADATA
*/
#endif


//[EndFile] You can add extra defines here...
//[/EndFile]
