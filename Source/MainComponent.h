/*
  ==============================================================================

    David Temming
    Elon University

  ==============================================================================
*/

#pragma once

#include "../JuceLibraryCode/JuceHeader.h"

//==============================================================================
/*
    This component lives inside our window, and this is where you should put all
    your controls and content.
*/
class MainComponent   : public AudioAppComponent,
                        public ChangeListener,
                        private Timer
{
public:
    //==============================================================================
    MainComponent();
    ~MainComponent();

    //==============================================================================
    void prepareToPlay (int samplesPerBlockExpected, double sampleRate) override;
    void getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill) override;
    void releaseResources() override;

    //==============================================================================
    void paint (Graphics& g) override;
    void resized() override;

private:
    //==============================================================================
    // PRIVATE MEMBER VARIABLES
    enum TransportState
    {
        Stopped,
        Starting,
        originalPlaying,
        additivePlaying,
        Pausing,
        Paused,
        Stopping
    };
    
    TransportState state;
    
    void openButtonClicked();
    void processButtonClicked();
    void playButtonClicked();
    void stopButtonClicked();
    void additivePlayButtonClicked();
    void additiveStopButtonClicked();
    
    void averageNormalize(AudioBuffer<float>* input, AudioBuffer<float>* output);
    void peakNormalize(AudioBuffer<float>* input, AudioBuffer<float>* output);
    void transportStateChanged(TransportState newState, AudioTransportSource* transport, TextButton* play, TextButton* stop);
    void thumbNailChanged();
    void paintIfNoFileLoaded(Graphics& g, const Rectangle<int>& thumbnailBounds);
    void paintIfFileLoaded(Graphics& g, const Rectangle<int>& thumbnailBounds);
    void paintIfNotResynthesized(Graphics& g, const Rectangle<int>& thumbnailBounds);
    void paintIfResynthesized(Graphics& g, const Rectangle<int>& thumbnailBounds);
    void timerCallback();
    void sumOfSquareDiff(AudioBuffer<float>* originalSample, AudioBuffer<float>* additiveSample);
    void writeFile(AudioBuffer<float> buffer);
    void printBuffer(AudioBuffer<float>* buffer);
    
    //virtual function that has to be implemented when ChangeListener is inherited
    void changeListenerCallback (ChangeBroadcaster *source) override;
    
    AudioFormatManager formatManager;
    std::unique_ptr<AudioFormatReaderSource> readerSource;
    AudioTransportSource originalTransport;
    AudioThumbnailCache inputThumbnailCache;
    AudioThumbnail inputThumbnail;
    AudioThumbnailCache outputThumbnailCache;
    AudioThumbnail outputThumbnail;
    bool resynthesized;
    String fileName;
    AudioBuffer<float> *inputSampleBuffer = NULL;
    AudioBuffer<float> *outputSampleBuffer = NULL;
    
    TextButton openButton;
    TextButton processButton;
    TextButton playButton;
    TextButton stopButton;
    TextButton additivePlayButton;
    TextButton additiveStopButton;
    Label fileNameDisplay;
    
    int countOfSamples;
    
    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR (MainComponent)
};
