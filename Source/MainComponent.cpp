/*
  ==============================================================================

    David Temming
    Elon University

    References:
    Hussain, Rafat. C Implementation of Wavelet Transform (DWT,SWT and MODWT), October 22, 2019. [Online], Available: https://github.com/rafat/wavelib

  ==============================================================================
*/

#include "MainComponent.h"
#include "wavelib.h"

#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>



using namespace std;

//==============================================================================
MainComponent::MainComponent() : state(Stopped), inputThumbnailCache(5), inputThumbnail(512, formatManager, inputThumbnailCache), outputThumbnailCache(5), outputThumbnail(512, formatManager, outputThumbnailCache),  resynthesized(false), fileName(NULL), openButton("Open"), processButton("Resynthesize"), playButton("Play"),  stopButton("Stop"),additivePlayButton("Play (disabled)"),additiveStopButton("Stop (disabled)"), fileNameDisplay("No file selected"), countOfSamples(0)
{
    // Make sure you set the size of the component after
    // you add any child components.
    setSize (1000, 330);
    
    openButton.onClick = [this] { openButtonClicked(); };
    addAndMakeVisible(&openButton);
    
    playButton.onClick = [this] {playButtonClicked(); };
    playButton.setColour(TextButton::buttonColourId, Colours::green);
    playButton.setEnabled(true);
    addAndMakeVisible(&playButton);
    
    stopButton.onClick = [this] {stopButtonClicked(); };
    stopButton.setColour(TextButton::buttonColourId, Colours::red);
    stopButton.setEnabled(false);
    addAndMakeVisible(&stopButton);
    
    additivePlayButton.onClick = [this] {additivePlayButtonClicked(); };
    additivePlayButton.setColour(TextButton::buttonColourId, Colours::green);
    additivePlayButton.setEnabled(true);
    addAndMakeVisible(&additivePlayButton);
    
    additiveStopButton.onClick = [this] {additiveStopButtonClicked(); };
    additiveStopButton.setColour(TextButton::buttonColourId, Colours::red);
    additiveStopButton.setEnabled(false);
    addAndMakeVisible(&additiveStopButton);
    
    processButton.onClick = [this] {processButtonClicked(); };
    processButton.setColour(TextButton::buttonColourId, Colours::purple);
    processButton.setEnabled(true);
    addAndMakeVisible(&processButton);
    
    fileNameDisplay.setText("No file selected", dontSendNotification);
    fileNameDisplay.setEnabled(true);
    addAndMakeVisible(&fileNameDisplay);
    
    //allows load of .wav and .aiff files
    formatManager.registerBasicFormats();
    
    //adding change listeners to the transports and thumbnails
    originalTransport.addChangeListener(this);
    inputThumbnail.addChangeListener(this);
    outputThumbnail.addChangeListener(this);
    
    // Some platforms require permissions to open input channels so request that here
    if (RuntimePermissions::isRequired (RuntimePermissions::recordAudio)
        && ! RuntimePermissions::isGranted (RuntimePermissions::recordAudio))
    {
        RuntimePermissions::request (RuntimePermissions::recordAudio,
                                     [&] (bool granted) { if (granted)  setAudioChannels (2, 2); });
    }
    else
    {
        // Specifying the number of input and output channels to open
        setAudioChannels (2, 2);
    }
}

MainComponent::~MainComponent()
{
    // This shuts down the audio device and clears the audio source.
    delete inputSampleBuffer;
    delete outputSampleBuffer;
    shutdownAudio();
}

//==============================================================================
void MainComponent::prepareToPlay (int samplesPerBlockExpected, double sampleRate)
{
    //called on the audio thread, not the GUI thread.
    originalTransport.prepareToPlay(samplesPerBlockExpected, sampleRate);
}

void MainComponent::timerCallback()
{
    repaint();
}

void MainComponent::sumOfSquareDiff(AudioBuffer<float>* originalSample, AudioBuffer<float>* additiveSample)
{
    //initialize the sum at zero
    double sumLeft = 0;
    double sumRight = 0;
    
    //calculating sum
    for(int i = 0; i < originalSample->getNumSamples(); i++){
        double originalLeft = originalSample->getSample(0, i);
        double additiveLeft = additiveSample->getSample(0, i);
        double originalRight = originalSample->getSample(1, i);
        double additiveRight = additiveSample->getSample(1, i);

        sumLeft = sumLeft + pow((originalLeft - additiveLeft), 2);
        sumRight = sumRight + pow((originalRight - additiveRight), 2);
    }
    
    //printing SSD results to the console
    DBG("Sum of Squared Differences in Left Channel = " + std::to_string(sumLeft));
    DBG("Sum of Squared Differences in Right Channel = " + std::to_string(sumRight));
    DBG("Sum of Squared Differences in Both Channels = " + std::to_string(sumLeft + sumRight));
}

void MainComponent::openButtonClicked()
{
    //file selection
    FileChooser chooser("Choose a .wav file", File("/Users/davidtemming/ResearchTestFiles"), "*.wav; *.aiff", true, false);

    //if the user chooses a file
    if (chooser.browseForFileToOpen())
    {
        File myFile;
        myFile = chooser.getResult();
        
        //reading the file
        AudioFormatReader* reader = formatManager.createReaderFor(myFile);

        if(reader != nullptr)
        {
            //set the fileNameDisplay label name
            fileName = myFile.getFileName();
            fileNameDisplay.setText(fileName, dontSendNotification);
            
            //preparing the file to play
            std::unique_ptr<AudioFormatReaderSource> tempSource (new AudioFormatReaderSource (reader, true));
            
            //initializing the inputSampleBuffer
            if(inputSampleBuffer != NULL){
                DBG("ABOUT TO DELETE INPUT SAMPLE BUFFER");
                delete inputSampleBuffer;
            }
            inputSampleBuffer = new AudioBuffer<float>(2, (int)reader->lengthInSamples);
            inputSampleBuffer->clear();
            
            //reading the data into the input sample buffer
            reader->read (inputSampleBuffer,
                          0,
                          (int) reader->lengthInSamples,
                          0,
                          true,
                          true);
            
            originalTransport.setSource(tempSource.get());
            transportStateChanged(Stopped, &originalTransport, &playButton, &stopButton);
            
            //set the inputThumbnail source
            inputThumbnail.setSource(new FileInputSource(myFile));
            
            //clear the outputThumbnail source
            if(outputSampleBuffer != NULL){
                delete outputSampleBuffer;
            }
            outputThumbnail.clear();
            resynthesized = false;
        
            //release tempSource
            readerSource.reset(tempSource.release());
        }
    }
}

void MainComponent::printBuffer(AudioBuffer<float>* buffer){
        DBG("");
        for(int i = 0; i < buffer->getNumSamples(); i++){
            DBG("LEFT SAMPLE" + std::to_string(i) + ": " + std::to_string(inputSampleBuffer->getSample(0, i)));
            DBG("RIGHT SAMPLE" + std::to_string(i) + ": " + std::to_string(inputSampleBuffer->getSample(1, i)));
        }
}

void MainComponent::processButtonClicked()
{
    
    //Setting wavelib parameters
    int i, N, J, subscale, a0, iter, nd, k;
    double *inpLeft, *oupLeft, *inpRight, *oupRight;
    double dt, dj, s0, param, mn, mn2;
    double td, tn, den, num, recon_mean, recon_var;
    cwt_object wt;
    cwt_object wt2;

    char const *wave = "morlet";
    char const *type = "pow";

    N = inputSampleBuffer->getNumSamples();
    param = 6.0;
    subscale = 4;
    dt = .25;
    //s0 = .0032;
    s0 = dt; // initial scale value
    dj = 1.0 / (double)subscale;
    J = 11 * subscale; // Total Number of scales
    a0 = 2;//power
    
    i = 0;
    
    //creating cwt objects
    wt = cwt_init(wave, param, N, dt, J);
    wt2 = cwt_init(wave, param, N, dt, J);
    
    //allocating space
    inpLeft = (double*)malloc(sizeof(double)* N);
    oupLeft = (double*)malloc(sizeof(double)* N);
    inpRight = (double*)malloc(sizeof(double)* N);
    oupRight = (double*)malloc(sizeof(double)* N);
    
    //filling wavelib input ararys from JUCE input buffers
    for (i = 0; i < N; i++) {
        inpLeft[i] = inputSampleBuffer->getSample(0, i);
        inpRight[i] = inputSampleBuffer->getSample(1, i);
    }
    
    //setting CWT scales
    setCWTScales(wt, s0, dj, type, a0);
    setCWTScales(wt2, s0, dj, type, a0);
    
    //performing the CWT
    cwt(wt, inpLeft);
    cwt(wt2, inpRight);
    
    //WAVELIB SUMMARY OUTPUT
    printf("\n MEAN LEFT%g \n", wt->smean);
    printf("\n MEAN RIGHT%g \n", wt2->smean);
    
    mn = 0.0;
    mn2 = 0.0;
    
    for (i = 0; i < N; i++) {
        mn += sqrt(wt->output[i].re * wt->output[i].re + wt->output[i].im * wt->output[i].im);
        mn2 += sqrt(wt2->output[i].re * wt2->output[i].re + wt2->output[i].im * wt2->output[i].im);
    }
    
    //LEFT CHANNEL
    DBG("LEFT SUMMARY");
    cwt_summary(wt);
    printf("\n abs mean %g \n", mn / N);
    
    printf("\n\n");
    printf("Let CWT w = w(j, n/2 - 1) where n = %d\n\n", N);
    nd = N/2 - 1;
    
    printf("%-15s%-15s%-15s%-15s \n","j","Scale","Period","ABS(w)^2");
    for(k = 0; k < wt->J;++k) {
        iter = nd + k * N;
        printf("%-15d%-15lf%-15lf%-15lf \n",k,wt->scale[k],wt->period[k],
               wt->output[iter].re * wt->output[iter].re + wt->output[iter].im * wt->output[iter].im);
    }
    
    //RIGHT CHANNEL
    DBG("RIGHT SUMMARY");
    cwt_summary(wt2);
    printf("\n abs mean %g \n", mn2 / N);
    
    printf("\n\n");
    printf("Let CWT w = w(j, n/2 - 1) where n = %d\n\n", N);
    nd = N/2 - 1;
    
    printf("%-15s%-15s%-15s%-15s \n","j","Scale","Period","ABS(w)^2");
    for(k = 0; k < wt->J;++k) {
        iter = nd + k * N;
        printf("%-15d%-15lf%-15lf%-15lf \n",k,wt2->scale[k],wt2->period[k],
               wt2->output[iter].re * wt2->output[iter].re + wt2->output[iter].im * wt2->output[iter].im);
    }
    
    //LEFT CHANNEL
    icwt(wt, oupLeft);
    
    num = den = recon_var = recon_mean = 0.0;
    printf("\n\n");
    printf("Signal Reconstruction\n");
    printf("%-15s%-15s%-15s \n","i","Input(i)","Output(i)");

    for (i = N - 10; i < N; i++) {
        printf("%-15d%-15lf%-15lf \n", i,inpLeft[i] , oupLeft[i]);
    }

    for (i = 0; i < N; i++) {
        //printf("%g %g \n", oup[i] ,inp[i] - wt->smean);
        td = inpLeft[i] ;
        tn = oupLeft[i] - td;
        num += (tn * tn);
        den += (td * td);
        recon_mean += oupLeft[i];
    }

    recon_var = sqrt(num / N);
    recon_mean /= N;

    printf("\nRMS Error %g \n", sqrt(num) / sqrt(den));
    printf("\nVariance %g \n", recon_var);
    printf("\nMean %g \n", recon_mean);
    
    //RIGHT CHANNEL COPY
    icwt(wt2, oupRight);
    num = den = recon_var = recon_mean = 0.0;
    printf("\n\n");
    printf("Signal Reconstruction\n");
    printf("%-15s%-15s%-15s \n","i","Input(i)","Output(i)");

    for (i = N - 10; i < N; i++) {
        printf("%-15d%-15lf%-15lf \n", i,inpRight[i] , oupRight[i]);
    }

    for (i = 0; i < N; i++) {
        //printf("%g %g \n", oup[i] ,inp[i] - wt->smean);
        td = inpRight[i];
        tn = oupRight[i] - td;
        num += (tn * tn);
        den += (td * td);
        recon_mean += oupRight[i];
    }

    recon_var = sqrt(num / N);
    recon_mean /= N;

    printf("\nRMS Error %g \n", sqrt(num) / sqrt(den));
    printf("\nVariance %g \n", recon_var);
    printf("\nMean %g \n", recon_mean);
    //END OF WAVELIB SUMMARY OUTPUT
    
    //creating 2d arrays for the modulus and phase
    double modulusLeft, modulusRight;
    double phaseLeft, phaseRight;
    double realLeft, realRight;
    double imaginaryLeft, imaginaryRight;
    
    double** modulusDataLeft = new double*[J];
    double** phaseDataLeft = new double*[J];
    double** modulusDataRight = new double*[J];
    double** phaseDataRight = new double*[J];
    
    for (int scales = 0; scales < J; scales++) {
        modulusDataLeft[scales] = new double[N];
        phaseDataLeft[scales] = new double[N];
        modulusDataRight[scales] = new double[N];
        phaseDataRight[scales] = new double[N];
    }
    
    //filling the 2d arrays for modulus and phase
    for (int scales = 0; scales < J; scales++) {
        for (int sampleInScale = 0; sampleInScale < N; sampleInScale++) {
            
            //Compute the imaginary and real parts
            realLeft = wt->output[(scales * N) + sampleInScale].re;
            imaginaryLeft = wt->output[(scales * N) + sampleInScale].im;
            realRight = wt2->output[(scales * N) + sampleInScale].re;
            imaginaryRight = wt2->output[(scales * N) + sampleInScale].im;
            
            // Compute the modulus and phase
            modulusLeft = sqrt(pow(realLeft, 2) + pow(imaginaryLeft, 2));
            modulusRight = sqrt(pow(realRight, 2) + pow(imaginaryRight, 2));
            //phase = atan(imaginary / real);
            phaseLeft = atan2(imaginaryLeft, realLeft);
            phaseRight = atan2(imaginaryRight, realRight);
            
            // Place the modulus and phase into the appropriate slot of the data array
            modulusDataLeft[scales][sampleInScale] = modulusLeft;
            phaseDataLeft[scales][sampleInScale] = phaseLeft;
            modulusDataRight[scales][sampleInScale] = modulusRight;
            phaseDataRight[scales][sampleInScale] = phaseRight;
        }
    }
    
    //creating buffer for the output
    outputSampleBuffer = new AudioBuffer<float>(2, inputSampleBuffer->getNumSamples());
    outputSampleBuffer->clear();
    
    //computing final amplitudes and filling the output buffer
    for(int j = 0; j < J; j++){
        for (i = 0; i < N; i++) {
            double currentSampleLeft = outputSampleBuffer->getSample(0, i);
            double currentSampleRight = outputSampleBuffer->getSample(1, i);

            if(j == 0){
                currentSampleLeft = 0.0;
                currentSampleRight = 0.0;
            }
            currentSampleLeft += modulusDataLeft[j][i] * cos(phaseDataLeft[j][i]);
            currentSampleRight += modulusDataRight[j][i] * cos(phaseDataRight[j][i]);
            
            outputSampleBuffer->setSample(0, i, currentSampleLeft);
            outputSampleBuffer->setSample(1, i, currentSampleRight);
        }
    }
    
    for(int i = 0; i < N; i++){
        outputSampleBuffer->setSample(0, i, outputSampleBuffer->getSample(0, i));
        outputSampleBuffer->setSample(1, i, outputSampleBuffer->getSample(1, i));

    }
    
    averageNormalize(inputSampleBuffer, outputSampleBuffer);
    //peakNormalize(inputSampleBuffer, outputSampleBuffer);
    
    //freeing objects
    free(inpLeft);
    free(oupLeft);
    free(inpRight);
    free(oupRight);
    cwt_free(wt);
    cwt_free(wt2);
    
    //performing SSD on the input and output
    sumOfSquareDiff(inputSampleBuffer, outputSampleBuffer);
    
    DBG("number of samples before writing the file: " + std::to_string(outputSampleBuffer->getNumSamples()));
    
    //writing out a file
    writeFile(*outputSampleBuffer);
    resynthesized = true;
    
    // Clean up
    for (int scales = 0; scales < J; scales++) {
        delete[] modulusDataLeft[scales];
        delete[] phaseDataLeft[scales];
        delete[] modulusDataRight[scales];
        delete[] phaseDataRight[scales];
    }
    delete[] modulusDataLeft;
    delete[] phaseDataLeft;
    delete[] modulusDataRight;
    delete[] phaseDataRight;
}

void MainComponent::averageNormalize(AudioBuffer<float>* input, AudioBuffer<float>* output){
    DBG("**** AVERAGE NORMALIZATION LEFT AND RIGHT");
    double avgLeftInput = 0.0;
    double avgLeftOutput = 0.0;
    double avgRightInput = 0.0;
    double avgRightOutput = 0.0;
    
    for(int i = 0; i < inputSampleBuffer->getNumSamples(); i++){
        avgLeftInput += abs(input->getSample(0, i));
        avgLeftOutput += abs(output->getSample(0, i));
        avgRightInput += abs(input->getSample(1, i));
        avgRightOutput += abs(output->getSample(1, i));
    }
    avgLeftInput = avgLeftInput / input->getNumSamples();
    avgLeftOutput = avgLeftOutput /  input->getNumSamples();
    avgRightInput = avgRightInput / input->getNumSamples();
    avgRightOutput = avgRightOutput /  input->getNumSamples();
    
    DBG("AVGS");
    DBG(avgLeftInput);
    DBG(avgLeftOutput);
    
    for(int i = 0; i < input->getNumSamples(); i++){
        output->setSample(0, i, output->getSample(0, i) * (avgLeftInput/avgLeftOutput));
        output->setSample(1, i, output->getSample(1, i) * (avgRightInput/avgRightOutput));
    }
}

void MainComponent::peakNormalize(AudioBuffer<float>* input, AudioBuffer<float>* output){
    DBG("**** PEAK NORMALIZATION LEFT AND RIGHT");
    double inputMaxLeft = 0.0;
    double inputMaxRight = 0.0;
    double inputMaxLocLeft = 0;
    double inputMaxLocRight = 0;
    double outputMaxLeft = 0.0;
    double outputMaxRight = 0.0;
    double outputMaxLocLeft = 0;
    double outputMaxLocRight = 0;
    for(int i = 0; i < input->getNumSamples(); i++){
        double currentInputLeft = input->getSample(0, i);
        double currentOutputLeft = output->getSample(0, i);
        double currentInputRight = input->getSample(1, i);
        double currentOutputRight = output->getSample(1, i);
    
        if(currentInputLeft > inputMaxLeft){
            inputMaxLeft = input->getSample(0, i);
            inputMaxLocLeft = i;
        }
        if(currentOutputLeft > outputMaxLeft){
            outputMaxLeft = output->getSample(0, i);
            outputMaxLocLeft = i;
        }
        if(currentInputRight > inputMaxRight){
            inputMaxRight = input->getSample(1, i);
            inputMaxLocRight = i;
        }
        if(currentOutputRight > outputMaxRight){
            outputMaxRight = output->getSample(1, i);
            outputMaxLocRight = i;
        }
    }
    
    DBG("**** PEAK NORMALIZATION LEFT AND RIGHT");
    DBG(inputMaxLocLeft);
    DBG(outputMaxLocLeft);
    DBG(inputMaxLocRight);
    DBG(outputMaxLocRight);

    for(int i = 0; i < input->getNumSamples(); i++){
        output->setSample(0, i, output->getSample(0, i) * (inputMaxLeft/outputMaxLeft));
        output->setSample(1, i, output->getSample(1, i) * (inputMaxRight/outputMaxRight));
    }
}

void MainComponent::writeFile(AudioBuffer<float> buffer){
    WavAudioFormat format;
    std::unique_ptr<AudioFormatWriter> writer;
    File OutputFile = File("/Users/davidtemming/ResearchOutputFiles/Additive_" + fileName);
    OutputFile.deleteFile();
    
    DBG("number of samples in the buffer: " + std::to_string(buffer.getNumSamples()));

    writer.reset (format.createWriterFor (new FileOutputStream (OutputFile),
                                          44100,
                                          buffer.getNumChannels(),
                                          24,
                                          {},
                                          0));
    
    if (writer != nullptr){
        writer->writeFromAudioSampleBuffer (buffer, 0, buffer.getNumSamples());
    }
}

void MainComponent::playButtonClicked()
{
    if((state == Stopped) || (state == Paused)){
        transportStateChanged(Starting, &originalTransport, &playButton, &stopButton);
    } else if(state == originalPlaying){
        transportStateChanged(Pausing, &originalTransport, &playButton, &stopButton);
    }
}

void MainComponent::stopButtonClicked()
{
    if(state == Paused){
        transportStateChanged(Stopped, &originalTransport, &playButton, &stopButton);
    } else{
        transportStateChanged(Stopping, &originalTransport, &playButton, &stopButton);
    }
}

void MainComponent::additivePlayButtonClicked(){
//CURRENTLY HAS NO FUNCTION
}

void MainComponent::additiveStopButtonClicked(){
//CURRENTLY HAS NO FUNCTION
}

void MainComponent::transportStateChanged(TransportState newState, AudioTransportSource* transport, TextButton* play, TextButton* stop)
{
    if(newState != state)
    {
        state = newState;
        
        switch (state) {
            case Stopped:
                play->setButtonText("Play");
                stop->setButtonText("Stop");
                stop->setEnabled(false);
                play->setEnabled(true);
                transport->setPosition(0.0);
                break;
                
            case Starting:
                transport->start();
                break;
                
            case originalPlaying:
                play->setButtonText("Pause");
                stop->setButtonText("Stop");
                stop->setEnabled(true);
                break;
            
            case additivePlaying:
                play->setButtonText("Pause");
                stop->setButtonText("Stop");
                stop->setEnabled(true);
                break;
                
            case Pausing:
                transport->stop();
                break;
            
            case Paused:
                play->setButtonText("Resume");
                stop->setButtonText("Return to Zero");
                break;
                
            case Stopping:
                transport->stop();
                break;
        }
    }
}

void MainComponent::changeListenerCallback (ChangeBroadcaster *source)
{
    if (source == &originalTransport)
    {
        if (originalTransport.isPlaying()){
            transportStateChanged(originalPlaying, &originalTransport, &playButton, &stopButton);
        } else if((state == Stopping) || (state == originalPlaying)){
            transportStateChanged(Stopped, &originalTransport, &playButton, &stopButton);
        }else if(Pausing == state){
            transportStateChanged(Paused, &originalTransport, &playButton, &stopButton);
        }
    }
    
    if(source == &inputThumbnail){
        thumbNailChanged();
    }
    if(source == &outputThumbnail){
        thumbNailChanged();
    }

}

void MainComponent::thumbNailChanged(){
    repaint();
}

void MainComponent::getNextAudioBlock (const AudioSourceChannelInfo& bufferToFill)
{
    if(readerSource.get() == nullptr)
    {
        bufferToFill.clearActiveBufferRegion();
        return;
    }
    originalTransport.getNextAudioBlock(bufferToFill);
}

void MainComponent::releaseResources()
{
    originalTransport.releaseResources();
}

//==============================================================================
void MainComponent::paint (Graphics& g)
{
    // (Our component is opaque, so we must completely fill the background with a solid colour)
    g.fillAll (getLookAndFeel().findColour (ResizableWindow::backgroundColourId));
    
    Rectangle<int> inputThumbnailBounds (10, 90, getWidth() / 2 - 20, getHeight() - 210);
    Rectangle<int> outputThumbnailBounds (10 + getWidth() / 2, 90, getWidth() / 2 - 20, getHeight() - 210);
    
    if(inputThumbnail.getNumChannels() == 0){
        paintIfNoFileLoaded (g, inputThumbnailBounds);
    } else{
        paintIfFileLoaded (g, inputThumbnailBounds);
    }
    
    if(resynthesized != true){
        paintIfNotResynthesized(g, outputThumbnailBounds);
    } else{
         paintIfResynthesized(g, outputThumbnailBounds);
    }
}

void MainComponent::paintIfNoFileLoaded (Graphics& g, const Rectangle<int>& thumbnailBounds)
{
    g.setColour (Colours::darkgrey);
    g.fillRect (thumbnailBounds);
    g.setColour (Colours::white);
    g.drawFittedText ("No File Loaded", thumbnailBounds, Justification::centred, 1.0f);
}

void MainComponent::paintIfNotResynthesized(Graphics& g, const Rectangle<int>& thumbnailBounds)
{
    g.setColour (Colours::darkgrey);
    g.fillRect (thumbnailBounds);
    g.setColour (Colours::white);
    g.drawFittedText ("No Resynthesis to Display", thumbnailBounds, Justification::centred, 1.0f);
}

void MainComponent::paintIfFileLoaded (Graphics& g, const Rectangle<int>& thumbnailBounds)
{
    g.setColour (Colours::white);
    g.fillRect (thumbnailBounds);
    
    g.setColour (Colours::red);                                     // [8]
    
    inputThumbnail.drawChannels (g,                                      // [9]
                            thumbnailBounds,
                            0.0,                                    // start time
                            inputThumbnail.getTotalLength(),             // end time
                            1.0f);                                  // vertical zoom
}

void MainComponent::paintIfResynthesized(Graphics& g, const Rectangle<int>& thumbnailBounds)
{
    outputThumbnail.setSource(new FileInputSource(File("/Users/davidtemming/ResearchOutputFiles/Additive_" + fileName)));


    g.setColour (Colours::white);
    g.fillRect (thumbnailBounds);
    
    g.setColour (Colours::blue);                                     // [8]
    
    outputThumbnail.drawChannels (g,                                      // [9]
                                 thumbnailBounds,
                                 0.0,                                    // start time
                                 outputThumbnail.getTotalLength(),             // end time
                                 1.0f);                                  // vertical zoom
}

void MainComponent::resized()
{
    openButton.setBounds(20, 10, getWidth() / 2 - 40, 30);
    processButton.setBounds(20, 50, getWidth() / 2 - 40, 30);
    playButton.setBounds(20, 220, getWidth() / 2 - 40, 30);
    stopButton.setBounds(20, 260, getWidth() / 2 - 40, 30);
    additivePlayButton.setBounds(20 + getWidth() / 2, 220, getWidth() / 2 - 40, 30);
    additiveStopButton.setBounds(20 + getWidth() / 2, 260, getWidth() / 2 - 40, 30);
    fileNameDisplay.setBounds(10, getHeight() - 30, getWidth() - 20, 30);
}

