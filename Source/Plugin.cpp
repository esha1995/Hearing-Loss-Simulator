#include <math.h>
#include <stdio.h>
#include <string>
#include <vector>
#include "FMODAPI.h"
#include <JuceHeader.h>
#include "DSPFilters.h"

extern "C"
{
    F_EXPORT FMOD_DSP_DESCRIPTION* F_CALL FMODGetDSPDescription();
}

// ==================== //
// CALLBACK DEFINITIONS //
// ==================== //
FMOD_RESULT Create_Callback                     (FMOD_DSP_STATE *dsp_state);
FMOD_RESULT Release_Callback                    (FMOD_DSP_STATE *dsp_state);
FMOD_RESULT Reset_Callback                      (FMOD_DSP_STATE *dsp_state);
FMOD_RESULT Read_Callback                       (FMOD_DSP_STATE *dsp_state, float *inbuffer, float *outbuffer, unsigned int length, int inchannels, int *outchannels);
FMOD_RESULT Process_Callback                    (FMOD_DSP_STATE *dsp_state, unsigned int length, const FMOD_DSP_BUFFER_ARRAY *inbufferarray, FMOD_DSP_BUFFER_ARRAY *outbufferarray, FMOD_BOOL inputsidle, FMOD_DSP_PROCESS_OPERATION op);
FMOD_RESULT SetPosition_Callback                (FMOD_DSP_STATE *dsp_state, unsigned int pos);
FMOD_RESULT ShouldIProcess_Callback             (FMOD_DSP_STATE *dsp_state, FMOD_BOOL inputsidle, unsigned int length, FMOD_CHANNELMASK inmask, int inchannels, FMOD_SPEAKERMODE speakermode);

FMOD_RESULT SetFloat_Callback                   (FMOD_DSP_STATE *dsp_state, int index, float value);
FMOD_RESULT SetInt_Callback                     (FMOD_DSP_STATE *dsp_state, int index, int value);
FMOD_RESULT SetBool_Callback                    (FMOD_DSP_STATE *dsp_state, int index, FMOD_BOOL value);
FMOD_RESULT SetData_Callback                    (FMOD_DSP_STATE *dsp_state, int index, void *data, unsigned int length);
FMOD_RESULT GetFloat_Callback                   (FMOD_DSP_STATE *dsp_state, int index, float *value, char *valuestr);
FMOD_RESULT GetInt_Callback                     (FMOD_DSP_STATE *dsp_state, int index, int *value, char *valuestr);
FMOD_RESULT GetBool_Callback                    (FMOD_DSP_STATE *dsp_state, int index, FMOD_BOOL *value, char *valuestr);
FMOD_RESULT GetData_Callback                    (FMOD_DSP_STATE *dsp_state, int index, void **data, unsigned int *length, char *valuestr);

FMOD_RESULT SystemRegister_Callback             (FMOD_DSP_STATE *dsp_state);
FMOD_RESULT SystemDeregister_Callback           (FMOD_DSP_STATE *dsp_state);
FMOD_RESULT SystemMix_Callback                  (FMOD_DSP_STATE *dsp_state, int stage);

// ==================== //
//      PARAMETERS      //
// ==================== //

// set all parameters in ENUM + a NUM_PARAMS which will know how many parameters the program has

enum
{
    HEARINGLOSSL = 0,
    HEARINGLOSSR,
    SMEARLEVELL,
    SMEARLEVELR,
    TEMPDISRUPTL,
    TEMPDISRUPTR,
    ONOFF,
    EXPANSIONL,
    EXPANSIONR,
    MIX,
    NUM_PARAMS
};

// create parameters as FMOD_DSP_PARAMETER DESC
static FMOD_DSP_PARAMETER_DESC hearingLossL, hearingLossR, smearLevelL, smearLevelR, tempDisruptL, tempDisruptR, onOff, expansionL, expansionR, mix;

// create a list with NUM_PARAM elements to get parameters
FMOD_DSP_PARAMETER_DESC* PluginsParameters[NUM_PARAMS] =
{
    &hearingLossL,
    &hearingLossR,
    &smearLevelL,
    &smearLevelR,
    &tempDisruptL,
    &tempDisruptR,
    &onOff,
    &expansionL,
    &expansionR,
    &mix
};


// ==================== //
//     SET CALLBACKS    //
// ==================== //

FMOD_DSP_DESCRIPTION PluginCallbacks =
{
    FMOD_PLUGIN_SDK_VERSION,    // version
    "Emils Hearing Loss Simulator",     // name
    0x00010000,                 // plugin version
    1,                          // no. input buffers
    1,                          // no. output buffers
    Create_Callback,            // create
    Release_Callback,           // release
    Reset_Callback,             // reset
    Read_Callback,              // read
    Process_Callback,           // process
    SetPosition_Callback,       // setposition
    NUM_PARAMS,                             // no. parameter
    PluginsParameters,                          // pointer to parameter descriptions
    SetFloat_Callback,          // Set float
    SetInt_Callback,            // Set int
    SetBool_Callback,           // Set bool
    SetData_Callback,           // Set data
    GetFloat_Callback,          // Get float
    GetInt_Callback,            // Get int
    GetBool_Callback,           // Get bool
    GetData_Callback,           // Get data
    ShouldIProcess_Callback,    // Check states before processing
    0,                          // User data
    SystemRegister_Callback,    // System register
    SystemDeregister_Callback,  // System deregister
    SystemMix_Callback          // Mixer thread exucute / after execute
};

char const* HEARINGLOSS_NAMES[3] = {"OFF", "MILD", "MODERATE"};
char const* SMEAR_NAMES[3] = {"OFF", "LOW", "HIGH"};
char const* TEMP_DISRUPT_NAMES[2] = {"OFF", "ON"};

// assign the parameters with a name and value, to be controlled from FMOD
extern "C"
{
    F_EXPORT FMOD_DSP_DESCRIPTION* F_CALL FMODGetDSPDescription ()
    {
        FMOD_DSP_INIT_PARAMDESC_FLOAT(hearingLossL, "LOSS L", "", "", 0.0f, 1.0f, 0.0f);
        FMOD_DSP_INIT_PARAMDESC_FLOAT(hearingLossR, "LOSS R", "", "", 0.0f, 1.0f, 0.0f);
        FMOD_DSP_INIT_PARAMDESC_FLOAT(smearLevelL, "SMEAR L", "", "", 0.0f, 1.0f, 0.0f);
        FMOD_DSP_INIT_PARAMDESC_FLOAT(smearLevelR, "SMEAR R", "", "", 0.0f, 1.0f, 0.0f);
        FMOD_DSP_INIT_PARAMDESC_FLOAT(tempDisruptL, "TEMP L", "", "", 0.0f, 1.0f, 0.0f);
        FMOD_DSP_INIT_PARAMDESC_FLOAT(tempDisruptR, "TEMP R", "", "", 0.0f, 1.0f, 0.0f);
        FMOD_DSP_INIT_PARAMDESC_BOOL(onOff, "ON", "", "ON OR OFF", true, TEMP_DISRUPT_NAMES);
        FMOD_DSP_INIT_PARAMDESC_FLOAT(expansionL, "RLG L", "", "", 0.0f, 1.0f, 0.0f);
        FMOD_DSP_INIT_PARAMDESC_FLOAT(expansionR, "RLG R", "", "", 0.0f, 1.0f, 0.0f);
        FMOD_DSP_INIT_PARAMDESC_FLOAT(mix, "MIX", "", "SIGNAL MIX", 0.0f, 1.0f, 1.0f);
        return &PluginCallbacks;
    }

}

// ==================== //
//     PLUGIN CLASS     //
// ==================== //

// in the plugin class you do all of your DSP stuff

enum class HearingLoss
{
    OFF = 0,
    MILD,
    MODERATE,
    SEVERE
};

class Plugin
{
public:
    Plugin() {};
    void Prepare (FMOD_DSP_STATE *dsp_state, int numChannels, int numSamples);
    void Release (FMOD_DSP_STATE *dsp_state);
    void Process (float* inbuffer, float* outbuffer, unsigned int length, int channels, FMOD_DSP_STATE *dsp_state);
    void Init(FMOD_DSP_STATE *dsp_state);
    void setParameterFloat(int index, float value);
    void getParameterFloat(int index, float* value);
    void setParamterInt(int index, int value);
    void setParameterBool(int index, bool value);
    
private:
    void SetJuceBufferWithData(juce::AudioBuffer<float>* buffer, float* inbuffer, int numChannels, int numSamples);
    void FillOutputBuffer(juce::AudioBuffer<float>* buffer, juce::AudioBuffer<float>* initBuffer, float *outbuffer,  int numChannels, int numSamples);
    void SumBuffersTogether(juce::AudioBuffer<float>* buffer1, juce::AudioBuffer<float>* buffer2, juce::AudioBuffer<float>* outputBuffer);
    void ApplySpectralSmearing(juce::AudioBuffer<float> *highBuffer, bool left);
    void PrepareAudiogram();
    void PrepareExpansion();
    void ProcessAudiogram(juce::AudioBuffer<float>* leftBuffer, juce::AudioBuffer<float>* rightBuffer);
    void UpdateHearingLoss(bool left, HearingLoss hearingLoss);
    void UpdateSmearingLevel(bool left, int value);
    void ApplyTemporalDisruption(juce::AudioBuffer<float>* buffer, bool left);
    void ApplyExpansion(juce::AudioBuffer<float>* buffer, bool left);
    
    bool hearingLossOn = true;
    bool expansionL = true;
    bool expansionR = true;
    float hearingLossMix = 1.0f;
    float tempDisruptMixL = 0.2f;
    float tempDisruptMixR = 0.2f;
    int sampleRate = 44100;
    int numChannels = 0;
    int numSamples = 0;
    bool prepared = false;
    bool smearL = false;
    bool smearR = false;
    bool tempDisruptL = false;
    bool tempDisruptR = false;
    
    DSP::Butterworth<DSP::LOWPASS,8>* lowpassL;
    DSP::Butterworth<DSP::HIGHPASS,8>* highpassL;
    DSP::Butterworth<DSP::LOWPASS,8>* lowpassR;
    DSP::Butterworth<DSP::HIGHPASS,8>* highpassR;
    DSP::Butterworth<DSP::LOWPASS,4>* smearLowpassL;
    DSP::Butterworth<DSP::LOWPASS,4>* smearLowpassR;
    DSP::Butterworth<DSP::LOWPASS,8>* tempLowpassL;
    DSP::Butterworth<DSP::LOWPASS,8>* tempLowpassR;
    
    HearingLoss hearingLossL = HearingLoss::OFF;
    HearingLoss hearingLossR = HearingLoss::OFF;
    
    std::vector<dsp::ProcessorDuplicator <dsp::IIR::Filter<float>, dsp::IIR::Coefficients <float>>*> audiogramFiltersL;
    std::vector<dsp::ProcessorDuplicator <dsp::IIR::Filter<float>, dsp::IIR::Coefficients <float>>*> audiogramFiltersR;
    
    juce::dsp::FFT* fftL;
    juce::dsp::FFT* fftR;
    juce::dsp::FFT* inverseFftL;
    juce::dsp::FFT* inverseFftR;
    
    juce::dsp::Compressor<float>* expL;
    juce::dsp::Compressor<float>* expR;
    
    const int fftOrder = 10;
};


void Plugin::Init(FMOD_DSP_STATE *dsp_state) {
    FMOD_DSP_GETSAMPLERATE(dsp_state, &sampleRate);
    lowpassL = new DSP::Butterworth<DSP::LOWPASS,8>();
    highpassL = new DSP::Butterworth<DSP::HIGHPASS,8>();
    lowpassR = new DSP::Butterworth<DSP::LOWPASS,8>();
    highpassR = new DSP::Butterworth<DSP::HIGHPASS,8>();
    smearLowpassL = new DSP::Butterworth<DSP::LOWPASS,4>();
    smearLowpassR = new DSP::Butterworth<DSP::LOWPASS,4>();
    tempLowpassL = new DSP::Butterworth<DSP::LOWPASS,8>();
    tempLowpassR = new DSP::Butterworth<DSP::LOWPASS,8>();
    
    lowpassL->setup(static_cast<double>(sampleRate), static_cast<double>(1750));
    highpassL->setup(static_cast<double>(sampleRate), static_cast<double>(1750));
    lowpassR->setup(static_cast<double>(sampleRate), static_cast<double>(1750));
    highpassR->setup(static_cast<double>(sampleRate), static_cast<double>(1750));
    tempLowpassL->setup(static_cast<double>(sampleRate), static_cast<double>(1750));
    tempLowpassR->setup(static_cast<double>(sampleRate), static_cast<double>(1750));
    smearLowpassL->setup(static_cast<double>(sampleRate), static_cast<double>(100));
    smearLowpassR->setup(static_cast<double>(sampleRate), static_cast<double>(100));
    
    int fftOrder = 9;
    if(numSamples == 512)
        fftOrder = 9;
    if(numSamples == 1024)
        fftOrder = 10;
    if(numSamples == 2048)
        fftOrder = 11;
    
    fftL = new juce::dsp::FFT(fftOrder);
    fftR = new juce::dsp::FFT(fftOrder);
    inverseFftL = new juce::dsp::FFT(fftOrder);
    inverseFftR = new juce::dsp::FFT(fftOrder);
}

void Plugin::Prepare(FMOD_DSP_STATE *dsp_state, int numChannels, int numSamples)
{
    if(prepared)
        return;
    
    this->numChannels = numChannels;
    this->numSamples = numSamples;

    PrepareAudiogram();
    PrepareExpansion();
    
    prepared = true;
}

void Plugin::Process(float *inbuffer, float *outbuffer, unsigned int numSamples, int numChannels, FMOD_DSP_STATE *dsp_state)
{
    if(!hearingLossOn)
    {
        *outbuffer = *inbuffer;
        return;
    }
    
    // MAKE JUCE AUDIO BUFFER
    juce::AudioBuffer<float> audioBuffer(2, numSamples);
    SetJuceBufferWithData(&audioBuffer, inbuffer, 2, numSamples);
    juce::AudioBuffer<float> initBuffer;
    initBuffer.makeCopyOf(audioBuffer);
    
    // SPLIT TO LEFT AND RIGHT
    juce::AudioBuffer<float> leftBuffer, rightBuffer;
    leftBuffer.makeCopyOf(audioBuffer); leftBuffer.clear(1, 0, numSamples);
    rightBuffer.makeCopyOf(audioBuffer); rightBuffer.clear(0, 0, numSamples);
    
    //HIGH FREQUENCY ATTENUATION
    ProcessAudiogram(&leftBuffer, &rightBuffer);
    
    // SPLIT SIGNAL TO LOW AND HIGH
    juce::AudioBuffer<float> lowBufferL, highBufferL, lowBufferR, highBufferR;
    // SPLIT LEFT SIGNAL ->
    lowBufferL.makeCopyOf(leftBuffer); highBufferL.makeCopyOf(leftBuffer);
    lowpassL->process(&lowBufferL); highpassL->process(&highBufferL);
    // SPLIT RIGHT SIGNAL ->
    lowBufferR.makeCopyOf(rightBuffer); highBufferR.makeCopyOf(rightBuffer);
    lowpassR->process(&lowBufferR); highpassR->process(&highBufferR);

    // PERFORM SPECTRAL SMEARING ON HIGH SIGNAL
    if(smearL) ApplySpectralSmearing(&highBufferL, true);
    if(smearR) ApplySpectralSmearing(&highBufferR, false);
    
    // PERFORM RLG / EXPANSION ON HIGH SIGNAL
    if(expansionL) ApplyExpansion(&highBufferL, true);
    if(expansionR) ApplyExpansion(&highBufferR, false);

    // PERFORM TEMPORAL DISRUPTION ON LOW SIGNAL
    if(tempDisruptL) ApplyTemporalDisruption(&lowBufferL, true);
    if(tempDisruptR) ApplyTemporalDisruption(&lowBufferR, false);

     
    // SUM EVERYTHING TOGETHER
    juce::AudioBuffer<float> outBufferL(numChannels, numSamples);
    juce::AudioBuffer<float> outBufferR(numChannels, numSamples);
    juce::AudioBuffer<float> outBuffer(numChannels, numSamples);
    SumBuffersTogether(&lowBufferL, &highBufferL, &outBufferL);
    SumBuffersTogether(&lowBufferR, &highBufferR, &outBufferR);
    SumBuffersTogether(&outBufferL, &outBufferR, &outBuffer);
    
    //FILL OUTBUFFER
    FillOutputBuffer(&outBuffer, &initBuffer, outbuffer, numChannels, numSamples);
}

void Plugin::Release(FMOD_DSP_STATE *dsp_state)
{
}

// gets called when a float parameter is changed ->
void Plugin::setParameterFloat(int index, float value) {
    float divide = 1.0f/3.0f;
    float divideLoss = 1.0f/4.0f;
    HearingLoss hearingLossL = HearingLoss::OFF;
    HearingLoss hearingLossR  = HearingLoss::OFF;
    int smearLevelL, smearLevelR = 0;
    
    switch (index) {
        case MIX:
            hearingLossMix = value;
            break;
        case HEARINGLOSSL:
            if(value < divideLoss) hearingLossL = HearingLoss::OFF;
            if(value > divideLoss && value < divideLoss * 2) hearingLossL = HearingLoss::MILD;
            if(value > divideLoss * 2 && value < divideLoss * 3) hearingLossL = HearingLoss::MODERATE;
            if(value > divideLoss * 3) hearingLossL = HearingLoss::SEVERE;
            UpdateHearingLoss(true, hearingLossL);
            break;
        case HEARINGLOSSR:
            if(value < divideLoss) hearingLossR = HearingLoss::OFF;
            if(value > divideLoss && value < divideLoss * 2) hearingLossR = HearingLoss::MILD;
            if(value > divideLoss * 2 && value < divideLoss * 3) hearingLossR = HearingLoss::MODERATE;
            if(value > divideLoss * 3) hearingLossR = HearingLoss::SEVERE;
            UpdateHearingLoss(false, hearingLossR);
            break;
        case SMEARLEVELL:
            if(value < divide) smearLevelL = 0;
            if(value > divide && value < divide * 2) smearLevelL = 1;
            if(value > divide * 2) smearLevelL = 2;
            UpdateSmearingLevel(true, smearLevelL);
            break;
        case SMEARLEVELR:
            if(value < divide) smearLevelR = 0;
            if(value > divide && value < divide * 2) smearLevelR = 1;
            if(value > divide * 2) smearLevelR = 2;
            UpdateSmearingLevel(false, smearLevelR);
            break;
        case TEMPDISRUPTL:
            if(value > 0.0099f) tempDisruptL = true;
            else tempDisruptL = false;
            tempDisruptMixL = value;
            break;
        case TEMPDISRUPTR:
            if(value > 0.0099f) tempDisruptR = true;
            else tempDisruptR = false;
            tempDisruptMixR = value;
            break;
        case EXPANSIONL:
            if(value > 0.0099f) expansionL = true;
            else expansionL = false;
            break;
        case EXPANSIONR:
            if(value > 0.0099f) expansionR = true;
            else expansionR = false;
            break;
        default:
            break;
    }
}

void Plugin::SetJuceBufferWithData(juce::AudioBuffer<float>* buffer, float *inbuffer, int numChannels, int numSamples) {

    for (int s = 0; s < numSamples; s++)
    {
        for (int ch = 0; ch < numChannels; ch++)
        {
            buffer->setSample(ch, s, *inbuffer++);
        }
    }
}

void Plugin::FillOutputBuffer(juce::AudioBuffer<float>* buffer, juce::AudioBuffer<float>* initBuffer, float *outbuffer,  int numChannels, int numSamples) {
    for (int s = 0; s < numSamples; s++)
    {
        for (int ch = 0; ch < numChannels; ch++)
        {
            *outbuffer++ = (hearingLossMix * buffer->getSample(ch, s)) + ((1.0f - hearingLossMix) * initBuffer->getSample(ch, s));
        }
    }
}

void Plugin::SumBuffersTogether(juce::AudioBuffer<float> *buffer1, juce::AudioBuffer<float> *buffer2, juce::AudioBuffer<float> *outputBuffer) { 
    for(int ch = 0; ch < numChannels; ++ch)
    {
        for(int s = 0; s < numSamples; ++s)
        {
            float sample = buffer1->getSample(ch, s) + buffer2->getSample(ch, s);
            outputBuffer->setSample(ch, s, sample);
        }
    }
}

void Plugin::ApplySpectralSmearing(juce::AudioBuffer<float> *highBuffer, bool left)
{
    // make lowpassed noise buffer:
    Random random;
    juce::AudioBuffer<float> noiseBuffer(1, numSamples);
    for(int s = 0; s < numSamples; ++s)
    {
        float noise = (random.nextFloat() * 2.0f - 1.0f) * 0.5f;
        noiseBuffer.setSample(0, s, noise);
    }
    if(left) smearLowpassL->process(&noiseBuffer);
    else smearLowpassR->process(&noiseBuffer);

    // multiply the channel data with the noise buffer
    float* channelData = highBuffer->getWritePointer(left ? 0 : 1);
    for(int s = 0; s < numSamples; ++s)
    {
        // the spectral smearing output ->
        float output = 20 * channelData[s] * noiseBuffer.getSample(0, s);
        
        // mix with original signal ->
        channelData[s] = 0.5f * output + 0.5f * channelData[s];
    }
}

void Plugin::PrepareAudiogram() {
    //SETUP ATTENUATION FILTERS ->
    audiogramFiltersL.resize(6);
    audiogramFiltersR.resize(6);
    
    dsp::ProcessSpec spec;
    spec.sampleRate = sampleRate;
    spec.maximumBlockSize = numSamples;
    spec.numChannels = numChannels;
    
    float frequencies[6] = {250.0f,500.0f,1000.0f,2000.0f,4000.0f,8000.0f};
    float qValues[6] = {1.0f,1.0f,1.0f,1.5f,2.0f,2.5f};
    float offGain[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};

    for(int i = 0; i < audiogramFiltersL.size(); ++i)
    {
        audiogramFiltersL[i] = new juce::dsp::ProcessorDuplicator <dsp::IIR::Filter<float>, dsp::IIR::Coefficients <float>>(dsp::IIR::Coefficients<float>::makePeakFilter(static_cast<double>(sampleRate), frequencies[i], qValues[i], static_cast<float>(Decibels::decibelsToGain(offGain[i]))));
        
        audiogramFiltersR[i] = new juce::dsp::ProcessorDuplicator <dsp::IIR::Filter<float>, dsp::IIR::Coefficients <float>>(dsp::IIR::Coefficients<float>::makePeakFilter(static_cast<double>(sampleRate), frequencies[i], qValues[i], static_cast<float>(Decibels::decibelsToGain(offGain[i]))));
        
        audiogramFiltersL[i]->prepare(spec);
        audiogramFiltersL[i]->reset();
        
        audiogramFiltersR[i]->prepare(spec);
        audiogramFiltersR[i]->reset();
    }
}

void Plugin::ProcessAudiogram(juce::AudioBuffer<float> *leftBuffer, juce::AudioBuffer<float> *rightBuffer) {
    
    juce::dsp::AudioBlock<float> leftBlock(*leftBuffer);
    juce::dsp::AudioBlock<float> rightBlock(*rightBuffer);
    juce::dsp::ProcessContextReplacing<float> leftContext(leftBlock);
    juce::dsp::ProcessContextReplacing<float> rightContext(rightBlock);
    
    for(int i = 0; i < audiogramFiltersL.size(); ++i)
    {
        if(hearingLossL != HearingLoss::OFF) audiogramFiltersL[i]->process(leftContext);
        if(hearingLossR != HearingLoss::OFF) audiogramFiltersR[i]->process(rightContext);
    }
}

void Plugin::setParamterInt(int index, int value) { 
    switch (index) {

        default:
            break;
    }
}

void Plugin::UpdateHearingLoss(bool left, HearingLoss hearingLoss) {
    
    float frequencies[6] = {250.0f,500.0f,1000.0f,2000.0f,4000.0f,8000.0f};
    float qValues[6] = {1.0f,1.0f,1.0f,1.5f,2.0f,2.5f};
    float mildGain[6] = {0.0f, -5.0f, -15.0f, -15.0f, -20.0f, -25.0f};
    float modGain[6] = {-5.0f, -10.0f, -20.0f, -20.0f, -30.0f, -45.0f};
    float sevGain[6] = {-10.0f, -15.0f, -25.0f, -25.0f, -40.0f, -60.0f};
    float offGain[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    
    float* useGain;
    
    switch (hearingLoss)
    {
        case HearingLoss::OFF:
            useGain = offGain;
            break;
        case HearingLoss::MILD:
            useGain = mildGain;
            break;
        case HearingLoss::MODERATE:
            useGain = modGain;
            break;
        case HearingLoss::SEVERE:
            useGain = sevGain;
            break;
    }
    
    if(left)
    {
        hearingLossL = hearingLoss;
        for(int i = 0; i < audiogramFiltersL.size(); ++i)
        {
            *audiogramFiltersL[i]->state = *dsp::IIR::Coefficients<float>::makePeakFilter(static_cast<double>(sampleRate), frequencies[i], qValues[i], static_cast<float>(Decibels::decibelsToGain(useGain[i])));
        }
    }
    else
    {
        hearingLossR = hearingLoss;
        for(int i = 0; i < audiogramFiltersR.size(); ++i)
        {
            *audiogramFiltersR[i]->state = *dsp::IIR::Coefficients<float>::makePeakFilter(static_cast<double>(sampleRate), frequencies[i], qValues[i], static_cast<float>(Decibels::decibelsToGain(useGain[i])));
        }
    }
     
}

void Plugin::UpdateSmearingLevel(bool left, int value) {
    double freqs[3] = {0, 100, 200};
    if(left)
    {
        if(value == 0)
        {
            smearL = false;
            return;
        }
        
        smearL = true;
        smearLowpassL->setup(static_cast<double>(sampleRate), freqs[value]);
    }
    else
    {
        if(value == 0)
        {
            smearR = false;
            return;
        }
        
        smearR = true;
        smearLowpassR->setup(static_cast<double>(sampleRate), freqs[value]);
    }
}

void Plugin::ApplyTemporalDisruption(juce::AudioBuffer<float> *buffer, bool left) {

    const float* inputBuffer = buffer->getReadPointer(left ? 0 : 1);
    float* outputBuffer = buffer->getWritePointer(left ? 0 : 1);
    
    std::vector<std::complex<float>> fftInput(numSamples);
    std::vector<std::complex<float>> fftBuffer(numSamples);

    for(int s = 0; s < numSamples; ++s)
    {
        fftInput[s] = std::complex<float>(inputBuffer[s], 0.0f);
        fftBuffer[s] = std::complex<float>(0.0f, 0.0f);
    }
    
    if(left) fftL->perform(fftInput.data(), fftBuffer.data(), false);
    else fftR->perform(fftInput.data(), fftBuffer.data(), false);
    
    float pi = juce::MathConstants<float>::pi;
    Random random;
    
    for(int bin = 0; bin < numSamples; bin++)
    {
        float magnitude = std::abs(fftBuffer[bin]);
        float phase = std::arg(fftBuffer[bin]);
        float randTheta = (-1*pi / 2) + ((pi / 2) * random.nextFloat()); // random theta
        float mixTheta = phase-(randTheta * (left ? tempDisruptMixL : tempDisruptMixR)); // SHIFT RANDOMLY DEPENDING ON THE TEMP DISRUPT MIX
        fftBuffer[bin] = std::polar(magnitude, mixTheta);
    }
    
    if(left) inverseFftL->perform(fftBuffer.data(), fftInput.data(), true);
    else inverseFftR->perform(fftBuffer.data(), fftInput.data(), true);

    // mix the signals to make it more subtle
    for(int s = 0; s < numSamples; ++s)
    {
        outputBuffer[s] = fftInput[s].real();
    }
         
    
    // remove higher frequencies
   if(left) tempLowpassL->process(buffer);
   else tempLowpassR->process(buffer);
}

void Plugin::setParameterBool(int index, bool value) { 
    switch (index) {
        case ONOFF:
            hearingLossOn = value;
            break;
        default:
            break;
    }
}

void Plugin::PrepareExpansion() { 
    dsp::ProcessSpec spec;
    spec.sampleRate = sampleRate;
    spec.maximumBlockSize = numSamples;
    spec.numChannels = numChannels;
    
    expL = new juce::dsp::Compressor<float>();
    expR = new juce::dsp::Compressor<float>();
    
    expL->prepare(spec);
    expR->prepare(spec);
    
    expL->setRatio(0.6f);
    expL->setAttack(0.3f);
    expL->setRelease(0.2f);
    expL->setThreshold(-50.0f);
    
    expR->setRatio(0.6f);
    expR->setAttack(0.3f);
    expR->setRelease(0.2f);
    expR->setThreshold(-50.0f);
}

void Plugin::ApplyExpansion(juce::AudioBuffer<float>* buffer, bool left)
{
    juce::dsp::AudioBlock<float> block(*buffer);
    juce::dsp::ProcessContextReplacing<float> context(block);
    if(left) expL->process(context);
    else expR->process(context);
}

















// ======================= //
// CALLBACK IMPLEMENTATION //
// ======================= //

FMOD_RESULT Create_Callback                     (FMOD_DSP_STATE *dsp_state)
{
    // create our plugin class and attach to fmod
    Plugin* state = (Plugin* )FMOD_DSP_ALLOC(dsp_state, sizeof(Plugin));
    dsp_state->plugindata = state;
    if (!dsp_state->plugindata)
    {
        return FMOD_ERR_MEMORY;
    }
    state->Init(dsp_state);
    return FMOD_OK;
}

FMOD_RESULT Release_Callback                    (FMOD_DSP_STATE *dsp_state)
{
    // release our plugin class
    Plugin* state = (Plugin* )dsp_state->plugindata;
    state->Release(dsp_state);
    FMOD_DSP_FREE(dsp_state, state);
    return FMOD_OK;
}

FMOD_RESULT Reset_Callback                      (FMOD_DSP_STATE *dsp_state)
{
    return FMOD_OK;
}

FMOD_RESULT Read_Callback                       (FMOD_DSP_STATE *dsp_state, float *inbuffer, float *outbuffer, unsigned int length, int inchannels, int *outchannels)
{
    return FMOD_OK;
}

FMOD_RESULT Process_Callback                    (FMOD_DSP_STATE *dsp_state, unsigned int length, const FMOD_DSP_BUFFER_ARRAY *inbufferarray, FMOD_DSP_BUFFER_ARRAY *outbufferarray, FMOD_BOOL inputsidle, FMOD_DSP_PROCESS_OPERATION op)
{
    Plugin* state = (Plugin* )dsp_state->plugindata;
    
    switch (op) {
        case FMOD_DSP_PROCESS_QUERY:
            if (outbufferarray && inbufferarray)
            {
            outbufferarray[0].bufferchannelmask[0] = inbufferarray[0].bufferchannelmask[0];
            outbufferarray[0].buffernumchannels[0] = inbufferarray[0].buffernumchannels[0];
            outbufferarray[0].speakermode       = inbufferarray[0].speakermode;
            }
            
            
            
            if (inputsidle)
            {
                return FMOD_ERR_DSP_DONTPROCESS;
            }
            
            // QUERY is before process is run first time, here we call prepare function of plugin to give number of channels, number of samples, sample rate
            state->Prepare(dsp_state, outbufferarray[0].buffernumchannels[0], length);
            
            break;
            
        case FMOD_DSP_PROCESS_PERFORM:
            
            if (inputsidle)
            {
                return FMOD_ERR_DSP_DONTPROCESS;
            }
            
            // actually process
            state->Process(inbufferarray[0].buffers[0], outbufferarray[0].buffers[0], length, outbufferarray[0].buffernumchannels[0], dsp_state);
            
            return FMOD_OK;
            break;
    }
    
    return FMOD_OK;
}

FMOD_RESULT SetPosition_Callback                (FMOD_DSP_STATE *dsp_state, unsigned int pos)
{
    return FMOD_OK;
}

FMOD_RESULT ShouldIProcess_Callback             (FMOD_DSP_STATE *dsp_state, FMOD_BOOL inputsidle, unsigned int length, FMOD_CHANNELMASK inmask, int inchannels, FMOD_SPEAKERMODE speakermode)
{
    if (inputsidle)
    {
        return FMOD_ERR_DSP_DONTPROCESS;
    }
    return FMOD_OK;
}

FMOD_RESULT SetFloat_Callback                   (FMOD_DSP_STATE *dsp_state, int index, float value)
{
    Plugin* state = (Plugin* )dsp_state->plugindata;
    state->setParameterFloat(index, value);
    return FMOD_OK;
}

FMOD_RESULT SetInt_Callback                     (FMOD_DSP_STATE *dsp_state, int index, int value)
{
    Plugin* state = (Plugin* )dsp_state->plugindata;
    state->setParamterInt(index, value);
    return FMOD_OK;
}

FMOD_RESULT SetBool_Callback                    (FMOD_DSP_STATE *dsp_state, int index, FMOD_BOOL value)
{
    Plugin* state = (Plugin * )dsp_state->plugindata;
    state->setParameterBool(index, (bool)value);
    return FMOD_OK;
}

FMOD_RESULT SetData_Callback                    (FMOD_DSP_STATE *dsp_state, int index, void *data, unsigned int length)
{
    return FMOD_OK;
}

FMOD_RESULT GetFloat_Callback                   (FMOD_DSP_STATE *dsp_state, int index, float *value, char *valuestr)
{
    Plugin* state = (Plugin* )dsp_state->plugindata;
    return FMOD_OK;
}

FMOD_RESULT GetInt_Callback                     (FMOD_DSP_STATE *dsp_state, int index, int *value, char *valuestr)
{
    return FMOD_OK;
}

FMOD_RESULT GetBool_Callback                    (FMOD_DSP_STATE *dsp_state, int index, FMOD_BOOL *value, char *valuestr)
{
    return FMOD_OK;
}

FMOD_RESULT GetData_Callback                    (FMOD_DSP_STATE *dsp_state, int index, void **data, unsigned int *length, char *valuestr)
{
    return FMOD_OK;
}

FMOD_RESULT SystemRegister_Callback             (FMOD_DSP_STATE *dsp_state)
{
    return FMOD_OK;
}

FMOD_RESULT SystemDeregister_Callback           (FMOD_DSP_STATE *dsp_state)
{
    return FMOD_OK;
}

FMOD_RESULT SystemMix_Callback                  (FMOD_DSP_STATE *dsp_state, int stage)
{
    return FMOD_OK;
}
