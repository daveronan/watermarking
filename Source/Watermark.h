/*
  ==============================================================================

    Watermark.h
    Created: 6 May 2016 5:57:12pm
    Author:  David

  ==============================================================================
*/

#ifndef WATERMARK_H_INCLUDED
#define WATERMARK_H_INCLUDED


#include "../JuceLibraryCode/JuceHeader.h"


#include "fletcher.h"
#include "watermark.h"


#define	MAX_CHAR			1024		// max Hex numbers to be embedded in 
//#define	NFREQ				2048
#define	NFREQ				2048
#define	WINDOWSIZE			11.14558	// window size in seconds
//#define	WINDOWSIZE			12	// window size in seconds
#define	NCHMAX				2			// max # of input channels
#define FWMIN				2000		// min freq of a chip subband in Hz
#define	FWMAX				7200		// max freq of a chip subband in Hz
#define	BITSPERWINDOW		4			// bits exored in a window
#define	FRAMESPERWINDOW		24			// window size in frames
#define	OFFSET				4			// chip magnitude change in dB
#define CHIPSPERBLOCK	120			// chips added in a block [60 to 100]
#define SPEEDUP				0.01		// resilience to dynamic freq stretch
#define	INVMAXX				(1.0/((float) LONG_MAX)) 
#define NWATERMARKS			16			// #CCIbits = log2(NWATERMARKS)
#define	NONE				(-1)		// value for no detection
#define MINAMPLITUDE		3.0517578E-5


#define CF				10			// lowpass filter on the cepstrum
#define DECISIONBAR		7.75   		// if (nc>) then watermark detected
#define PREDECISIONBAR	0.50		// if (nc>) then watermark detected
#define PM				3			// cepstrum amplitude clip 
// Pre-echo control parameters

#define NSEC				8			// subblocks per block
#define ETHR				1E-8		// to stabilize energy ratio 
#define	ERLIM				150.		// max energy ratio of two subbblocks

#define	DB(x)			pow(10.0,(x)/20.0)
// Set the state of internal RN generator
// A better set of seeds should be defined later

// Watermark detection parameters

#define DBCUT				-4.0		// how many dBs to cut
#define SRTIME				7			// how many time scales to search 
#define SRFREQ				7			// how many freq scales to search
#define TIME_RESILIENCE		10.			// resilience to time scale in %
#define FREQ_RESILIENCE		5.			// resilience to freq shift in %
#define FDNOISEFLOOR		1e-4		// min subband energy to correlate

// Search steps

#define BASIC_STEP			0.3			// basic search step in frames 
// COMMENT: load_buffer() is written such that integrationarea=2*BASIC_STEP 
// #define SEARCHSTEP			0.13931973	// search step in time
#define	iSTUBBORNESS		WINDOWSIZE/2
#define	sSTUBBORNESS		WINDOWSIZE*3/4
#define	STUBBORNESS			WINDOWSIZE/2
#define xSTUBBORNESS		WINDOWSIZE*3/8
#define sPROGRESSSTEP		WINDOWSIZE*3/4
#define uPROGRESSSTEP		WINDOWSIZE*2/4

#define SPECTRUM_NORMALIZATION

#define	INITIALIZE_SEED(seed, i)	srand((71912317*(i*(i+2)*(seed+1))+(i+4779)*317*(seed))%15991)

	// an audio clip
	typedef struct {
		int				bitsloaded;
		unsigned int	xCCI[MAX_CHAR];
		unsigned int	xLOAD[MAX_CHAR];
	} WMBITS;

	typedef struct {
	   int fstart[CHIPSPERBLOCK];	// subband lower limit for SS bit
	   int fend[CHIPSPERBLOCK];		// subband upper limit for SS bit
	   int fmiddle[CHIPSPERBLOCK];	// pointers to center of subbands
	} SSBANDS;

	// Structures for detection

	typedef struct {				// circular buffer for MCLT data
		float		*buffer;
		short int	*ht;
		long		length;
		long 		pointer;
	} CIRCULAR_BUFFER;

	// Structure for SS bits detection

	typedef struct {
	   int fmiddle[CHIPSPERBLOCK][SRFREQ];	// pointers to center of subbands
	   int cbe[CHIPSPERBLOCK][SRFREQ];
	   int cbs[CHIPSPERBLOCK][SRFREQ];
	} SSDBANDS;

	// Structure that contains the correlation values
	typedef struct {
		int		card[BITSPERWINDOW][2];
		float   corr[BITSPERWINDOW][2];
		float 	sumsquares;
	} PAYLOAD;

	typedef struct {
		float	dtime;
		float 	nc;
		int		cci;
		int		load;
		float 	payload[BITSPERWINDOW];
	} RESULT;

	typedef struct {
		int		cci;
		int		load;

	} WATERMARKFOUNDINFO;


class Watermark
{
public:

	Watermark()
	{
	}
	~Watermark()
	{

	}

	void readWMBits(std::string watremarktoapply);
	void resetBuffer(int blocks2load);
	void addWaterMark(AudioSampleBuffer &inputbuffer);
	void detectWaterMark(AudioSampleBuffer &inputbuffer);
	void initialize(float samplerate, int numchannels, int buffersize);
	static std::string stringToHex(const std::string& word);
	std::string hexToString(const std::string& input) const;
	std::string getDebugInfo(void){return _debugInfo;};
	std::string getUIUpdate();
	void cleanUp();

private:

	//Application of Watermark
	void genSubIndex(SSBANDS *ssbands, long fs, float scalefactor) const;
	float warpFreq(float x, float a) const;
	static void createPermutation (int seed, int (&permuatation)[CHIPSPERBLOCK][BITSPERWINDOW]);
	void getChips(int seed, int cci, int usage);
	static int getWMFlag(float *xx);
	

	std::string _debugInfo;
	float _sampleRate;
	int _numChannels;
	float _bitSize;
	int _bufferSize;
	float _nBlocks;
	float _freqMin;
	float _freqMax;
	float* _nfreqBufferApply;

	float	_blocktlength;		// how many seconds per block
	float	_frametlength;		// how many seconds per frame
	float	_frametime;			// how many seconds within frame
	float	_windowtime;			// how many seconds within the window
	long	_window;				// window counter
	long	_frame;				// frame counter
	bool	_firstrun;

	float	_xx[NCHMAX][NFREQ];	// time-domain samples
	float	_Xc[NCHMAX][NFREQ];		// MCLT coefficients, cosine part
	float	_Xs[NCHMAX][NFREQ];		// MCLT coefficients, sine part
	float	_h[NFREQ];				// MCLT window
	float	_ya[NCHMAX][3*NFREQ/2];	// MCLT internal buffers - analysis
	float	_ys[NCHMAX][3*NFREQ/2];	// MCLT internal buffers - synthesis

	SSBANDS	_ssBands;
	WMBITS _wmbits;
	int		_permutationAdd[CHIPSPERBLOCK][BITSPERWINDOW];
	// chips of the watermark
	int		_watermarkAdd[CHIPSPERBLOCK][FRAMESPERWINDOW][NWATERMARKS];
	float	_chip[NFREQ][FRAMESPERWINDOW];	// watermarking chips, one per MCLT component
								// (used for insertion of CCI btis)

	/////////////////////////////////////////////////////////////////////////////////////////////////////
	//Detect

	void createDecoder (long fs);
	void createPointers ();
	void initResults();
	void loadBuffer(AudioSampleBuffer &inputbuffer, CIRCULAR_BUFFER *cb, long blocks2load, int Nch, long storageblock, long fs);
	static void cepstrumFiltering(float  *inout, float  *cepstrum);
	int computeCorrelationsNormalized (CIRCULAR_BUFFER *cb, long storageblock, int *timeINDEX, int *freqINDEX);
	void createWatermark (int seed, int (&watermark)[CHIPSPERBLOCK][FRAMESPERWINDOW][NWATERMARKS]) const;
	void ComputeMeanCorrelations(int* timeINDEX, int* freqINDEX, int counter, float maximum, float wmean, float wstd);
	void AddSumOfCoeffs(float blockcorr, int blockcard, long i, int frame, int timeindex, int freqindex, int bit, long w);
	void intializeCorrelationSums(int timeindex, int freqindex);
	void performCorrelationTest(int* timeINDEX, int* freqINDEX, int& counter, float temp, float& maximum, long w, int freqindex, int timeindex);
	void pointToDifferentSubBandsCorrelation(CIRCULAR_BUFFER* cb, long k, long p1, float blockcorr, int blockcard, long i, int frame, int timeindex, int freqindex, int bit);
	static int signOfNum(float x);
	float correctAmplitudeToGreaterThan90db(float amplitude) const;
	
	// Buffers for cepstrum filtering and audibility testing
	float       _Ht[NFREQ];				// hearing threshold
	float		_bufA[NFREQ];			// x(w) of the current 3-block
	short int	_htA[NFREQ];				// audibility of the current 3-block

	float* _nfreqBufferDetect;

	float	_xxWat[NCHMAX][NFREQ];	// time-domain samples
	float	_XcWat[NFREQ];		// MCLT coefficients, cosine part
	float	_XsWat[NFREQ];		// MCLT coefficients, sine part
	float	_hWat[NFREQ];				// MCLT window
	float	_yaWat[3*NFREQ/2];	// MCLT internal buffers - analysis
	float	_ybWat[3*NFREQ/2];	// MCLT internal buffers - synthesis

	bool	_bufferFull;
	AudioSampleBuffer _accumBuffer;
	int     _bufferIndex;

	float	_dbcut, _dbnoisefloor, _fdnoisefloor;

	int		_freqINDEX, _timeINDEX;

	float _nDetectBlocks;
	float _blocksize;
	long _blocksperwindow;
	long _blocksperbit;
	long _blocksperframe;
	long _framesize;
	long _searchstep;
	long _framesperbit;
	long _storageblock;
	long _currentwindow;
	long _blocks2load;	// various parameters of the scheme
	float _windowtimeDetect;

	// Watermarks and pointers
	//-------------------------------------------------
	// chips of the watermark
	int			_watermarkDetect[CHIPSPERBLOCK][FRAMESPERWINDOW][NWATERMARKS];
	// pointers to blocks for different timescales
	long		_pointers[FRAMESPERWINDOW][SRTIME];
	// pointers to subbands for SS bits
	SSDBANDS	_ssdecoder;
	// starting and ending frequency magnitude
	int			_startfreq, _endfreq;

	int		_permutationDetect[CHIPSPERBLOCK][BITSPERWINDOW];

	//Circular buffer for detection
	CIRCULAR_BUFFER	_cb;

	// program management variables
	int				_state; 
	float			_stubborness;

	// statistics variables
	float			_temp, _dthres;
	int				_counter, _wdetected, _usagebits; 

	// individual correlations for each payload bit
	PAYLOAD		_payload[SRTIME][SRFREQ][NWATERMARKS];
	// circular buffer of last three correlation tests
	RESULT		_buffer[3];
	// best case correlation test
	RESULT		_bestcase;
	// pointer to end of the circular buffer with results
	int			_pbuff;
	//Keeps track of the watermark letters.
	std::vector<WATERMARKFOUNDINFO> _waterMarkFound;

};

#endif  // WATERMARK_H_INCLUDED
