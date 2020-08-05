/*
  ==============================================================================

    Watermark.cpp
    Created: 6 May 2016 5:57:12pm
    Author:  David

  ==============================================================================
*/

#include "watermark.h"
#include "fxform.h"
#include "hthres.h"
#include "waveio.h"
#include <string>
#include <sstream>

void Watermark::readWMBits(std::string watremarktoapply)
{
	int				i, result = 0;
	unsigned char	temp;

	_wmbits.bitsloaded = static_cast<int>(watremarktoapply.length());

	for (i = 0; i < MAX_CHAR; i++) 
	{
		_wmbits.xCCI[i] =  0;
		_wmbits.xLOAD[i] = 0;
	}

	for (i = 0; i < _wmbits.bitsloaded; i++) 
	{
		temp = watremarktoapply[i];
		if (temp >= 48 && temp <= 57)
		{
			result = temp - 48;
		}
		else if (temp >= 65 && temp <= 70)
		{
			result = temp - 55;
		}
		else if (temp >= 97 && temp <= 102)
		{
			result = temp - 87;
		}

		if (i % 2 == 0)
		{
			_wmbits.xCCI[i>>1] =  result;
		}
		else
		{
			_wmbits.xLOAD[i>>1] = result;
		}
			
	}
		
	// Report values
	printf("WM bits: CCI  = ");
	for (i = 0; i < (_wmbits.bitsloaded + 1) >> 1; i++)
		printf("%X", _wmbits.xCCI[i]);
	printf("\nWM bits: LOAD = ");
	for (i = 0; i < _wmbits.bitsloaded >> 1; i++)
		printf("%X", _wmbits.xLOAD[i]);
	printf("\n");
}

void Watermark::resetBuffer(int blocks2load)
{
	_accumBuffer.setSize(2, NFREQ * blocks2load * _searchstep, 0, 1, 1);
	_nDetectBlocks = blocks2load * _searchstep;

	for(int i = 0; i < _accumBuffer.getNumChannels(); i++)
	{
		float* bufferData = _accumBuffer.getWritePointer(i);

		for(int j =0; j < NFREQ * blocks2load  * _searchstep; j++)
		{
			bufferData[j] = 0;
		}
	}
}

void Watermark::initialize(float samplerate, int numchannels, int buffersize)
{
	_sampleRate = samplerate;
	_numChannels = numchannels;

	_bufferSize = buffersize;
	_firstrun = true;
	
	_freqMin = static_cast<float>((FWMIN*2.0/_sampleRate) * static_cast<float>(NFREQ));
	_freqMax = static_cast<float>((FWMAX*2.0/_sampleRate) * static_cast<float>(NFREQ));

	// Initialize counters
	_blocktlength = NFREQ / static_cast<float>(_sampleRate);
	_frametlength = static_cast<float> (WINDOWSIZE / FRAMESPERWINDOW); 

	_window = 0;
	_frame = 0;
	_frametime = 0.0;
	_windowtime = 0.0;

		//	Initialize MCLT window
	mlt_sine_window(_h, NFREQ);
		
	//// Initialize subband pointers; for insertion, no frequency shift
	genSubIndex(&_ssBands, _sampleRate, 1); 

		// Get first set of watermarking chips
	createPermutation(_window, _permutationAdd);
	getChips(_window, _wmbits.xCCI[_window], _wmbits.xLOAD[_window]);

	for(int i=0; i < NCHMAX; i++)
	{
		for(int j=0; j < NFREQ; j++)
		{
			_xx[i][j] = 0;
			_Xc[i][j] = 0;
			_Xs[i][j] = 0;
		}

		for(int j=0; j < 3 * NFREQ / 2; j++)
		{
			_ya[i][j] = 0;
			_ys[i][j] = 0;
		}
	}
	
	////////////////////////////Detect

		// constants
	_dbcut = pow(10, (DBCUT)/10.);
	_dbnoisefloor = 20.*log10(FDNOISEFLOOR);
	_fdnoisefloor = FDNOISEFLOOR * FDNOISEFLOOR;

	for(int i=0; i < NCHMAX; i++)
	{
		for(int j=0; j < NFREQ; j++)
		{
			_xxWat[i][j] = 0;
			_XcWat[j] = 0;
			_XsWat[j] = 0;
			_Ht[j] = 0.;
		}

		for(int j=0; j < 3 *NFREQ / 2; j++)
		{
			_yaWat[j] = 0;
			_ybWat[j] = 0;
		}
	}

	// blocksize = block duration in seconds
	_blocksize = NFREQ / static_cast<float>(_sampleRate); 	 // length (secs) of each block
	// blocksperwindow = blocks per window
	_blocksperwindow = static_cast<long>(WINDOWSIZE/_blocksize);
	// blocksperbit = blocks per payload bit
	_blocksperbit = static_cast<long>(_blocksperwindow) / BITSPERWINDOW;
	// blocksperframe = blocks per frame
	_blocksperframe = static_cast<long>(_blocksperwindow) / FRAMESPERWINDOW;
	// framesize = frame duration in seconds
	_framesize = _blocksize * _blocksperframe;
	// simple search step
	_searchstep = static_cast<long>((_blocksperframe * BASIC_STEP));
	// number of frames per bit
	_framesperbit = static_cast<long>(FRAMESPERWINDOW/BITSPERWINDOW);
	
	// secret generation
	_currentwindow = 0;
	createPermutation(_currentwindow, _permutationDetect);
	createWatermark(_currentwindow, _watermarkDetect);
	createDecoder(_sampleRate);
	createPointers();
	
	mlt_sine_window(_hWat, NFREQ);

	// circular buffer initialization
	// size of a block in freq magnitudes
	_storageblock = _endfreq - _startfreq + 1;
	// size in blocks of the circular buffer
	_cb.length = static_cast<long>(ceil(_blocksperwindow * (1. + TIME_RESILIENCE / 100.) / _searchstep));
	// pointer to a current begining of the circular buffer
	_cb.pointer = 0;
	// size of memory to allocate
	long mem_i = _cb.length * _storageblock;
	// memory allocation

	_cb.buffer = new float[mem_i];
	_cb.ht = new short int[mem_i];

	// reset data structures
	_blocks2load = _cb.length;
	resetBuffer(_blocks2load);

	_windowtimeDetect = 0.0;

	_bufferIndex = 0;
	_bufferFull = false;

	// init result management strucutres
	initResults();
	// init other variables
	_stubborness = iSTUBBORNESS;
	_wdetected = 0;
	_state = 0;

}

void Watermark::initResults()
{
	_pbuff = 0;
	_bestcase.dtime = -1000.;
	_bestcase.nc = -1000.;		
	_bestcase.cci = NONE;
	_bestcase.load = NONE;
	for (int i =  0; i < BITSPERWINDOW; i++) 
	{
		_bestcase.payload[i] = 0.;
	}
	for (int i = 0; i < 3; i++) 
	{
		_buffer[i].dtime = -1000.;	
		_buffer[i].nc = -1000.;		
		_buffer[i].cci = NONE;	
		_buffer[i].load = NONE;
		for (int j =  0; j < BITSPERWINDOW; j++) 
		{
			_buffer[i].payload[j] = 0;
		}
	}
	for (int i = 0; i < NFREQ; i++)
	{
		_bufA[i] = 0.;
		_htA[i] = 0;
	}
}

void Watermark::loadBuffer(AudioSampleBuffer& inputbuffer, CIRCULAR_BUFFER* cb, long blocks2load, int Nch, long storageblock, long fs)
{
	for (int i = 0; i < inputbuffer.getNumChannels(); i++)
	{
		_accumBuffer.copyFrom(i, _bufferIndex, inputbuffer, i, 0, inputbuffer.getNumSamples());
	}

	_bufferIndex += inputbuffer.getNumSamples();

	if(_bufferIndex + inputbuffer.getNumSamples() > _accumBuffer.getNumSamples())
	{
		_bufferIndex = 0;
		_bufferFull = true;		
	}

	if(_bufferFull)
	{
		// preload data to circular buffer		
		for (int i = cb->pointer, z = 0; i < cb->pointer + blocks2load; i++, z++) 
		{
			// define the location of the circular buffer start
			long index = storageblock * (i % cb->length);
			// accept the previous subsum and reset the future one;
			for (long j = _startfreq; j <= _endfreq; j++) 
			{
				cb->ht[index + j -_startfreq] =  _htA[j];
				_htA[j] = 0;

				cb->buffer[index + j -_startfreq] = _bufA[j];
				_bufA[j] = 0.;
			}

			// load blocks of freq magnitudes
			for (long k = 0; k < _searchstep; k++) 
			{
				for (int ch = 0; ch < Nch; ch++) 
				{	
					if(((z * _searchstep * NFREQ) + (k * NFREQ) + NFREQ) <= _accumBuffer.getNumSamples())
					{
						_nfreqBufferDetect = _accumBuffer.getWritePointer(ch, ((z * _searchstep * NFREQ) + (k * NFREQ)));
					}
					else
					{
						_nfreqBufferDetect = nullptr;
					}

					if (_nfreqBufferDetect != nullptr)
					{				
				
						for (int j = 0; j < NFREQ; j++) 
						{
							_nfreqBufferDetect[j] = correctAmplitudeToGreaterThan90db(_nfreqBufferDetect[j]);
							_xxWat[ch][j] = _nfreqBufferDetect[j];
						}
			
						// Compute direct MCLT, results in vectors Xc and Xs
						fmclt(_XcWat, _XsWat, _xxWat[ch], _yaWat, _hWat, NFREQ);
						// Compute freq magnitudes & block energy            
					
						float	en = 0;
						float	P = 0.1;
						float	a = 0.005;
						float	b = 0.03;
						float	ga, gb, ec, gain;

						for (long j = 0; j < NFREQ; j++) 
						{
							_xxWat[ch][j] = _XcWat[j] * _XcWat[j] + _XsWat[j] * _XsWat[j];
							en += _xxWat[ch][j];
						}

						en /= NFREQ; // en now contains block energy
						gain = 1;

						if (en < P) 
						{
							ga = (P - b) / (P - a);
							gb = P * (b - a) / (P - a);
							if (en > a) 
							{
								ec = ga * en + gb;
							} 
							else 
							{
								ec = (b / a) * en;
							}
							// ec has desired block energy (after dynamic range compression)
							gain = ec/en;
							// printf("Block energy = %g\n", en);
							for (long j = 0; j < NFREQ; j++) 
							{
								_xxWat[ch][j] *= gain;
							}
						}
					
						// Mark inaudible freq magnitudes
						// caution: the power of the signal is used 
						// for faster hearing threshold computation.
						hthres(_Ht, _xxWat[ch], NFREQ, fs); 
						// perform the comparison of audibility and go to dB domain
						for (long j = 0; j < NFREQ; j++) 
						{
							if (_xxWat[ch][j] > _Ht[j] * _dbcut && _xxWat[ch][j] > _fdnoisefloor)
							{
								_Ht[j] = 1.;
							}
							else 
							{
								_Ht[j] = 0.;
							}

							_xxWat[ch][j] = 10.* log10(_xxWat[ch][j]);
						}
						// cepstrum filtering of the signal in dB
						cepstrumFiltering(&_xxWat[ch][0], &_ybWat[0]);
						// copying the data into buffers
						int dave2 =0;
						for (long j = _startfreq; j <= _endfreq; j++) 
						{
							if (_Ht[j] > 0.) 
							{
								cb->buffer[index + j -_startfreq] += static_cast<float>(_xxWat[ch][j]);
								cb->ht[index + j -_startfreq] += static_cast<short int>(_Ht[j]);
								_bufA[j] += static_cast<float>(_xxWat[ch][j]);
								_htA[j] += static_cast<short int>(_Ht[j]);
								dave2++;
							}
						}
						int stop;
					} // both channels processed.
				}
			}
		}
	int dave2=0;
	}
	
}

void Watermark::detectWaterMark(AudioSampleBuffer& inputbuffer)
{
	loadBuffer(inputbuffer, &_cb, _blocks2load, _numChannels, _storageblock, _sampleRate);

	if(_bufferFull)
	{
		_bufferFull = false;
		// audio clip buffer loaded and processed
		// when data is loaded the beginning of window starts here
		_cb.pointer = (_cb.pointer + _blocks2load) % _cb.length;
		printf("Time=%3.3fsec ", _windowtimeDetect);

		_counter = computeCorrelationsNormalized(&_cb, _storageblock, &_timeINDEX, &_freqINDEX);

		// once again counter stores the bits embedded by selecting one out of NWATERMARKS
		// the sign of the partial correlations for the max overall normalized correlation
		// indicates the payload embedded by XORing the watermark.
		// now lets see what are the bits stored in the payload.
		_usagebits = 0;
		for (int bit = 0; bit < BITSPERWINDOW; bit++) 
		{
			_temp = _payload[_timeINDEX][_freqINDEX][_counter].corr[bit][1]/_payload[_timeINDEX][_freqINDEX][_counter].card[bit][1];
			_temp -= _payload[_timeINDEX][_freqINDEX][_counter].corr[bit][0]/_payload[_timeINDEX][_freqINDEX][_counter].card[bit][0];
			
			if (_temp < 0)
			{
				_usagebits = _usagebits | (0x1 << bit);
			}
		}

		// now lets store the results to the result buffer
		// the result buffer stores the results of the last three correlation tests
		// watermark is claimed tobe detected if on all of these tests 
		// normalized correlations consistently indicated watermark existence
		_pbuff = (_pbuff + 1) % 3;
		_buffer[_pbuff].dtime = _windowtimeDetect;
		_buffer[_pbuff].nc = _payload[_timeINDEX][_freqINDEX][_counter].sumsquares;
		_buffer[_pbuff].cci = _counter;
		_buffer[_pbuff].load = _usagebits & 0xF;

		for (long bit = 0; bit < BITSPERWINDOW; bit++) 
		{
			_temp = _payload[_timeINDEX][_freqINDEX][_counter].corr[bit][1]/_payload[_timeINDEX][_freqINDEX][_counter].card[bit][1];
			_temp -= _payload[_timeINDEX][_freqINDEX][_counter].corr[bit][0]/_payload[_timeINDEX][_freqINDEX][_counter].card[bit][0];
			_buffer[_pbuff].payload[bit] =_temp;
		}

		// print results

		char debug [100];
		sprintf(debug, "Time=%3.3fsec  [NC=%8.4f WM=%X%X]\n",_windowtimeDetect, _payload[_timeINDEX][_freqINDEX][_counter].sumsquares, _counter, _usagebits);
		std::string strdebug = debug;

		// check for detection
		_dthres = static_cast<float>(DECISIONBAR + 1e-5*(static_cast<float>(rand() & 8191) - 4095.5)); // Dinei's trick
		_temp = _buffer[0].nc + _buffer[1].nc + _buffer[2].nc;
		if (_temp / 3 > _dthres) 
		{
			_wdetected = 1;
		}
		else 
		{
			_wdetected = 0;
		}

		// watermark has been found before. 
		if (_state == 0) 
		{
			// watermark found for the first time
			if (_wdetected == 1 && _buffer[0].cci == _buffer[1].cci && _buffer[1].cci == _buffer[2].cci) 
			{
				_state = 1;
				_bestcase.dtime = _windowtimeDetect;
				_bestcase.nc = _buffer[0].nc + _buffer[1].nc + _buffer[2].nc;
				_bestcase.cci = _buffer[(_pbuff + 2) % 3].cci;
				_bestcase.load = _buffer[(_pbuff + 2) % 3].load;
				for (long bit = 0; bit < BITSPERWINDOW; bit++) 
				{
					_bestcase.payload[bit] = _buffer[0].payload[bit] + _buffer[1].payload[bit] + _buffer[2].payload[bit];
				}
				_blocks2load = 1;
				resetBuffer(_blocks2load);
				_stubborness = -1e+20;
			} 
			else if (_stubborness < 0) 
			{
				// watermark not found in current window
				// change window
				_state = 0; // state is SEARECH
				if (_currentwindow == 0) 
				{
					_blocks2load = 1;
					resetBuffer(_blocks2load);
					_stubborness = sSTUBBORNESS;
				} 
				else 
				{
					_blocks2load = static_cast<long>(floor(uPROGRESSSTEP/(_searchstep * _blocksize)));
					resetBuffer(_blocks2load);
					_stubborness = STUBBORNESS;
				}
				//print not found anything
				printf("Window=%d Watermark Not Detected.\t\t\t\t\n", _currentwindow);
				_currentwindow++;
				createPermutation(_currentwindow, _permutationDetect);
				createWatermark(_currentwindow, _watermarkDetect);
			} 
			else 
			{
				// watermark not found 
				// stay in teh same window
				_state = 0; // state stays SEARCH
				_blocks2load = 1;
				resetBuffer(_blocks2load);
				_stubborness -= _searchstep * _blocksize * static_cast<float>(_blocks2load);
			}
		} 
		else if (_state == 1) 
		{
			if (_wdetected == 0 || !(_buffer[0].cci == _buffer[1].cci && _buffer[1].cci == _buffer[2].cci)) 
			{
				// lost synch with watermark - end of watermark 
				_state = 0; 
				printf("Window=%d Watermark Detected [NC=%8.4f WM=%X%X]\t\t\t\n", _currentwindow, _bestcase.nc / 3, _bestcase.cci, _bestcase.load);

				WATERMARKFOUNDINFO wmfound;
				wmfound.cci = _bestcase.cci;
				wmfound.load = _bestcase.load;
				_waterMarkFound.push_back(wmfound);

				_blocks2load = static_cast<long>(floor((sPROGRESSSTEP)/(_searchstep * _blocksize)));
				resetBuffer(_blocks2load);

				_stubborness = xSTUBBORNESS;
				_currentwindow++;

				createPermutation(_currentwindow, _permutationDetect);
				createWatermark(_currentwindow, _watermarkDetect);
				initResults();
			} 
			else 
			{
				// watermark was detected and it is detected
				// in this test as well
				_state = 1; // still DETECTING
				_blocks2load = 1;
				resetBuffer(_blocks2load);
				if (_buffer[0].nc + _buffer[1].nc + _buffer[2].nc > _bestcase.nc) 
				{
					_bestcase.nc = _buffer[0].nc + _buffer[1].nc + _buffer[2].nc;
					_bestcase.dtime = _windowtimeDetect;
					_bestcase.cci = _buffer[(_pbuff + 2) % 3].cci;
					_bestcase.load = _buffer[(_pbuff + 2) % 3].load;

					for (long bit = 0; bit < BITSPERWINDOW; bit++) 
					{
						_bestcase.payload[bit] = _buffer[0].payload[bit] + _buffer[1].payload[bit] + _buffer[2].payload[bit];
					}
				}
			}
		} 

		_windowtimeDetect += _blocksize * _blocks2load * _searchstep;
	}
}

std::string Watermark::stringToHex(const std::string& word)
{
	static const char* const lut = "0123456789ABCDEF";
    size_t len = word.length();

    std::string output;
    output.reserve(2 * len);
    for (size_t i = 0; i < len; ++i)
    {
        const unsigned char c = word[i];
        output.push_back(lut[c >> 4]);
        output.push_back(lut[c & 15]);
    }
    return output;
}

std::string Watermark::hexToString(const std::string& input) const
{

    static const char* const lut = "0123456789ABCDEF";
    size_t len = input.length();
    if (len & 1) throw std::invalid_argument("odd length");

    std::string output;
    output.reserve(len / 2);
    for (size_t i = 0; i < len; i += 2)
    {
        char a = input[i];
        const char* p = std::lower_bound(lut, lut + 16, a);
        if (*p != a) throw std::invalid_argument("not a hex digit");

        char b = input[i + 1];
        const char* q = std::lower_bound(lut, lut + 16, b);
        if (*q != b) throw std::invalid_argument("not a hex digit");

        output.push_back(((p - lut) << 4) | (q - lut));
    }    return output;

}

void Watermark::addWaterMark(AudioSampleBuffer &inputbuffer)
{
	float   mul, invmul;		// signal gain due to watermarking
	mul = DB(OFFSET);
	invmul = 1 / mul;

	float	wmgain;				// gain factor for each watermark

	int		nblknw[NCHMAX];		// how many blocks were not watermarked per channel

	for (int ch = 0; ch < _numChannels; ch++)
	{
		nblknw[ch] = 0;
	}

	//// Initialize subband pointers; for insertion, no frequency shift
	genSubIndex(&_ssBands, _sampleRate, 1); 

	_nBlocks = std::ceil(static_cast<float> ((_bufferSize + NFREQ - 1) / NFREQ));

	for (int block = 0; block < _nBlocks; block++) 
	{	
		// Copy sample data to signal vector xx and compute MCLT
		for (int ch = 0; ch < _numChannels; ch++) 
		{
			if((block * NFREQ) + NFREQ <= inputbuffer.getNumSamples())
			{
				_nfreqBufferApply = inputbuffer.getWritePointer(ch, block * NFREQ);
			}
			else
			{
				_nfreqBufferApply = nullptr;
			}

			if (_nfreqBufferApply != nullptr)
			{				
				
				for (int i = 0 ; i < NFREQ ; i++) 
				{
					_xx[ch][i] = _nfreqBufferApply[i];
				}
				
				 // signal MCLT computation
				fmclt(_Xc[ch], _Xs[ch], _xx[ch], _ya[ch], _h, NFREQ);

				 //Compute watermark flag; set to 1 if block is
				 //"good enough " (no pre-echo) to be watermarked
				 //Modify frequencies only up to FWMAX.
				 //All channels are watermarked with the same chips
				if (getWMFlag(_xx[ch]) == 1.0f) 
				{
					for (int freq = _freqMin; freq <= _freqMax; freq++) 
					{		
						if (_chip[freq][_frame] == 1.0f) 
						{
							wmgain = mul;
						}
						else
						{ 
							wmgain = invmul;
						}

						_Xc[ch][freq] *= wmgain;
						_Xs[ch][freq] *= wmgain;
					}
				} 
				else
				{
					nblknw[ch] += 1;
				}

    			//Inverse MCLT and convert signal vector xx to sample data,
				fimclt(_Xc[ch], _Xs[ch], _xx[ch], _ys[ch], _h, NFREQ);

				if(_firstrun)
				{
					_firstrun = false;
				}
				else
				{
					for (int i = 0 ; i < NFREQ ; i++) 
					{
						_nfreqBufferApply[i]= _xx[ch][i];
					}	
				}
			}
		}

		//// Update frame & window pointers
		_frametime  += _blocktlength;
		_windowtime += _blocktlength;

		// Check if current block should finish current frame
		if (_frametime > _frametlength) 
		{
			// Yes, this was the last block of a frame
			_frametime -= _frametlength; // update counting for next frame
			_frame++; // increment frame counter
		}

		if (_windowtime > WINDOWSIZE) 
		{
			_windowtime -= (float) WINDOWSIZE;
			_window++;	// increment window counter
			_frame = 0;	// reset frame counter

			// Get next set of spread-spectrum chips
			// and mext ramp pattern
			createPermutation(_window, _permutationAdd);
			getChips(_window, _wmbits.xCCI[_window], _wmbits.xLOAD[_window]);
		}
	}

	for (int ch = 0; ch < _numChannels; ch++) 
	{
		if (nblknw[ch] > 0) 
		{
			float dave = nblknw[ch] * 100.0 / _nBlocks;
			//printf("In channel %2d, %4.1f%% of blocks were not watermarked.\n", ch+1, nblknw[ch] * 100.0 / _nBlocks);
			_debugInfo = "In channel %2d, %4.1f%% of blocks were not watermarked.\n", ch+1, nblknw[ch] * 100.0 / _nBlocks;
		}
	}

}

void Watermark::genSubIndex(SSBANDS* ssbands, long fs, float scalefactor) const
{
	int		bit;
	int		fc;
	float	freq, x, y, a;

	// First define bits for SS subbands, on which the chips carrying
	// the CCI bits will be inserted
	//
	// The band from FWCCIMIN to FWMAX is divided in CHIPSPERBLOCK subbands
	//
	//	 FWCCIMIN                                     FWMAX
	//	    |---------|---------|--- . . . ---|---------|
	//          :         :                       :
	//        bit 0     bit 1             bit CHIPSPERBLOCK-1
	//
	// The bands are not uniformly distributed. The parameter SPEEDUP controls
	// the subband with distortion: higgher subbands should be wider as we
	// increase speedup.
	// The main trick: figure out how much the last subband width needs to
	// be amplified; that's the parameter g

    // first determine a from max gradient of warp
	a = 1.7 * log(1 + SPEEDUP * CHIPSPERBLOCK) + 0.5 * (SPEEDUP * CHIPSPERBLOCK); 
	fc = static_cast<int>((FWMIN * scalefactor * 2.0 / fs) * static_cast<float>(NFREQ));
	bit = 0;
    // start at index corresponding to FWCCIMIN
	ssbands->fstart[0] = fc; 
	while (bit < CHIPSPERBLOCK) {
		x = (bit + 1.0) / CHIPSPERBLOCK;
		y = warpFreq(x, a);
        // analog frequency, in Hz
		freq = FWMIN + (FWMAX-FWMIN) * y; 
        // convert to index
		fc = (int) (freq * scalefactor * (2.0 / fs) * NFREQ); 
		ssbands->fend[bit] = fc;
		ssbands->fmiddle[bit] = (ssbands->fstart[bit] + fc) / 2;
		// ??? ssbands->width[bit] = ssbands->fend[bit] - ssbands->fstart[bit] + 1;
		bit++;
		if (bit == CHIPSPERBLOCK) break;
		ssbands->fstart[bit] = fc + 1;
	}
}

float Watermark::warpFreq(float x, float a) const
{
	// ------------------------------------------------------------
	// Routine to warp the frequency scale, to support frequency scale changes
	return((exp(a*x) - 1) / (exp(a) - 1) );

}

void Watermark::createPermutation(int seed, int (&permuatation)[CHIPSPERBLOCK][BITSPERWINDOW])
{
	int count, temp;
	srand(seed);
	for (int bit = 0; bit < CHIPSPERBLOCK; bit++) 
	{
		for (int i = 0; i < BITSPERWINDOW; i++)
		{	
			permuatation[bit][i] = -1;
		}
		for (int i = 0; i < BITSPERWINDOW; i++) 
		{
			temp = rand() % (BITSPERWINDOW - i);
			count = 0;
			for (int j = 0; j < BITSPERWINDOW; j++) 
			{
				if (permuatation[bit][j] == -1) 
				{
					if (count == temp) 
					{
						permuatation[bit][j] = i;
						break;
					} 
					else
					{ 
						count++;
					}
				}
			}
		}
	}
}

void Watermark::getChips(int seed, int cci, int usage)
{
	
	int		bit, i, j, k, freq, rows = 0;
	int		ssbit = 0, frame;
	int		toembed[BITSPERWINDOW];
	int		framesperbit;

	// how the watermark is created?
	// we have hardwired two lookup tables of rows for
	// chess watermarks of length 2 and four bits.
	// we create the larger 6 and 8 bit wide tables  
	// by crossconcatenating the smaller ones.

	// look up tables of size two and four blocks in a row
	int		lut2[2][2] = {   1, 0, 
							 0, 1};
	int		lut4[6][4] = {	1, 1, 0, 0, 
							1, 0, 1, 0, 
							1, 0, 0, 1, 
							0, 0, 1, 1, 
							0, 1, 0, 1, 
							0, 1, 1, 0}; 
	// look up tables of size 8 and 6 blocks in a row
	int		lut8[36][8];
	int		lut6[12][6];
	int		p_lut = 0;

	// create the lookup tables for watermark creation
	// lookup tables of size 6 and 8 are made from 
	// lookup tables fo size 2 and four
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++) 
			for (k = 0; k < 8; k++) 
				if (k < 4) lut8[i*6+j][k] = lut4[i][k];
				else lut8[i*6+j][k] = lut4[j][k-4];
	for (i = 0; i < 6; i++)
		for (j = 0; j < 2; j++) 
			for (k = 0; k < 6; k++) 
				if (k < 4) lut6[i*2+j][k] = lut4[i][k];
				else lut6[i*2+j][k] = lut2[j][k-4];

	// Prepare seed for the random number generator; it
	// should be a function of the CCI bits and the input seed.
	// For now, we just add the 4-bit value of CCI to the seed.

	// Initialize RN generator with the highest bit of cci
	// there are 16 different watermarks 0 <= i <= 15.
	i = cci & 0xF;  
	INITIALIZE_SEED(seed, i);

	// encoding the usage and cci bits
	for (i = 0; i < BITSPERWINDOW; i++) {
		toembed[i] = (usage & (0x1 << i)) >> i;
	}

	// period equals the number of chip blocks per PAYLOADBIT frame
	framesperbit = (int) FRAMESPERWINDOW / BITSPERWINDOW;
	// depending on period the appropriate lookup table is selected
	// first we select its dimension (number of different rows)
	if (framesperbit == 6) rows = 12;
	else if (framesperbit == 8) rows = 36;
	else if (framesperbit == 4) rows = 6;
	// else exit(-1); // this case is an error
	for (bit = 0; bit < CHIPSPERBLOCK; bit++) {
		for (frame = 0; frame < FRAMESPERWINDOW; frame++) {
			// point to different row in the lookup table
			// after each period.
			if (frame % framesperbit == 0) {
				p_lut = static_cast<int>((rand()*rows)/RAND_MAX);
			}
			// get the actual bit to be embedded at [bit,frame]
			if (framesperbit == 6) ssbit = lut6[p_lut][frame%6];
			else if (framesperbit == 8) ssbit = lut8[p_lut][frame%8];
			else if (framesperbit == 4) ssbit = lut4[p_lut][frame%4];
			// define chip from ssbit
			for (freq = _ssBands.fstart[bit]; freq <= _ssBands.fend[bit]; freq++) {
				// --- embedding the actual chip  or toembed = 0x0.
				// chip[freq][frame] = ssbit;
				// --- xoring with the communication channel (toembed) without permutations.
				// chip[freq][frame] = ssbit ^ toembed[(int) frame/period];
				// --- xoring + permutations of the communication channel
				// for each subband (bit from CHIPSPERBLOCK).
				float dave = ssbit ^ toembed[_permutationAdd[bit][(int) frame/framesperbit]];
				_chip[freq][frame] = ssbit ^ toembed[_permutationAdd[bit][(int) frame/framesperbit]];
			}
		}
	}
}

int Watermark::getWMFlag(float* xx)
{
	int scale = NFREQ / NSEC;
	float	er, max_energy, min_energy, temp;

	// The case for not watermarking a block. If abrupt changes in energy are
	// detected in a block, watermark is not inserted to avoid pre-echoes

	max_energy = 0;
	min_energy = static_cast<float>(1e+20);

	for (int i = 0; i < NSEC; i++) 
	{
		temp = 0;
		for (int j = 0; j < scale; j++) 
		{
			temp += xx[i*scale + j] * xx[i*scale + j];
		}
		temp /= scale;
		
		if (temp > max_energy)
		{
			max_energy = temp;
		}

		if (temp < min_energy)
		{
			min_energy = temp;
		}
	}

	er = max_energy/(min_energy + ETHR);

	if (er < ERLIM) return(1); else return(0);
}

void Watermark::createDecoder(long fs)
{
	
	SSBANDS s;
	int i, bit;
	float scalefactor;
    for (i = 0; i < SRFREQ; i++) {
		// computing the frequency shift scaling factor
		scalefactor = 1.-FREQ_RESILIENCE/100.+2.*i*FREQ_RESILIENCE/(100.*(SRFREQ-1));
		genSubIndex(&s, fs, scalefactor);
		for (bit = 0; bit < CHIPSPERBLOCK; bit++) {
			// compute the middle frequency of each subband where
			// a single bit is embedded
		    _ssdecoder.fmiddle[bit][i] = s.fmiddle[bit];
			// determine the range within this subband where
			// the detection is performed
			switch (s.fend[bit] - s.fstart[bit] + 1) {
				case  0: _ssdecoder.cbs[bit][i] =  0; _ssdecoder.cbe[bit][i] = 0; break;
				case  1: _ssdecoder.cbs[bit][i] =  0; _ssdecoder.cbe[bit][i] = 0; break;
				case  2: _ssdecoder.cbs[bit][i] =  0; _ssdecoder.cbe[bit][i] = 0; break;
				case  3: _ssdecoder.cbs[bit][i] =  0; _ssdecoder.cbe[bit][i] = 0; break;
				case  4: _ssdecoder.cbs[bit][i] =  0; _ssdecoder.cbe[bit][i] = 1; break;
				case  5: _ssdecoder.cbs[bit][i] = -1; _ssdecoder.cbe[bit][i] = 1; break;
				case  6: _ssdecoder.cbs[bit][i] =  0; _ssdecoder.cbe[bit][i] = 1; break;
				case  7: _ssdecoder.cbs[bit][i] = -1; _ssdecoder.cbe[bit][i] = 1; break;
				case  8: _ssdecoder.cbs[bit][i] = -1; _ssdecoder.cbe[bit][i] = 2; break;
				case  9: _ssdecoder.cbs[bit][i] = -2; _ssdecoder.cbe[bit][i] = 2; break;
				case 10: _ssdecoder.cbs[bit][i] = -2; _ssdecoder.cbe[bit][i] = 3; break;
				case 11: _ssdecoder.cbs[bit][i] = -2; _ssdecoder.cbe[bit][i] = 2; break;
				case 12: _ssdecoder.cbs[bit][i] = -2; _ssdecoder.cbe[bit][i] = 3; break;
				default: _ssdecoder.cbs[bit][i] = -3; _ssdecoder.cbe[bit][i] = 3;
			}
		}
	}
	// startfreq determines the lowest coeff that takes part in watermarking
	_startfreq = _ssdecoder.fmiddle[0][0] + _ssdecoder.cbs[0][0];;
	// endfreq is the highest coeff that takes part in watermarking
	_endfreq = _ssdecoder.fmiddle[CHIPSPERBLOCK-1][SRFREQ-1] + _ssdecoder.cbe[CHIPSPERBLOCK-1][SRFREQ-1];
}

void Watermark::createPointers()
{
	float scalefactor;
	for (int timeindex = 0; timeindex < SRTIME; timeindex++) 
	{
		// compute the time scaling factor for each test
		scalefactor = 1. - TIME_RESILIENCE / 100. +2. * timeindex * TIME_RESILIENCE / (100. * (SRTIME - 1));

		for (int frame = 0; frame < FRAMESPERWINDOW; frame++) 
		{
			// create the pointers which point to the center of each time interval
			// where a single block of chips is used. for different scaling factors the
			// set of pointers is more or less diffused.
			_pointers[frame][timeindex] = static_cast<long>((1 + (0.2 + frame) * scalefactor * _blocksperwindow / FRAMESPERWINDOW)) / _searchstep;
			// pointers[][] point to the first block within a frame to be integrated.
			// the scheme assumes that 60% of blocks within a frame will be integrated.
			// another assumption is that the basic step is 3 blocks!!!
			// therefore we have 20% of the frame at both sides of the integration area.
			// we get the integer value of the pointer by computing [(x+1)/3]
			// thats how 01->0, 234->1 567->2 etc
			// NOTE: for MCLT=2048 the only thing that should change is 3->6 as the basic
			// step is 6 
		}
	}
}

void Watermark::cepstrumFiltering(float* inout, float* cepstrum)
{
	int i;
	// translating the signal into cepstrum
	// Rico: new cepstral filtering via DST-IV
	fdstiv(inout, NFREQ);
	// lowpass filtering
	for (i = 0; i < CF; i++)
	{
	 	inout[i] = 0;
	}
	// peak removal. it is important to notice that
	// abs(inout[i]) > PM works very bad. cannot explain.
	for (i = CF; i < NFREQ; i++)
	{
		if (inout[i] > PM)
		{
			inout[i] = PM;
		}
	}
	/*
	// highpass filtering. removed due to addition of noise.
	for (i = CM; i < Nbands; i++)
		inout[i] = 0;
	*/
	// translating the signal back to freq
	fdstiv(inout, NFREQ);
	// remove artificial spectral peak from DST-IV
	for (i = 0; i < 4; i++)
	{
	 	inout[i] = 0;
	}
	// dst
}

void Watermark::createWatermark(int seed, int (&watermark)[CHIPSPERBLOCK][FRAMESPERWINDOW][NWATERMARKS]) const
{
	int		i, j, k, rows, frame, bit, p_lut;

	// how the watermark is created?
	// we have hardwired two lookup tables of rows for
	// chess watermarks of length 2 and four bits.
	// we create the larger 6 and 8 bit wide tables  
	// by crossconcatenating the smaller ones.

	// look up tables of size two and four blocks in a row
	int		lut2[2][2] = {  1, 0, 
							0, 1};
	int		lut4[6][4] = {	1, 1, 0, 0, 
							1, 0, 1, 0, 
							1, 0, 0, 1, 
							0, 0, 1, 1, 
							0, 1, 0, 1, 
							0, 1, 1, 0}; 
	// look up tables of size 8 and 6 blocks in a row
	int		lut8[36][8];
	int		lut6[12][6];

	// create the lookup tables for watermark creation
	// lookup tables of size 6 and 8 are made from 
	// lookup tables fo size 2 and four
	for (i = 0; i < 6; i++)
		for (j = 0; j < 6; j++) 
			for (k = 0; k < 8; k++) 
				if (k < 4) lut8[i*6+j][k] = lut4[i][k];
				else lut8[i*6+j][k] = lut4[j][k-4];
	for (i = 0; i < 6; i++)
		for (j = 0; j < 2; j++) 
			for (k = 0; k < 6; k++) 
				if (k < 4) lut6[i*2+j][k] = lut4[i][k];
				else lut6[i*2+j][k] = lut2[j][k-4];

	// depending on period the appropriate lookup table is selected
	// first we select its dimension (number of different rows)
	if (_framesperbit == 6)
	{
		rows = 12;
	}
	else if (_framesperbit == 8)
	{
		rows = 36;
	}
	else if (_framesperbit == 4)
	{
		rows = 6;
	}
	// else exit(-1); // this case is an error

	// create the watermarks
	for (i = 0; i < NWATERMARKS; i++) 
	{
		// all possible watermarks are created. there are 
		// NWATERMARKS=2^cci different watermarks.
		INITIALIZE_SEED(seed, i);

		for (bit = 0; bit < CHIPSPERBLOCK; bit++) 
		{
			for (frame = 0; frame < FRAMESPERWINDOW; frame++) 
			{
				// point to different row in the lookup table
				// after each period.
				if (frame % _framesperbit == 0) 
				{
					p_lut = static_cast<int>((rand()*rows)/RAND_MAX);
				}
				// copy the bits from the selected row of the lookup table
				// into the watermark structure.
				if (_framesperbit == 6)
				{
					watermark[bit][frame][i] = lut6[p_lut][frame%6];
				}
				else if (_framesperbit == 8)
				{
					watermark[bit][frame][i] = lut8[p_lut][frame%8];
				}
				else if (_framesperbit == 4)
				{
					watermark[bit][frame][i] = lut4[p_lut][frame%4];
				}
			}
		}
	}
}

int Watermark::signOfNum(float x)
{
	int retValue = 0;
	if (x >= 0)
		retValue = 1;
	if (x < 0)
		retValue = -1;

	return retValue;
}

float Watermark::correctAmplitudeToGreaterThan90db(float amplitude) const
{
	float retValue = amplitude;

	if(abs(amplitude) < MINAMPLITUDE && amplitude != 0.0)
	{
		if(signOfNum(amplitude) == -1)
		{
			retValue = MINAMPLITUDE * -1.;
		}
		else if(signOfNum(amplitude) == 1)
		{
			retValue = MINAMPLITUDE;
		}
	}

	return retValue;
}

std::string Watermark::getUIUpdate()
{
	std::string waterMarkString;

	for(int i = 0; i < _waterMarkFound.size(); i++)
	{	
		std::stringstream ss;
		ss << std::uppercase << std::hex << _waterMarkFound[i].cci;	
		ss << std::uppercase << std::hex << _waterMarkFound[i].load;
		std::string HexString = hexToString(ss.str());
		waterMarkString.append(HexString);
		waterMarkString.append("");
	}

	return waterMarkString;
}

void Watermark::cleanUp()
{
	//Clean up everything that was new.

	if(_cb.buffer != nullptr)
	{
		_cb.buffer = nullptr;
		delete _cb.buffer;
	}

	if(_cb.ht != nullptr)
	{
		_cb.ht = nullptr;
		delete _cb.ht;
	}
}

int Watermark::computeCorrelationsNormalized(CIRCULAR_BUFFER* cb, long storageblock, int* timeINDEX, int* freqINDEX)
{
	int			counter;
	long        k,p1;
	float		blockcorr, temp;
	int			blockcard;
	float 		maximum, wmean, wstd;

	// initializing the data structure that holds the correlation sums
	for (int timeindex = 0; timeindex < SRTIME; timeindex++) 
	{
		for (int freqindex = 0; freqindex < SRFREQ; freqindex++) 
		{
			for (long w = 0; w < NWATERMARKS; w++) 
			{
			 	_payload[timeindex][freqindex][w].sumsquares = 0.;

				for (long i = 0; i < BITSPERWINDOW; i++) 
				{
					_payload[timeindex][freqindex][w].corr[i][0] = 0.;
					_payload[timeindex][freqindex][w].corr[i][1] = 0.;
					_payload[timeindex][freqindex][w].card[i][0] = 0;
					_payload[timeindex][freqindex][w].card[i][1] = 0;
				}
			}
		}
	}

	// computing all correlations
	// the buffer is processed payloadbit by payloadbit
	for (long i = 0; i < BITSPERWINDOW; i++) 
	{
		// and then for each frame that is a part of the window
		// where the curent payloadbit is spread.
		for (int frame = static_cast<int>(i) * _framesperbit; frame < static_cast<int>(1 + i) * _framesperbit; frame++) 
		{
			// for each time scaling we point to a different subset of blocks
			for (int timeindex = 0; timeindex < SRTIME; timeindex++) 
			{
				// p1 is equal to the first block of storageblock coeffs
				// for the current window of the audio clip
				p1 = cb->pointer + _pointers[frame][timeindex];
				// p1 = p1%cb->length
				if (p1 >= cb->length) 
				{
					p1 -= cb->length;
				}
				// pointing at storageblocks
				p1 = p1 * storageblock;
				// for each frequency scaling 
				for (int freqindex = 0; freqindex < SRFREQ; freqindex++) 
				{
					// point to different subbands. for each subband (bitperblock)
					for (int bit = 0; bit < CHIPSPERBLOCK; bit++) 
					{
						// this is the sum/card for all coeffs that are part of one
						// bit of the watermark
						blockcorr = 0.;
						blockcard = 0;
						// point to the first coeff of the current bit
						k = p1 + _ssdecoder.fmiddle[bit][freqindex] - _startfreq;
						// prevent overflow of the circular buffer.
						// if (k >= cb->length * storageblock) 
						// 	k -= cb->length * storageblock;
						// for each coeff in the current frame and current subband
						for (long kk = _ssdecoder.cbs[bit][freqindex]; kk <= _ssdecoder.cbe[bit][freqindex]; kk++) 
						{
							// if audible add it up
							if (cb->ht[k+kk] > 0) 
							{
								blockcorr += cb->buffer[k+kk];
								blockcard += cb->ht[k+kk];
							}
						}
						// add the sum of coeffs to corr[bit][1] or corr[bit][0]
						// depending on whether the watermark bit is one or zero.
						for (long w = 0; w < NWATERMARKS; w++) 
						{
							// if there are permutations we have to make an extra
							// memory access to put the right sum into the right payload bit.
							_payload[timeindex][freqindex][w].corr[_permutationDetect[bit][i]][_watermarkDetect[bit][frame][w]] += blockcorr;
							_payload[timeindex][freqindex][w].card[_permutationDetect[bit][i]][_watermarkDetect[bit][frame][w]] += blockcard;
						} 
					}
				}
			}
		}
	}

	// find the largest correlation among all watermarks.
	// memorize the watermark (counter) and the attackers scaling
	// (timeINDEX and freqINDEX)
	counter = -1;
	maximum = 0;
	*timeINDEX = -1;
	*freqINDEX = -1;
	// computing sum of correlation squares
	for (long w = 0; w < NWATERMARKS; w++) 
	{
		for (int freqindex = 0; freqindex < SRFREQ; freqindex++) 
		{
			for (int timeindex = 0; timeindex < SRTIME; timeindex++) 
			{
				for (int bit = 0; bit < BITSPERWINDOW; bit++) 
				{
					// correlation test. note that here we compute the
					// ssum(1)/card(1) - sum(0)/card(0) test.
					// it is much faster and as accurate as 
					// the normalized standard corr test.
					// i have a matlab code thatgivesbetter results than
					// the standard corr test (itis in the matlab files that
					// i have sent to you.
					temp = _payload[timeindex][freqindex][w].corr[bit][1]/_payload[timeindex][freqindex][w].card[bit][1];
					temp -= _payload[timeindex][freqindex][w].corr[bit][0]/_payload[timeindex][freqindex][w].card[bit][0];
					// compute squares of correlations
					// payload[timeindex][freqindex][w].sumsquares += temp*temp/(4*OFFSET*OFFSET);
					// or compute abs of correlations
					_payload[timeindex][freqindex][w].sumsquares += static_cast<float>(fabs(temp)/(2*OFFSET));
					// ZEROMEAN computation - subtract the mean AVE line625 insteadof line 623.
					// payload[timeindex][freqindex][w].sumsquares += fabs(temp-ave)/(2*OFFSET);
					// or compute the intact correlations (payload has to be 0x0
					// in order for this option to be used
					// payload[timeindex][freqindex][w].sumsquares += temp/(2*OFFSET);
				}

				_payload[timeindex][freqindex][w].sumsquares /= BITSPERWINDOW;

				// search for the maximum correlation
				if (_payload[timeindex][freqindex][w].sumsquares > maximum) 
				{
					counter = w;
					maximum = _payload[timeindex][freqindex][w].sumsquares;
					*timeINDEX = timeindex;
					*freqINDEX = freqindex;
				}
			}
		}
	}

	// Normalized correlation
	// Computing the mean of all partial correlations
	wmean = 0.;
	for (long w = 0; w < NWATERMARKS; w++)
	{
		for (int timeindex = 0; timeindex < SRTIME; timeindex++)
		{
			for (int freqindex = 0; freqindex < SRFREQ; freqindex++)
			{
				wmean += _payload[timeindex][freqindex][w].sumsquares;
			}
		}
	}
	wmean /= NWATERMARKS * SRFREQ * SRTIME;

	// Computing the standard deviation of all partial correlations
	wstd = 0.;
	for (long w = 0; w < NWATERMARKS; w++)
	{
		for (int timeindex = 0; timeindex < SRTIME; timeindex++)
		{
			for (int freqindex = 0; freqindex < SRFREQ; freqindex++)
			{
				wstd += (_payload[timeindex][freqindex][w].sumsquares - wmean)*(_payload[timeindex][freqindex][w].sumsquares - wmean);
			}
		}
	}
	wstd = static_cast<float>(sqrt(wstd / (NWATERMARKS * SRFREQ * SRTIME)));

	// normalizing the maximal correlation value in order to reflect
	// the standard deviation of computed normalizedcorrelations for all tests 
	_payload[*timeINDEX][*freqINDEX][counter].sumsquares = (maximum - wmean) / wstd;

	return(counter);
}
