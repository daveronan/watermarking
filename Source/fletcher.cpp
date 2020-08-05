#include "hthres.h"
#include "watermark.h"

extern "C"
{
// Lookup tables
	float		fm_lookup[NFREQ];		// precomputed fletcher munson limits

	void create_fletcher_munson (long fs) 
	{
		int i;
		for(i = 0; i < NFREQ; i++) 
			fm_lookup[i] = thrabs((0.5 + i) * static_cast<float>(fs) * 0.5 / static_cast<float>(NFREQ));
	}
}