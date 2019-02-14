#include"rttimer.h"

#include<string.h>

#if _WIN32
#  include<Windows.h>
#else
#  define UNIX_CLOCK_GETTIME
#  include<unistd.h>

#  include<time.h>
#  include<errno.h>
#endif // _WIN32

//=============================================================================
// Func: _get_raw_time
// Desc: Get the raw time from the system timer.
//-----------------------------------------------------------------------------
static double _get_raw_time(rttimer_p tmr)
{
#if _WIN32
	LARGE_INTEGER li;

	QueryPerformanceCounter(&li);
	return li.QuadPart / tmr->perf_freq;
#elif UNIX_CLOCK_GETTIME
	struct timespec ts;

	if(clock_gettime(CLOCK_MONOTONIC_RAW, &ts) == 0)
	{
		return ts.tv_sec + (time_real_t)ts.tv_nsec * 0.000000001;
	}
	else if(clock_gettime(CLOCK_MONOTONIC, &ts) == 0)
	{
		return ts.tv_sec + (time_real_t)ts.tv_nsec * 0.000000001;
	}
	else
	{
		return (double)clock() / tmr->perf_freq;
	}
#endif
}

void rttimer_init(rttimer_p tmr, int init_pause)
{
#if _WIN32
	LARGE_INTEGER li;

	memset(tmr, 0, sizeof *tmr);

	QueryPerformanceFrequency(&li);
	tmr->perf_freq = (time_real_t)li.QuadPart;
#elif UNIX_CLOCK_GETTIME
	tmr->perf_freq = CLOCKS_PER_SEC;
#endif

	tmr->is_paused = init_pause;
	tmr->start_time = tmr->pause_time = _get_raw_time(tmr);
}

void rttimer_start(rttimer_p tmr)
{
	if(tmr->is_paused)
	{
		time_real_t start_time = _get_raw_time(tmr);
		tmr->is_paused = 0;
		tmr->start_time += start_time - tmr->pause_time;
	}
}

void rttimer_pause(rttimer_p tmr)
{
	if(!tmr->is_paused)
	{
		tmr->is_paused = 1;
		tmr->pause_time = _get_raw_time(tmr);
	}
}

double rttimer_gettime(rttimer_p tmr)
{
	if(tmr->is_paused)
		return tmr->pause_time - tmr->start_time;
	else
		return _get_raw_time(tmr) - tmr->start_time;
}

void rttimer_settime(rttimer_p tmr, time_real_t timeval)
{
	if(tmr->is_paused)
		tmr->start_time = tmr->pause_time - timeval;
	else
		tmr->start_time = _get_raw_time(tmr) - timeval;
}

void rttimer_cleanup(rttimer_p tmr)
{
	memset(tmr, 0, sizeof *tmr);
}