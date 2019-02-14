//=============================================================================
// Module: rttimer
// Desc: Real Time Timer, designed for precise
// Author: 0xAA55
//-----------------------------------------------------------------------------

#ifndef _RTTIMER_H_
#define _RTTIMER_H_

#include<stdint.h>

typedef double time_real_t, *time_real_p;

typedef struct rttimer_s
{
	time_real_t start_time;
	time_real_t pause_time;
	int is_paused;

	time_real_t perf_freq;
}rttimer_t, *rttimer_p;

//=============================================================================
// Func: rttimer_init
// Desc: Initialize the timer
//-----------------------------------------------------------------------------
void rttimer_init(rttimer_p tmr, int init_pause);

//=============================================================================
// Func: rttimer_start
// Desc: Start or resume the timer
//-----------------------------------------------------------------------------
void rttimer_start(rttimer_p tmr);

//=============================================================================
// Func: rttimer_pause
// Desc: Pause.
//-----------------------------------------------------------------------------
void rttimer_pause(rttimer_p tmr);

//=============================================================================
// Func: rttimer_gettime
// Desc: Get the current time in SECONDS.
//-----------------------------------------------------------------------------
double rttimer_gettime(rttimer_p tmr);

//=============================================================================
// Func: rttimer_settime
// Desc: Set the value of the timer.
//-----------------------------------------------------------------------------
void rttimer_settime(rttimer_p tmr, time_real_t timeval);

//=============================================================================
// Func: rttimer_cleanup
// Desc: Garbage collection
//-----------------------------------------------------------------------------
void rttimer_cleanup(rttimer_p tmr);


#endif // !_REALTIME_TIMER_