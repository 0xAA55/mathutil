




%macro Function 3-4

%1 %2 %3

%endmacro

Function real_t r_rnd (uint32_t *p_seed)
{
	return (real_t)(((*p_seed = *p_seed * 214013L + 2531011L) >> 16) & 0x7fff) / (real_t)RAND_MAX;
}
