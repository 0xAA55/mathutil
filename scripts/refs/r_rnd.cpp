%line 1+1 r_rnd.c

%line 7+1 r_rnd.c

Function real_t r_rnd (uint32_t *p_seed)
{
 return (real_t)(((*p_seed = *p_seed * 214013L + 2531011L) >> 16) & 0x7fff) / (real_t)RAND_MAX
}
