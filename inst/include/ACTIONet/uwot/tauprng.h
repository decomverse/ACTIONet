#ifndef UWOT_TAUPRNG_H
#define UWOT_TAUPRNG_H

class tau_prng {
  long long state0;
  long long state1;
  long long state2;
public:
  tau_prng(long long state0, long long state1, long long state2)
    : state0(state0), state1(state1), state2(state2) {}

  int operator()() {
    state0 = (((state0 & 4294967294LL) << 12) & 0xffffffff) ^
      ((((state0 << 13) & 0xffffffff) ^ state0) >> 19);
    state1 = (((state1 & 4294967288LL) << 4) & 0xffffffff) ^
      ((((state1 << 2) & 0xffffffff) ^ state1) >> 25);
    state2 = (((state2 & 4294967280LL) << 17) & 0xffffffff) ^
      ((((state2 << 3) & 0xffffffff) ^ state2) >> 11);

    return state0 ^ state1 ^ state2;
  }
};

#endif // UWOT_TAUPRNG_H
