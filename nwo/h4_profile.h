

enum h4_transitions_e {
  h4_TMM = 0,
  h4_TMI = 1,
  h4_TMD = 2,
  h4_TIM = 3,
  h4_TII = 4,
  h4_TID = 5, 
  h4_TDM = 6, 
  h4_TDI = 7,
  h4_TDD = 8
};
#define h4_NTRANSITIONS 9


typedef struct {
  int     M;
  float **t;   // transitions.     [0..M][0..8]
  float **e;   // match emissions. [0..M][0..K-1]


} H4_PROFILE;
