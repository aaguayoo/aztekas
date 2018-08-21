void RIEMANN();
void KH();
void JET();
void SPH_ACC();
void INI_CUSTOM();

void OUTFLOW(double *B);
void PERIODIC(double *B, int r, int l, int u, int d, int f, int b);
void REFLECTION(double *B, int r, int l, int u, int d, int f, int b);
void JET_LAUNCH(double *B);
void IN_OUT_BOUND(double *B);
void BOUND_CUSTOM();

void RESTART();
