/* This file was automatically generated by CasADi.
   The CasADi copyright holders make no ownership claim of its contents. */
#ifdef __cplusplus
extern "C" {
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

int UKF_Func(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem);
int UKF_Func_alloc_mem(void);
int UKF_Func_init_mem(int mem);
void UKF_Func_free_mem(int mem);
int UKF_Func_checkout(void);
void UKF_Func_release(int mem);
void UKF_Func_incref(void);
void UKF_Func_decref(void);
casadi_int UKF_Func_n_out(void);
casadi_int UKF_Func_n_in(void);
casadi_real UKF_Func_default_in(casadi_int i);
const char* UKF_Func_name_in(casadi_int i);
const char* UKF_Func_name_out(casadi_int i);
const casadi_int* UKF_Func_sparsity_in(casadi_int i);
const casadi_int* UKF_Func_sparsity_out(casadi_int i);
int UKF_Func_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w);
#ifdef __cplusplus
} /* extern "C" */
#endif
