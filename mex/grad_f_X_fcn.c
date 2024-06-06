/* This file was automatically generated by CasADi 3.6.5.
 *  It consists of: 
 *   1) content generated by CasADi runtime: not copyrighted
 *   2) template code copied from CasADi source: permissively licensed (MIT-0)
 *   3) user code: owned by the user
 *
 */
#ifdef __cplusplus
extern "C" {
#endif

/* How to prefix internal symbols */
#ifdef CASADI_CODEGEN_PREFIX
  #define CASADI_NAMESPACE_CONCAT(NS, ID) _CASADI_NAMESPACE_CONCAT(NS, ID)
  #define _CASADI_NAMESPACE_CONCAT(NS, ID) NS ## ID
  #define CASADI_PREFIX(ID) CASADI_NAMESPACE_CONCAT(CODEGEN_PREFIX, ID)
#else
  #define CASADI_PREFIX(ID) grad_f_X_fcn_ ## ID
#endif

#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef MATLAB_MEX_FILE
#include <mex.h>
#endif

#ifndef casadi_real
#define casadi_real double
#endif

#ifndef casadi_int
#define casadi_int long long int
#endif

/* Add prefix to internal symbols */
#define casadi_c0 CASADI_PREFIX(c0)
#define casadi_clear CASADI_PREFIX(clear)
#define casadi_copy CASADI_PREFIX(copy)
#define casadi_densify CASADI_PREFIX(densify)
#define casadi_f0 CASADI_PREFIX(f0)
#define casadi_fill CASADI_PREFIX(fill)
#define casadi_from_mex CASADI_PREFIX(from_mex)
#define casadi_project CASADI_PREFIX(project)
#define casadi_s0 CASADI_PREFIX(s0)
#define casadi_s1 CASADI_PREFIX(s1)
#define casadi_s2 CASADI_PREFIX(s2)
#define casadi_s3 CASADI_PREFIX(s3)
#define casadi_s4 CASADI_PREFIX(s4)
#define casadi_s5 CASADI_PREFIX(s5)
#define casadi_s6 CASADI_PREFIX(s6)
#define casadi_s7 CASADI_PREFIX(s7)
#define casadi_to_mex CASADI_PREFIX(to_mex)
#define casadi_trans CASADI_PREFIX(trans)

/* Symbol visibility in DLLs */
#ifndef CASADI_SYMBOL_EXPORT
  #if defined(_WIN32) || defined(__WIN32__) || defined(__CYGWIN__)
    #if defined(STATIC_LINKED)
      #define CASADI_SYMBOL_EXPORT
    #else
      #define CASADI_SYMBOL_EXPORT __declspec(dllexport)
    #endif
  #elif defined(__GNUC__) && defined(GCC_HASCLASSVISIBILITY)
    #define CASADI_SYMBOL_EXPORT __attribute__ ((visibility ("default")))
  #else
    #define CASADI_SYMBOL_EXPORT
  #endif
#endif

void casadi_clear(casadi_real* x, casadi_int n) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = 0;
  }
}

void casadi_copy(const casadi_real* x, casadi_int n, casadi_real* y) {
  casadi_int i;
  if (y) {
    if (x) {
      for (i=0; i<n; ++i) *y++ = *x++;
    } else {
      for (i=0; i<n; ++i) *y++ = 0.;
    }
  }
}

#define CASADI_CAST(x,y) ((x) y)

void casadi_densify(const casadi_real* x, const casadi_int* sp_x, casadi_real* y, casadi_int tr) {
  casadi_int nrow_x, ncol_x, i, el;
  const casadi_int *colind_x, *row_x;
  if (!y) return;
  nrow_x = sp_x[0]; ncol_x = sp_x[1];
  colind_x = sp_x+2; row_x = sp_x+ncol_x+3;
  casadi_clear(y, nrow_x*ncol_x);
  if (!x) return;
  if (tr) {
    for (i=0; i<ncol_x; ++i) {
      for (el=colind_x[i]; el!=colind_x[i+1]; ++el) {
        y[i + row_x[el]*ncol_x] = CASADI_CAST(casadi_real, *x++);
      }
    }
  } else {
    for (i=0; i<ncol_x; ++i) {
      for (el=colind_x[i]; el!=colind_x[i+1]; ++el) {
        y[row_x[el]] = CASADI_CAST(casadi_real, *x++);
      }
      y += nrow_x;
    }
  }
}

void casadi_project(const casadi_real* x, const casadi_int* sp_x, casadi_real* y, const casadi_int* sp_y, casadi_real* w) {
  casadi_int ncol_x, ncol_y, i, el;
  const casadi_int *colind_x, *row_x, *colind_y, *row_y;
  ncol_x = sp_x[1];
  colind_x = sp_x+2; row_x = sp_x + 2 + ncol_x+1;
  ncol_y = sp_y[1];
  colind_y = sp_y+2; row_y = sp_y + 2 + ncol_y+1;
  for (i=0; i<ncol_x; ++i) {
    for (el=colind_y[i]; el<colind_y[i+1]; ++el) w[row_y[el]] = 0;
    for (el=colind_x[i]; el<colind_x[i+1]; ++el) w[row_x[el]] = x[el];
    for (el=colind_y[i]; el<colind_y[i+1]; ++el) y[el] = w[row_y[el]];
  }
}

void casadi_trans(const casadi_real* x, const casadi_int* sp_x, casadi_real* y,
    const casadi_int* sp_y, casadi_int* tmp) {
  casadi_int ncol_x, nnz_x, ncol_y, k;
  const casadi_int* row_x, *colind_y;
  ncol_x = sp_x[1];
  nnz_x = sp_x[2 + ncol_x];
  row_x = sp_x + 2 + ncol_x+1;
  ncol_y = sp_y[1];
  colind_y = sp_y+2;
  for (k=0; k<ncol_y; ++k) tmp[k] = colind_y[k];
  for (k=0; k<nnz_x; ++k) {
    y[tmp[row_x[k]]++] = x[k];
  }
}

void casadi_fill(casadi_real* x, casadi_int n, casadi_real alpha) {
  casadi_int i;
  if (x) {
    for (i=0; i<n; ++i) *x++ = alpha;
  }
}

#ifdef MATLAB_MEX_FILE
casadi_real* casadi_from_mex(const mxArray* p, casadi_real* y, const casadi_int* sp, casadi_real* w) {
  casadi_int nrow, ncol, is_sparse, c, k, p_nrow, p_ncol;
  const casadi_int *colind, *row;
  mwIndex *Jc, *Ir;
  const double* p_data;
  if (!mxIsDouble(p) || mxGetNumberOfDimensions(p)!=2)
    mexErrMsgIdAndTxt("Casadi:RuntimeError",
      "\"from_mex\" failed: Not a two-dimensional matrix of double precision.");
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
  p_nrow = mxGetM(p);
  p_ncol = mxGetN(p);
  is_sparse = mxIsSparse(p);
  Jc = 0;
  Ir = 0;
  if (is_sparse) {
    Jc = mxGetJc(p);
    Ir = mxGetIr(p);
  }
  p_data = (const double*)mxGetData(p);
  if (p_nrow==1 && p_ncol==1) {
    casadi_int nnz;
    double v = is_sparse && Jc[1]==0 ? 0 : *p_data;
    nnz = sp[ncol];
    casadi_fill(y, nnz, v);
  } else {
    casadi_int tr = 0;
    if (nrow!=p_nrow || ncol!=p_ncol) {
      tr = nrow==p_ncol && ncol==p_nrow && (nrow==1 || ncol==1);
      if (!tr) mexErrMsgIdAndTxt("Casadi:RuntimeError",
                                 "\"from_mex\" failed: Dimension mismatch. "
                                 "Expected %d-by-%d, got %d-by-%d instead.",
                                 nrow, ncol, p_nrow, p_ncol);
    }
    if (is_sparse) {
      if (tr) {
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]+c*nrow]=0;
        for (c=0; c<p_ncol; ++c)
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[c+Ir[k]*p_ncol] = p_data[k];
        for (c=0; c<ncol; ++c)
          for (k=colind[c]; k<colind[c+1]; ++k) y[k] = w[row[k]+c*nrow];
      } else {
        for (c=0; c<ncol; ++c) {
          for (k=colind[c]; k<colind[c+1]; ++k) w[row[k]]=0;
          for (k=Jc[c]; k<(casadi_int) Jc[c+1]; ++k) w[Ir[k]]=p_data[k];
          for (k=colind[c]; k<colind[c+1]; ++k) y[k]=w[row[k]];
        }
      }
    } else {
      for (c=0; c<ncol; ++c) {
        for (k=colind[c]; k<colind[c+1]; ++k) {
          y[k] = p_data[row[k]+c*nrow];
        }
      }
    }
  }
  return y;
}

#endif

#define casadi_to_double(x) ((double) x)

#ifdef MATLAB_MEX_FILE
mxArray* casadi_to_mex(const casadi_int* sp, const casadi_real* x) {
  casadi_int nrow, ncol, c, k;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int nnz;
#endif
  const casadi_int *colind, *row;
  mxArray *p;
  double *d;
#ifndef CASADI_MEX_NO_SPARSE
  casadi_int i;
  mwIndex *j;
#endif /* CASADI_MEX_NO_SPARSE */
  nrow = *sp++;
  ncol = *sp++;
  colind = sp;
  row = sp+ncol+1;
#ifndef CASADI_MEX_NO_SPARSE
  nnz = sp[ncol];
  if (nnz!=nrow*ncol) {
    p = mxCreateSparse(nrow, ncol, nnz, mxREAL);
    for (i=0, j=mxGetJc(p); i<=ncol; ++i) *j++ = *colind++;
    for (i=0, j=mxGetIr(p); i<nnz; ++i) *j++ = *row++;
    if (x) {
      d = (double*)mxGetData(p);
      for (i=0; i<nnz; ++i) *d++ = casadi_to_double(*x++);
    }
    return p;
  }
#endif /* CASADI_MEX_NO_SPARSE */
  p = mxCreateDoubleMatrix(nrow, ncol, mxREAL);
  if (x) {
    d = (double*)mxGetData(p);
    for (c=0; c<ncol; ++c) {
      for (k=colind[c]; k<colind[c+1]; ++k) {
        d[row[k]+c*nrow] = casadi_to_double(*x++);
      }
    }
  }
  return p;
}

#endif

#ifndef CASADI_PRINTF
#ifdef MATLAB_MEX_FILE
  #define CASADI_PRINTF mexPrintf
#else
  #define CASADI_PRINTF printf
#endif
#endif

static const casadi_int casadi_s0[5] = {4, 1, 0, 1, 1};
static const casadi_int casadi_s1[4] = {0, 1, 3, 5};
static const casadi_int casadi_s2[6] = {4, 1, 0, 2, 1, 3};
static const casadi_int casadi_s3[5] = {4, 1, 0, 1, 3};
static const casadi_int casadi_s4[13] = {4, 4, 0, 2, 4, 5, 6, 0, 1, 1, 3, 2, 3};
static const casadi_int casadi_s5[13] = {4, 4, 0, 1, 3, 4, 6, 0, 0, 1, 2, 1, 3};
static const casadi_int casadi_s6[8] = {4, 1, 0, 4, 0, 1, 2, 3};
static const casadi_int casadi_s7[5] = {1, 1, 0, 1, 0};

static const casadi_real casadi_c0[4] = {1., 0., 1., 1.};

/* grad_f_X_fcn:(i0[4],i1,i2,i3,i4,i5,i6,i7,i8)->(o0[4x4,6nz]) */
static int casadi_f0(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem) {
  casadi_int i;
  casadi_real *rr, *ss;
  const casadi_int *cii;
  const casadi_real *cs;
  casadi_real *w0=w+4, *w1=w+10, w2, w4, w7, *w8=w+17, *w9=w+21, *w10=w+23, *w11=w+25;
  /* #0: @0 = zeros(4x4,6nz) */
  casadi_clear(w0, 6);
  /* #1: @1 = [1, 0, 1, 1] */
  casadi_copy(casadi_c0, 4, w1);
  /* #2: @2 = input[1][0] */
  w2 = arg[1] ? arg[1][0] : 0;
  /* #3: @3 = 00 */
  /* #4: @4 = 1 */
  w4 = 1.;
  /* #5: @5 = 00 */
  /* #6: @6 = 00 */
  /* #7: @7 = vertcat(@3, @4, @5, @6) */
  rr=(&w7);
  *rr++ = w4;
  /* #8: @7 = (@2*@7) */
  w7  = (w2*w7);
  /* #9: @8 = dense(@7) */
  casadi_densify((&w7), casadi_s0, w8, 0);
  /* #10: @1 = (@1+@8) */
  for (i=0, rr=w1, cs=w8; i<4; ++i) (*rr++) += (*cs++);
  /* #11: (@0[0, 1, 3, 5] = @1) */
  for (cii=casadi_s1, rr=w0, ss=w1; cii!=casadi_s1+4; ++cii, ++ss) rr[*cii] = *ss;
  /* #12: @7 = ones(4x1,1nz) */
  w7 = 1.;
  /* #13: @9 = project(@7) */
  casadi_project((&w7), casadi_s0, w9, casadi_s2, w);
  /* #14: @3 = 00 */
  /* #15: @5 = 00 */
  /* #16: @6 = 00 */
  /* #17: @7 = 1 */
  w7 = 1.;
  /* #18: @4 = vertcat(@3, @5, @6, @7) */
  rr=(&w4);
  *rr++ = w7;
  /* #19: @2 = (@2*@4) */
  w2 *= w4;
  /* #20: @10 = project(@2) */
  casadi_project((&w2), casadi_s3, w10, casadi_s2, w);
  /* #21: @9 = (@9+@10) */
  for (i=0, rr=w9, cs=w10; i<2; ++i) (*rr++) += (*cs++);
  /* #22: @10 = @9[:2] */
  for (rr=w10, ss=w9+0; ss!=w9+2; ss+=1) *rr++ = *ss;
  /* #23: (@0[2:6:2] = @10) */
  for (rr=w0+2, ss=w10; rr!=w0+6; rr+=2) *rr = *ss++;
  /* #24: @11 = @0' */
  casadi_trans(w0,casadi_s5, w11, casadi_s4, iw);
  /* #25: output[0][0] = @11 */
  casadi_copy(w11, 6, res[0]);
  return 0;
}

CASADI_SYMBOL_EXPORT int grad_f_X_fcn(const casadi_real** arg, casadi_real** res, casadi_int* iw, casadi_real* w, int mem){
  return casadi_f0(arg, res, iw, w, mem);
}

CASADI_SYMBOL_EXPORT int grad_f_X_fcn_alloc_mem(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT int grad_f_X_fcn_init_mem(int mem) {
  return 0;
}

CASADI_SYMBOL_EXPORT void grad_f_X_fcn_free_mem(int mem) {
}

CASADI_SYMBOL_EXPORT int grad_f_X_fcn_checkout(void) {
  return 0;
}

CASADI_SYMBOL_EXPORT void grad_f_X_fcn_release(int mem) {
}

CASADI_SYMBOL_EXPORT void grad_f_X_fcn_incref(void) {
}

CASADI_SYMBOL_EXPORT void grad_f_X_fcn_decref(void) {
}

CASADI_SYMBOL_EXPORT casadi_int grad_f_X_fcn_n_in(void) { return 9;}

CASADI_SYMBOL_EXPORT casadi_int grad_f_X_fcn_n_out(void) { return 1;}

CASADI_SYMBOL_EXPORT casadi_real grad_f_X_fcn_default_in(casadi_int i) {
  switch (i) {
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* grad_f_X_fcn_name_in(casadi_int i) {
  switch (i) {
    case 0: return "i0";
    case 1: return "i1";
    case 2: return "i2";
    case 3: return "i3";
    case 4: return "i4";
    case 5: return "i5";
    case 6: return "i6";
    case 7: return "i7";
    case 8: return "i8";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const char* grad_f_X_fcn_name_out(casadi_int i) {
  switch (i) {
    case 0: return "o0";
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* grad_f_X_fcn_sparsity_in(casadi_int i) {
  switch (i) {
    case 0: return casadi_s6;
    case 1: return casadi_s7;
    case 2: return casadi_s7;
    case 3: return casadi_s7;
    case 4: return casadi_s7;
    case 5: return casadi_s7;
    case 6: return casadi_s7;
    case 7: return casadi_s7;
    case 8: return casadi_s7;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT const casadi_int* grad_f_X_fcn_sparsity_out(casadi_int i) {
  switch (i) {
    case 0: return casadi_s4;
    default: return 0;
  }
}

CASADI_SYMBOL_EXPORT int grad_f_X_fcn_work(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 13;
  if (sz_res) *sz_res = 2;
  if (sz_iw) *sz_iw = 5;
  if (sz_w) *sz_w = 31;
  return 0;
}

CASADI_SYMBOL_EXPORT int grad_f_X_fcn_work_bytes(casadi_int *sz_arg, casadi_int* sz_res, casadi_int *sz_iw, casadi_int *sz_w) {
  if (sz_arg) *sz_arg = 13*sizeof(const casadi_real*);
  if (sz_res) *sz_res = 2*sizeof(casadi_real*);
  if (sz_iw) *sz_iw = 5*sizeof(casadi_int);
  if (sz_w) *sz_w = 31*sizeof(casadi_real);
  return 0;
}

#ifdef MATLAB_MEX_FILE
void mex_grad_f_X_fcn(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  casadi_int i;
  int mem;
  casadi_real w[49];
  casadi_int iw[5];
  const casadi_real* arg[13] = {0};
  casadi_real* res[2] = {0};
  if (argc>9) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"grad_f_X_fcn\" failed. Too many input arguments (%d, max 9)", argc);
  if (resc>1) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"grad_f_X_fcn\" failed. Too many output arguments (%d, max 1)", resc);
  if (--argc>=0) arg[0] = casadi_from_mex(argv[0], w, casadi_s6, w+18);
  if (--argc>=0) arg[1] = casadi_from_mex(argv[1], w+4, casadi_s7, w+18);
  if (--argc>=0) arg[2] = casadi_from_mex(argv[2], w+5, casadi_s7, w+18);
  if (--argc>=0) arg[3] = casadi_from_mex(argv[3], w+6, casadi_s7, w+18);
  if (--argc>=0) arg[4] = casadi_from_mex(argv[4], w+7, casadi_s7, w+18);
  if (--argc>=0) arg[5] = casadi_from_mex(argv[5], w+8, casadi_s7, w+18);
  if (--argc>=0) arg[6] = casadi_from_mex(argv[6], w+9, casadi_s7, w+18);
  if (--argc>=0) arg[7] = casadi_from_mex(argv[7], w+10, casadi_s7, w+18);
  if (--argc>=0) arg[8] = casadi_from_mex(argv[8], w+11, casadi_s7, w+18);
  --resc;
  res[0] = w+12;
  grad_f_X_fcn_incref();
  mem = grad_f_X_fcn_checkout();
  i = grad_f_X_fcn(arg, res, iw, w+18, mem);
  if (i) mexErrMsgIdAndTxt("Casadi:RuntimeError","Evaluation of \"grad_f_X_fcn\" failed.");
  grad_f_X_fcn_release(mem);
  grad_f_X_fcn_decref();
  if (res[0]) resv[0] = casadi_to_mex(casadi_s4, res[0]);
}
#endif

casadi_int main_grad_f_X_fcn(casadi_int argc, char* argv[]) {
  casadi_int j;
  casadi_real* a;
  const casadi_real* r;
  casadi_int flag;
  casadi_int iw[5];
  casadi_real w[49];
  const casadi_real* arg[13];
  casadi_real* res[2];
  arg[0] = w+0;
  arg[1] = w+4;
  arg[2] = w+5;
  arg[3] = w+6;
  arg[4] = w+7;
  arg[5] = w+8;
  arg[6] = w+9;
  arg[7] = w+10;
  arg[8] = w+11;
  res[0] = w+12;
  a = w;
  for (j=0; j<12; ++j) if (scanf("%lg", a++)<=0) return 2;
  flag = grad_f_X_fcn(arg, res, iw, w+18, 0);
  if (flag) return flag;
  r = w+12;
  for (j=0; j<6; ++j) CASADI_PRINTF("%g ", *r++);
  CASADI_PRINTF("\n");
  return 0;
}


#ifdef MATLAB_MEX_FILE
void mexFunction(int resc, mxArray *resv[], int argc, const mxArray *argv[]) {
  char buf[13];
  int buf_ok = argc > 0 && !mxGetString(*argv, buf, sizeof(buf));
  if (!buf_ok) {
    mex_grad_f_X_fcn(resc, resv, argc, argv);
    return;
  } else if (strcmp(buf, "grad_f_X_fcn")==0) {
    mex_grad_f_X_fcn(resc, resv, argc-1, argv+1);
    return;
  }
  mexErrMsgTxt("First input should be a command string. Possible values: 'grad_f_X_fcn'");
}
#endif
CASADI_SYMBOL_EXPORT int main(int argc, char* argv[]) {
  if (argc<2) {
    /* name error */
  } else if (strcmp(argv[1], "grad_f_X_fcn")==0) {
    return main_grad_f_X_fcn(argc-2, argv+2);
  }
  fprintf(stderr, "First input should be a command string. Possible values: 'grad_f_X_fcn'\nNote: you may use function.generate_input to create a command string.\n");
  return 1;
}
#ifdef __cplusplus
} /* extern "C" */
#endif
