#pragma once

#ifndef OFFC_DELPHI_EXPORTS
#define OFFC_DELPHI_API __declspec(dllexport)
#else
#define OFFC_DELPHI_API __declspec(dllimport)
#endif

OFFC_DELPHI_API void mm_ikj_seq(double *a, double *b, double *c, int m, int k, int n);
OFFC_DELPHI_API void mm_ijk_vec(double *a, double *b, double *c, int m, int k, int n);
OFFC_DELPHI_API void mm_ikj_par(double *a, double *b, double *c, int m, int k, int n, int t);
OFFC_DELPHI_API void mm_ikj_parvec(double *a, double *b, double *c, int m, int k, int n, int t);
OFFC_DELPHI_API void mm_strassen(double *a, double *b, double *c, int m, int k, int n);
