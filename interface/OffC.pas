unit OffC;

interface

uses
  System.SysUtils;

const
  DLLNAME = 'offc_delphi.dll';

procedure mm_ikj_seq(A, B, C: PDouble; M, K, N: Integer); cdecl; external DLLNAME name 'mm_ikj_seq';
procedure mm_ijk_vec(A, B, C: PDouble; M, K, N: Integer); cdecl; external DLLNAME name 'mm_ijk_vec';
procedure mm_ikj_par(A, B, C: PDouble; M, K, N, T: Integer); cdecl; external DLLNAME name 'mm_ikj_par';
procedure mm_ikj_parvec(A, B, C: PDouble; M, K, N, T: Integer); cdecl; external DLLNAME name 'mm_ikj_parvec';
procedure mm_strassen(A, B, C: PDouble; M, K, N, TH: Integer); cdecl; external DLLNAME name 'mm_strassen';

implementation

end.
