unit Bridge;

interface

uses
  System.SysUtils;

const
  DLLNAME = 'bridge_delphi.dll';

procedure MMSeqGustavson(A, B, C: PDouble; M, K, N: Integer); cdecl; external DLLNAME name 'mm_seq_gustavson';
procedure MMParGustavson(A, B, C: PDouble; M, K, N, T: Integer); cdecl; external DLLNAME name 'mm_par_gustavson';
procedure MMVec(A, B, C: PDouble; M, K, N: Integer); cdecl; external DLLNAME name 'mm_vec';

implementation

end.
