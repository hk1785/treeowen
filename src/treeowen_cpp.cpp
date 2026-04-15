// [[Rcpp::depends(Rcpp)]]
// [[Rcpp::plugins(cpp11)]]

#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <numeric>
#include <cstdint>
#include <limits>


// ── Package sentinel ─────────────────────────────────────────────────────────
// Called by .onLoad to confirm the shared library loaded successfully.
// Returns R NULL; never throws.
// [[Rcpp::export]]
SEXP treeowen_ping() {
  return R_NilValue;
}


// treeowen_cpp.cpp  —  Rcpp C++ accelerators for TreeOwen v5.12
// [[Rcpp::plugins(cpp11)]]
//
// ── v5.11 → v5.12 변경사항 ───────────────────────────────────────────────────
//
// [Fix-miss] missing 방향 처리 수정 (R + C++ 동시 적용)
//   - 기존: missL=0을 "miss 방향 없음" sentinel로 사용
//     → miss 방향이 yes/no child와 일치하는 경우(XGBoost/LightGBM 일반 케이스)
//       missL=0으로 저장되어 NA x_j를 가중평균으로 처리 → 효율성 공리 위반
//   - 수정: has_miss_dir 벡터 추가 (R logical → C++ uint8_t)
//     · has_miss_dir[i]=TRUE/1: missL[i]가 항상 유효한 miss 방향을 가리킴
//     · has_miss_dir[i]=FALSE/0: miss 방향 정보 없음 → 가중평균 유지
//   - 변경 범위:
//     · R  .prepare_dp_model: missL 저장 조건에서 != yl_row/nr_row 제거,
//                              has_miss_dir 벡터 생성 및 trees[[tt]]에 포함
//     · R  .eval_v_masks_forest_chunked_prepared: mrL!=0 → has_miss_dir[loc]
//     · C++ TreeDP: has_miss_dir 필드 추가
//     · C++ list_to_tree: has_miss_dir 읽기 추가
//     · C++ eval_one_mask, eval_one_flag, eval_one_tree_flag:
//           missL!=0 → has_miss_dir[loc-1] 조건으로 교체
//   - 수치 영향: NA feature 값이 없는 경우 동일 결과. NA가 있는 경우에만 차이.
//   - 기존(v5.10): v0/v1 계산 시 relevant tree만 delta 방식으로 재계산
//   - 수정(v5.11): eval_one_flag() 로 전체 forest full scan으로 복원
//   - 이유: treeowen_exact_enum은 XGBoost/LightGBM(depth=3, p 大) 환경에서
//     역인덱스 효과가 크나, 본 시뮬레이션 목적상 단순한 구현을 유지
//   - 주의: inner_shapley_exact_cpp, inner_mc_perm_cpp의 역인덱스는 유지
//     (FeatIndex 구조체 및 build_relevant_mask 제거, fidx 참조 제거는
//      inner_shapley_enum_one_group_cpp 범위에만 적용)
//   - build_relevant_mask() 헬퍼 함수도 이 함수만 사용하므로 함께 제거
//
//
// ── v5.10 → v5.11 변경사항 (유지) ───────────────────────────────────────────
//
// [Revert-C07] inner_shapley_enum_one_group_cpp: FeatIndex 역인덱스 제거
//
// ── v5.9 → v5.10 변경사항 (유지) ────────────────────────────────────────────
//
// [Fix-C07] inner_shapley_enum_one_group_cpp: FeatIndex 역인덱스 → v5.11에서 제거
//
// [Fix-C08] eval_v_forest_chunk_exact_cpp: R에서 미연결 → 제거 (dead export)
//   - R fallback 경로가 이 함수를 호출하지 않음 (직접 R 구현 사용)
//   - 불필요한 export 제거로 컴파일 크기 감소
//
// ── v5.8 → v5.9 변경사항 (유지) ─────────────────────────────────────────────
//
// [Fix-C01] inner_shapley_exact_cpp (m<=30): Gray code 합산 O(B) → delta 방식 O(|affected|)
//   - prev_total에서 delta만 갱신: total += (new_val - tree_root_cache[tt])
//   - Ranger 500 trees 기준 Gray code 합산 ~500배 비용 제거
//
// [Fix-C02] inner_mc_perm_cpp: contrib_abs reserve 부족 방어
//   - adaptive=TRUE 시 n_perm 이상 실행 가능 → reserve(n_perm*4+8)로 증가
//
// [Fix-C03] eff_maxj 통일: max(R전달 maxj, cache_ref.maxj) 보수적 처리
//
// [Fix-C04] eval_v_forest_chunk_exact_cpp: 명시적 maxj 파라미터 (→ v5.10에서 제거)
//
// [Fix-C05] list_to_tree: std::vector<bool> → std::vector<uint8_t>
//
// [Fix-C06] inner_shapley_exact_cpp: eff_maxj 범위 검사 통일
//
// ── v5.7.1 → v5.8 변경사항 (유지) ───────────────────────────────────────────
//
// [Perf-3] inner_mc_perm_cpp: adaptive SE 조기 종료 C++ 이식
// [Perf-4] TreesCache.max_nn 캐싱
//
// ── v5.7 → v5.7.1 버그 수정 (유지) ──────────────────────────────────────────
//
// [Fix-1] list_to_tree: seen 배열 OOB 수정
// [Fix-2] FeatIndex.build: fidx.maxj 명시적 전달
// [Fix-3] Gray code checkUserInterrupt
// [Fix-4] inner_mc_perm_cpp: od 루프 중복 계산 제거
// [Fix-5] T_bits.reserve 개선
// [Fix-7] .try_load_cpp() 정규식 → R 파일 수정
// [Fix-8] d_minus >= 30 overflow 방어
//
// ── 불변 사항 ─────────────────────────────────────────────────────────────────
//   - 수치 로직 전부 유지 (정확도 불변)
//   - 인덱스 규약: 1-based local index, fcol_local 0=leaf 유지
//   - 외부 export 함수: eval_v_forest_chunk_exact_cpp 제거 (v5.10)
// =============================================================================


using namespace Rcpp;


// =============================================================================
// §0  공통 헬퍼: popcount64 / portable_ctz64 / safe_shapley_w
// =============================================================================

static inline int popcount64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_popcountll(x);
#elif defined(_MSC_VER) && defined(_WIN64)
  return (int)__popcnt64(x);
#else
  x -= (x >> 1) & UINT64_C(0x5555555555555555);
  x  = (x & UINT64_C(0x3333333333333333)) + ((x >> 2) & UINT64_C(0x3333333333333333));
  x  = (x + (x >> 4)) & UINT64_C(0x0f0f0f0f0f0f0f0f);
  return (int)((x * UINT64_C(0x0101010101010101)) >> 56);
#endif
}

static inline int portable_ctz64(uint64_t x) {
  // x == 0 은 호출 전 보장되어야 함 (Gray code diff는 항상 nonzero)
#if defined(__GNUC__) || defined(__clang__)
  return __builtin_ctzll(x);
#elif defined(_MSC_VER)
  unsigned long idx = 0;
#  ifdef _WIN64
  _BitScanForward64(&idx, x);
#  else
  if (_BitScanForward(&idx, (unsigned long)(x & 0xFFFFFFFFULL))) return (int)idx;
  _BitScanForward(&idx, (unsigned long)(x >> 32)); return (int)idx + 32;
#  endif
  return (int)idx;
#else
  int cnt = 0;
  while (!(x & 1ULL)) { x >>= 1; ++cnt; }
  return cnt;
#endif
}

static inline double safe_shapley_w(int s, int m) {
  if (m <= 0 || s < 0 || s >= m) return 0.0;
  return std::exp(R::lgammafn(s + 1.0) + R::lgammafn(m - s) - R::lgammafn(m + 1.0));
}


// =============================================================================
// §1  트리 구조체 & 변환
// =============================================================================

struct TreeDP {
  std::vector<int>    post_local;
  std::vector<int>    fcol_local;   // 1-based feature col, 0 = leaf node
  std::vector<double> split_val;
  std::vector<double> pred_val;
  std::vector<int>    nkids;
  std::vector<int>    kid1, kid2;
  std::vector<double> p1, p2;
  std::vector<int>    yesL, noL, missL;
  std::vector<uint8_t> has_miss_dir; // [Fix-miss] 1 = 이 노드에 유효한 miss 방향 있음
  int                 root_local;
  bool                has_missing;
  // [II-3] 이 tree에서 split에 실제 사용된 feature (1-based, 중복 없음)
  std::vector<int>    used_feats;
};

static TreeDP list_to_tree(const List& tr) {
  TreeDP t;
  t.post_local  = as<std::vector<int>>(tr["post_local"]);
  t.fcol_local  = as<std::vector<int>>(tr["fcol_local"]);
  t.split_val   = as<std::vector<double>>(tr["split_local"]);
  t.pred_val    = as<std::vector<double>>(tr["pred_local"]);
  t.nkids       = as<std::vector<int>>(tr["nkids"]);
  t.kid1        = as<std::vector<int>>(tr["kid1"]);
  t.kid2        = as<std::vector<int>>(tr["kid2"]);
  t.p1          = as<std::vector<double>>(tr["p1"]);
  t.p2          = as<std::vector<double>>(tr["p2"]);
  t.yesL        = as<std::vector<int>>(tr["yesL"]);
  t.noL         = as<std::vector<int>>(tr["noL"]);
  t.missL       = as<std::vector<int>>(tr["missL"]);
  // [Fix-miss] has_miss_dir: R logical vector → uint8_t vector
  {
    LogicalVector hmd = tr["has_miss_dir"];
    t.has_miss_dir.resize(hmd.size());
    for (int i = 0; i < hmd.size(); ++i)
      t.has_miss_dir[i] = (hmd[i] == TRUE) ? 1 : 0;
  }
  t.root_local  = as<int>(tr["root_local"]);
  t.has_missing = as<bool>(tr["has_missing"]);

  // [II-3] used_feats 수집
  // [Fix-1] seen 배열은 실제 feature index 최대값 기준으로 할당해야 함.
  //         fcol_local.size() = node 수이고, feature index j는 최대 maxj (전체 feature 수)
  //         까지 될 수 있으므로 seen(fcol_local.size()) 로 하면 OOB.
  // [Fix-C05] vector<bool> → vector<uint8_t> 일관성 수정
  {
    int max_feat = 0;
    for (int j : t.fcol_local)
      if (j > max_feat) max_feat = j;
    if (max_feat > 0) {
      std::vector<uint8_t> seen(max_feat + 1, 0);  // [Fix-C05] bool → uint8_t
      for (int j : t.fcol_local) {
        if (j > 0 && !seen[j]) {
          seen[j] = 1;
          t.used_feats.push_back(j);
        }
      }
    }
  }
  return t;
}


// =============================================================================
// §2  FeatIndex: feature (1-based) → tree index (0-based) 역인덱스
//
// [Fix-2] build() 에 명시적으로 전달받은 maxj 사용.
//         trees_for() 범위 검사 유지 (방어적 반환).
// =============================================================================

struct FeatIndex {
  int maxj;  // 역인덱스 커버 범위: feature 1..maxj
  std::vector<std::vector<int>> feat_to_trees;

  FeatIndex() : maxj(0) {}

  void build(const std::vector<TreeDP>& trees, int mj) {
    maxj = mj;
    if (maxj <= 0) return;
    feat_to_trees.assign(maxj, std::vector<int>());
    for (int tt = 0; tt < (int)trees.size(); ++tt) {
      for (int f : trees[tt].used_feats) {
        const int fi = f - 1;  // 0-based
        if (fi >= 0 && fi < maxj)
          feat_to_trees[fi].push_back(tt);
      }
    }
  }

  // feature f (1-based)를 split에 사용하는 tree index 목록
  // [Fix-2] fi >= maxj인 경우 빈 벡터 반환 (방어적)
  const std::vector<int>& trees_for(int f) const {
    static const std::vector<int> empty;
    const int fi = f - 1;
    if (fi < 0 || fi >= maxj) return empty;
    return feat_to_trees[fi];
  }
};


// =============================================================================
// §3  TreesCache: XPtr 캐시 컨테이너
// =============================================================================

struct TreesCache {
  std::vector<TreeDP> trees;
  FeatIndex           fidx;    // [II-3]
  int                 maxj;    // [Fix-2] 전체 feature 최대 index (1-based)
  int                 max_nn;  // [Perf-4] 전체 트리 최대 노드 수 — 1회 계산 캐싱
};


// =============================================================================
// §4  XPtr 캐싱: prepare_trees_xptr / get_cache_ref / get_trees_ref
// =============================================================================

// [[Rcpp::export]]
SEXP prepare_trees_xptr(List trees_list) {
  const int n = trees_list.size();
  auto* cache = new TreesCache();
  cache->trees.resize(n);

  int maxj   = 0;
  int max_nn = 0;  // [Perf-4] 1회 계산
  for (int tt = 0; tt < n; ++tt) {
    cache->trees[tt] = list_to_tree(as<List>(trees_list[tt]));
    for (int f : cache->trees[tt].used_feats)
      if (f > maxj) maxj = f;
    const int nn = (int)cache->trees[tt].fcol_local.size();
    if (nn > max_nn) max_nn = nn;
  }

  cache->maxj   = maxj;
  cache->max_nn = max_nn;  // [Perf-4]
  // [Fix-2] fidx에 실제 최대 feature index를 명시적으로 전달
  cache->fidx.build(cache->trees, maxj);

  return Rcpp::XPtr<TreesCache>(cache, true);
}

static inline const TreesCache& get_cache_ref(SEXP xptr) {
  return *Rcpp::XPtr<TreesCache>(xptr);
}
static inline const std::vector<TreeDP>& get_trees_ref(SEXP xptr) {
  return get_cache_ref(xptr).trees;
}


// =============================================================================
// §5  핵심 DP 평가 함수
//
//  (A) eval_one_mask  — uint64_t bitmask 기반 (inner exact 청크 전용)
//  (B) eval_one_flag  — uint8_t flag 배열 기반 (MC / enum 전용, 전체 forest)
//  (C) eval_one_tree_flag — 단일 tree (역인덱스 skip용)
//
// 모두 [II-1] uint8_t 적용, 수치 로직 불변.
// =============================================================================

// (A) bitmask 기반 — eval_one_mask
static inline double eval_one_mask(
    const std::vector<TreeDP>& trees,
    const std::vector<double>& x_row,
    const std::vector<uint8_t>& outer_known,
    const std::vector<int>&    group_bitpos,
    uint64_t                   mask,
    std::vector<double>&       Vbuf,
    int                        maxj
) {
  double total = 0.0;
  for (const auto& tr : trees) {
    const int nn = (int)tr.fcol_local.size();
    if (nn <= 0) continue;
    std::fill(Vbuf.begin(), Vbuf.begin() + nn + 1, 0.0);

    for (int li = 0; li < (int)tr.post_local.size(); ++li) {
      const int loc = tr.post_local[li];
      if (loc < 1 || loc > nn) continue;
      const int j  = tr.fcol_local[loc - 1];
      const int nk = tr.nkids[loc - 1];

      if (j == 0) { Vbuf[loc] = tr.pred_val[loc - 1]; continue; }
      if      (nk == 0) Vbuf[loc] = 0.0;
      else if (nk == 1) Vbuf[loc] = Vbuf[tr.kid1[loc - 1]];
      else              Vbuf[loc] = tr.p1[loc-1] * Vbuf[tr.kid1[loc-1]]
                                  + tr.p2[loc-1] * Vbuf[tr.kid2[loc-1]];
      if (j < 1 || j > maxj) continue;

      bool is_known = (outer_known[j - 1] != 0);
      if (!is_known) {
        const int bp = group_bitpos[j - 1];
        if (bp > 0) is_known = ((mask >> (uint64_t)(bp - 1)) & 1ULL) == 1ULL;
      }
      if (!is_known) continue;

      const double xv = (j <= (int)x_row.size())
                      ? x_row[j - 1]
                      : std::numeric_limits<double>::quiet_NaN();
      int det = 0;
      if (std::isnan(xv)) {
        // [Fix-miss] has_miss_dir[loc-1]=1 이면 missL이 항상 유효한 miss 방향을 가리킴
        if (tr.has_missing && loc - 1 < (int)tr.has_miss_dir.size()
            && tr.has_miss_dir[loc - 1]) det = tr.missL[loc - 1];
      } else {
        det = (xv <= tr.split_val[loc - 1]) ? tr.yesL[loc - 1] : tr.noL[loc - 1];
      }
      if (det > 0) Vbuf[loc] = Vbuf[det];
    }
    total += Vbuf[tr.root_local];
  }
  return total;
}


// (B) flag 배열 기반 — eval_one_flag (전체 forest)
static inline double eval_one_flag(
    const std::vector<TreeDP>& trees,
    const std::vector<double>& x_row,
    const std::vector<uint8_t>& known_flag,
    std::vector<double>&       Vbuf,
    int                        maxj
) {
  double total = 0.0;
  for (const auto& tr : trees) {
    const int nn = (int)tr.fcol_local.size();
    if (nn <= 0) continue;
    std::fill(Vbuf.begin(), Vbuf.begin() + nn + 1, 0.0);

    for (int li = 0; li < (int)tr.post_local.size(); ++li) {
      const int loc = tr.post_local[li];
      if (loc < 1 || loc > nn) continue;
      const int j  = tr.fcol_local[loc - 1];
      const int nk = tr.nkids[loc - 1];

      if (j == 0) { Vbuf[loc] = tr.pred_val[loc - 1]; continue; }
      if      (nk == 0) Vbuf[loc] = 0.0;
      else if (nk == 1) Vbuf[loc] = Vbuf[tr.kid1[loc - 1]];
      else              Vbuf[loc] = tr.p1[loc-1] * Vbuf[tr.kid1[loc-1]]
                                  + tr.p2[loc-1] * Vbuf[tr.kid2[loc-1]];
      if (j < 1 || j > maxj) continue;
      if (!known_flag[j - 1]) continue;

      const double xv = (j <= (int)x_row.size())
                      ? x_row[j - 1]
                      : std::numeric_limits<double>::quiet_NaN();
      int det = 0;
      if (std::isnan(xv)) {
        // [Fix-miss] has_miss_dir[loc-1]=1 이면 missL이 항상 유효한 miss 방향을 가리킴
        if (tr.has_missing && loc - 1 < (int)tr.has_miss_dir.size()
            && tr.has_miss_dir[loc - 1]) det = tr.missL[loc - 1];
      } else {
        det = (xv <= tr.split_val[loc - 1]) ? tr.yesL[loc - 1] : tr.noL[loc - 1];
      }
      if (det > 0) Vbuf[loc] = Vbuf[det];
    }
    total += Vbuf[tr.root_local];
  }
  return total;
}


// (C) 단일 tree — eval_one_tree_flag (역인덱스 skip용)
static inline double eval_one_tree_flag(
    const TreeDP&               tr,
    const std::vector<double>&  x_row,
    const std::vector<uint8_t>& known_flag,
    std::vector<double>&        Vbuf,
    int                         maxj
) {
  const int nn = (int)tr.fcol_local.size();
  if (nn <= 0) return 0.0;
  std::fill(Vbuf.begin(), Vbuf.begin() + nn + 1, 0.0);

  for (int li = 0; li < (int)tr.post_local.size(); ++li) {
    const int loc = tr.post_local[li];
    if (loc < 1 || loc > nn) continue;
    const int j  = tr.fcol_local[loc - 1];
    const int nk = tr.nkids[loc - 1];

    if (j == 0) { Vbuf[loc] = tr.pred_val[loc - 1]; continue; }
    if      (nk == 0) Vbuf[loc] = 0.0;
    else if (nk == 1) Vbuf[loc] = Vbuf[tr.kid1[loc - 1]];
    else              Vbuf[loc] = tr.p1[loc-1] * Vbuf[tr.kid1[loc-1]]
                                + tr.p2[loc-1] * Vbuf[tr.kid2[loc-1]];
    if (j < 1 || j > maxj) continue;
    if (!known_flag[j - 1]) continue;

    const double xv = (j <= (int)x_row.size())
                    ? x_row[j - 1]
                    : std::numeric_limits<double>::quiet_NaN();
    int det = 0;
    if (std::isnan(xv)) {
      // [Fix-miss] has_miss_dir[loc-1]=1 이면 missL이 항상 유효한 miss 방향을 가리킴
      if (tr.has_missing && loc - 1 < (int)tr.has_miss_dir.size()
          && tr.has_miss_dir[loc - 1]) det = tr.missL[loc - 1];
    } else {
      det = (xv <= tr.split_val[loc - 1]) ? tr.yesL[loc - 1] : tr.noL[loc - 1];
    }
    if (det > 0) Vbuf[loc] = Vbuf[det];
  }
  return Vbuf[tr.root_local];
}


// =============================================================================
// §6  SECTION A — 유틸리티 export
// =============================================================================

// [[Rcpp::export]]
IntegerVector popcount_table_cpp(int max_mask) {
  if (max_mask < 0) return IntegerVector(0);
  IntegerVector pc(max_mask + 1, 0);
  for (int mm = 1; mm <= max_mask; ++mm) pc[mm] = pc[mm >> 1] + (mm & 1);
  return pc;
}

// [[Rcpp::export]]
NumericVector shapley_weights_cpp(int m) {
  if (m <= 0) return NumericVector(0);
  NumericVector w(m);
  for (int s = 0; s < m; ++s) w[s] = safe_shapley_w(s, m);
  return w;
}


// =============================================================================
// §7  SECTION B — known_cols 기반 단일 평가 (eval_v_knownset_cpp)
//
// [Fix-C03] eff_maxj = max(R 전달 maxj, cache_ref.maxj) 로 통일
//   - R의 maxj = dp_model$used_feat_cols의 max (실제 사용 feature 기준)
//   - cache의 maxj = prepare_trees_xptr에서 계산 (동일하거나 더 클 수 있음)
//   - 두 값 중 큰 쪽을 사용하여 known_flag 범위 보장
// =============================================================================

// [[Rcpp::export]]
double eval_v_knownset_cpp(
    SEXP          trees_xptr,
    NumericVector x_row,
    IntegerVector known_cols,  // 1-based
    int           maxj
) {
  const TreesCache& cache_ref = get_cache_ref(trees_xptr);
  const std::vector<TreeDP>& trees = cache_ref.trees;

  // [Fix-C03] eff_maxj: R 전달값과 cache 값 중 큰 쪽 사용
  const int eff_maxj = std::max(maxj, cache_ref.maxj);
  if (eff_maxj <= 0) return 0.0;

  // [Perf-4] cache에서 직접 참조 — 루프 계산 제거
  std::vector<double>  Vbuf(cache_ref.max_nn + 1, 0.0);
  std::vector<uint8_t> known_flag(eff_maxj, 0);

  std::vector<int> set_bits;
  set_bits.reserve(known_cols.size());
  for (int i = 0; i < (int)known_cols.size(); ++i) {
    const int c = known_cols[i] - 1;
    if (c >= 0 && c < eff_maxj && !known_flag[c]) {
      known_flag[c] = 1;
      set_bits.push_back(c);
    }
  }

  std::vector<double> xvec = as<std::vector<double>>(x_row);
  const double result = eval_one_flag(trees, xvec, known_flag, Vbuf, eff_maxj);

  for (int b : set_bits) known_flag[b] = 0;
  return result;
}


// =============================================================================
// §8  SECTION C — treeowen_exact_enum C++ 가속
//     (inner_shapley_enum_one_group_cpp)
//
// [Revert-C07] FeatIndex 역인덱스 제거 — v0/v1 모두 eval_one_flag full scan으로 복원
//   - v5.10의 relevant_T, relevant_work, base_T, work_cache 구조 제거
//   - build_relevant_mask() 헬퍼 함수 제거 (이 함수에서만 사용)
//   - fidx 참조 제거 (inner_shapley_exact_cpp / inner_mc_perm_cpp의 fidx는 유지)
//   - 수치 로직 불변: v0 = eval(T∪S), v1 = eval(T∪S∪{feat_i}) 동일
//
// 알고리즘:
//   group g에 대해 outer T (2^(K-1) 부분집합) × inner S (2^(mG-1) 부분집합) 이중 루프.
//   각 (T, S, feat_i) 조합에서 v0 = v(T∪S), v1 = v(T∪S∪{feat_i})를
//   eval_one_flag (전체 B trees full scan)으로 계산.
//   combined_flag를 flip-back 패턴으로 재사용하여 메모리 재할당 없음.
// =============================================================================

// [Revert-C07] build_relevant_mask 제거 — inner_shapley_enum_one_group_cpp만
// 사용하던 헬퍼였으므로 함께 삭제. inner_shapley_exact_cpp / inner_mc_perm_cpp
// 의 역인덱스(fidx)는 그대로 유지.

// [[Rcpp::export]]
NumericVector inner_shapley_enum_one_group_cpp(
    SEXP          trees_xptr,
    NumericVector x_row,
    int           maxj,
    IntegerVector group_feats_g,
    List          group_feats_others,
    int           K,
    int           p
) {
  NumericVector out(p, 0.0);
  const int mG = group_feats_g.size();
  if (mG == 0) return out;

  const TreesCache&          cache_ref = get_cache_ref(trees_xptr);
  const std::vector<TreeDP>& trees     = cache_ref.trees;
  // [Revert-C07] fidx 참조 제거 — 역인덱스 미사용
  const int B = (int)trees.size(); (void)B;

  const int eff_maxj = std::max(maxj, cache_ref.maxj);
  std::vector<double> Vbuf(cache_ref.max_nn + 1, 0.0);
  std::vector<double> xvec = as<std::vector<double>>(x_row);

  const int d_minus = (int)group_feats_others.size();
  std::vector<std::vector<int>> others(d_minus);
  for (int gg = 0; gg < d_minus; ++gg)
    others[gg] = as<std::vector<int>>(group_feats_others[gg]);

  std::vector<int> gf = as<std::vector<int>>(group_feats_g);

  // [I-8] Gminus_all: outer T 루프 밖 1회 precompute
  std::vector<std::vector<int>> Gminus_all(mG);
  for (int ii = 0; ii < mG; ++ii) {
    Gminus_all[ii].reserve(mG - 1);
    for (int jj = 0; jj < mG; ++jj)
      if (jj != ii) Gminus_all[ii].push_back(gf[jj]);
  }

  // [II-1] uint8_t flip-back 배열
  std::vector<uint8_t> T_feat_flag(eff_maxj, 0);
  std::vector<uint8_t> combined_flag(eff_maxj, 0);

  // [Fix-8] d_minus >= 30 overflow 방어
  if (d_minus >= 30)
    Rcpp::stop("inner_shapley_enum_one_group_cpp: d_minus=%d >= 30 (K too large for exact enum).",
               d_minus);
  const int64_t n_outer = INT64_C(1) << d_minus;

  for (int64_t tmask = 0; tmask < n_outer; ++tmask) {
    const int    t_size = popcount64((uint64_t)tmask);
    const double wT     = safe_shapley_w(t_size, K);

    // [Fix-5] T_bits reserve: 실제 세트 비트의 others 크기 합산
    std::vector<int> T_bits;
    {
      int total_feat = 0;
      for (int gg = 0; gg < d_minus; ++gg)
        if ((tmask >> gg) & 1) total_feat += (int)others[gg].size();
      T_bits.reserve(total_feat);
    }
    for (int gg = 0; gg < d_minus; ++gg) {
      if ((tmask >> gg) & 1) {
        for (int f : others[gg]) {
          const int fi = f - 1;
          if (fi >= 0 && fi < eff_maxj && !T_feat_flag[fi]) {
            T_feat_flag[fi] = 1;
            T_bits.push_back(fi);
          }
        }
      }
    }

    // T_bits를 combined_flag에 반영
    for (int b : T_bits) combined_flag[b] = 1;

    for (int ii = 0; ii < mG; ++ii) {
      const int feat_i    = gf[ii];
      const int feat_i_fi = feat_i - 1;
      const std::vector<int>& Gminus = Gminus_all[ii];
      const int mG_minus = (int)Gminus.size();
      const int n_inner  = 1 << mG_minus;

      for (int smask = 0; smask < n_inner; ++smask) {
        const int    s_size = popcount64((uint64_t)smask);
        const double wS     = safe_shapley_w(s_size, mG);

        // S_extra: smask에 해당하는 Gminus 내 추가 feature를 combined_flag에 반영
        std::vector<int> S_extra;
        S_extra.reserve(mG_minus);
        for (int jj = 0; jj < mG_minus; ++jj) {
          if ((smask >> jj) & 1) {
            const int fi = Gminus[jj] - 1;
            if (fi >= 0 && fi < eff_maxj && !combined_flag[fi]) {
              combined_flag[fi] = 1;
              S_extra.push_back(fi);
            }
          }
        }

        // [Revert-C07] v0 계산: combined_flag = T ∪ S → 전체 forest full scan
        const double v0 = eval_one_flag(trees, xvec, combined_flag, Vbuf, eff_maxj);

        // [Revert-C07] v1 계산: combined_flag = T ∪ S ∪ {feat_i} → 전체 forest full scan
        bool feat_i_added = false;
        if (feat_i_fi >= 0 && feat_i_fi < eff_maxj && !combined_flag[feat_i_fi]) {
          combined_flag[feat_i_fi] = 1;
          feat_i_added = true;
        }
        const double v1 = eval_one_flag(trees, xvec, combined_flag, Vbuf, eff_maxj);

        if (feat_i_added)     combined_flag[feat_i_fi] = 0;
        for (int b : S_extra) combined_flag[b]         = 0;

        if (feat_i_fi >= 0 && feat_i_fi < p)
          out[feat_i_fi] += wT * wS * (v1 - v0);
      }
    }

    // T_bits 제거 (combined_flag, T_feat_flag flip-back)
    for (int b : T_bits) {
      combined_flag[b] = 0;
      T_feat_flag[b]   = 0;
    }

    if ((tmask & INT64_C(255)) == INT64_C(255)) Rcpp::checkUserInterrupt();
  }

  return out;
}


// =============================================================================
// §9  (제거됨: eval_v_forest_chunk_exact_cpp)
//     [Fix-C08] R에서 실제로 호출되지 않는 dead export였으므로 v5.10에서 제거.
//     R fallback 경로는 직접 .eval_v_masks_forest_chunked_prepared()를 사용.
// =============================================================================


// =============================================================================
// §10  SECTION E — inner Shapley exact (inner_shapley_exact_cpp)
//      m <= 30: Gray code + 역인덱스  [II-4]
//      m 31-60: 청크 bit-insert trick
//
// [Fix-C01] m<=30 Gray code 경로: delta 방식 합산으로 O(B) → O(|affected|)
//   - prev_total 변수 도입: 이전 Gray code step의 합계를 보관
//   - affected tree만 재계산하고 그 delta를 prev_total에 더함
//   - tree_root_cache 갱신으로 이후 step에서도 delta 계산 가능
//   - 수치 결과: delta 합산 == 전체 합산 (동치 증명 가능)
//   - 성능: Ranger 500 trees → step당 ~500 덧셈 → ~|affected| 덧셈
//
// [Fix-C06] eff_maxj 범위 검사 통일
// =============================================================================

// [[Rcpp::export]]
NumericVector inner_shapley_exact_cpp(
    SEXP          trees_xptr,
    NumericVector x_row,
    IntegerVector feats_g,
    IntegerVector outer_known_cols,
    int           maxj,
    int           p,
    int           chunk_size
) {
  const int m = feats_g.size();
  NumericVector out(p, 0.0);
  if (m == 0) return out;
  if (m > 60) Rcpp::stop("inner_shapley_exact_cpp: m > 60 not supported");

  const TreesCache&          cache_ref = get_cache_ref(trees_xptr);
  const std::vector<TreeDP>& trees     = cache_ref.trees;
  const FeatIndex&           fidx      = cache_ref.fidx;   // [II-4]

  // [Fix-C03/C06] eff_maxj 통일
  const int eff_maxj = std::max(maxj, cache_ref.maxj);
  if (eff_maxj <= 0) return out;

  // [II-1] outer known flag (uint8_t)
  std::vector<uint8_t> ok_vec(eff_maxj, 0);
  for (int i = 0; i < (int)outer_known_cols.size(); ++i) {
    const int c = outer_known_cols[i] - 1;
    if (c >= 0 && c < eff_maxj) ok_vec[c] = 1;
  }

  // group_bitpos (m>30 청크 경로에서 eval_one_mask에 전달)
  std::vector<int> group_bitpos(eff_maxj, 0);
  std::vector<int> fg = as<std::vector<int>>(feats_g);
  for (int ii = 0; ii < m; ++ii) {
    const int f = fg[ii] - 1;
    if (f >= 0 && f < eff_maxj) group_bitpos[f] = ii + 1;
  }

  std::vector<double> xvec = as<std::vector<double>>(x_row);

  // [Perf-4] cache에서 직접 참조
  std::vector<double> Vbuf(cache_ref.max_nn + 1, 0.0);

  // Shapley weight 테이블 (inner group 크기 m 기준)
  std::vector<double> w(m);
  for (int s = 0; s < m; ++s) w[s] = safe_shapley_w(s, m);

  // ── m <= 30: Gray code + 역인덱스 + [Fix-C01] delta 합산 ──────────────────
  //
  // 알고리즘:
  //   i = 0..M-1 에서 cur_gray = i ^ (i>>1) 순서로 순회.
  //   diff = prev_gray ^ cur_gray 는 항상 1-bit.
  //   kf(bitmask) == cur_gray 가 항상 성립.
  //
  // [Fix-C01] delta 합산:
  //   prev_total: 이전 step의 forest 합계 (초기값 = vvals[0])
  //   affected tree만 재계산 → 각 tree의 delta를 prev_total에 더함
  //   tree_root_cache도 갱신 → 다음 step에서 delta 계산에 재사용
  if (m <= 30) {
    const int M = 1 << m;
    const int B = (int)trees.size();

    // base 상태 (inner feature 없음) tree root value 캐시
    std::vector<double> tree_root_cache(B, 0.0);
    for (int tt = 0; tt < B; ++tt)
      tree_root_cache[tt] = eval_one_tree_flag(trees[tt], xvec, ok_vec, Vbuf, eff_maxj);

    std::vector<double> vvals(M, 0.0);
    {
      // [Fix-C01] prev_total 초기화: vvals[0] = 전체 합계
      double s0 = 0.0;
      for (int tt = 0; tt < B; ++tt) s0 += tree_root_cache[tt];
      vvals[0] = s0;
    }

    // kf: ok_vec 복사본 (inner feature set/reset용)
    std::vector<uint8_t> kf(ok_vec.begin(), ok_vec.end());
    kf.resize(eff_maxj, 0);  // 혹시 ok_vec가 짧은 경우 방어

    uint64_t prev_gray  = 0ULL;
    double   prev_total = vvals[0];  // [Fix-C01] delta 합산을 위한 이전 합계 보관

    for (int i = 1; i < M; ++i) {
      const uint64_t cur_gray = (uint64_t)i ^ ((uint64_t)i >> 1);
      const uint64_t diff     = prev_gray ^ cur_gray;
      const int      bit_pos  = portable_ctz64(diff);
      const bool     adding   = ((cur_gray >> bit_pos) & 1ULL) == 1ULL;

      const int feat_1based = fg[bit_pos];   // 1-based
      const int feat_0based = feat_1based - 1;

      if (feat_0based >= 0 && feat_0based < eff_maxj)
        kf[feat_0based] = adding ? 1 : 0;

      // [Fix-C01] delta 방식: affected tree만 재계산 후 delta를 prev_total에 더함
      // 기존: 전체 B개 tree_root_cache 합산 → O(B) per step
      // 수정: affected tree의 (new_val - old_val) delta만 누적 → O(|affected|) per step
      const std::vector<int>& affected = fidx.trees_for(feat_1based);
      double total = prev_total;  // 이전 합계에서 시작
      for (int tt : affected) {
        const double new_val = eval_one_tree_flag(trees[tt], xvec, kf, Vbuf, eff_maxj);
        total                += (new_val - tree_root_cache[tt]);
        tree_root_cache[tt]  = new_val;  // 캐시 갱신 (다음 step에서 delta 계산에 사용)
      }

      vvals[cur_gray] = total;
      prev_total      = total;   // [Fix-C01] 다음 step을 위해 갱신
      // 정확성 근거: kf의 inner feature bitmask == cur_gray 항상 성립
      // → vvals[cur_gray] = v(subset=cur_gray) 로 mask 인덱스와 일치

      prev_gray = cur_gray;

      // [Fix-3] 65536 step마다 R 중단 신호 확인
      if ((i & 0xFFFF) == 0) Rcpp::checkUserInterrupt();
    }

    // kf 복원 (ok_vec 상태로)
    for (int ii = 0; ii < m; ++ii) {
      const int f0 = fg[ii] - 1;
      if (f0 >= 0 && f0 < eff_maxj) kf[f0] = ok_vec[f0];
    }

    // Shapley weight 적용 (vvals는 mask 인덱스 기준)
    std::vector<int> pc(M, 0);
    for (int mm = 1; mm < M; ++mm) pc[mm] = pc[mm >> 1] + (mm & 1);

    for (int ii = 0; ii < m; ++ii) {
      const int bit_i = 1 << ii;
      double shapley_acc = 0.0;
      for (int S = 0; S < M; ++S)
        if ((S & bit_i) == 0)
          shapley_acc += w[pc[S]] * (vvals[S | bit_i] - vvals[S]);
      const int f_out = fg[ii] - 1;
      if (f_out >= 0 && f_out < p) out[f_out] = shapley_acc;
    }

    return out;
  }

  // ── m 31..60: 청크 bit-insert trick ──────────────────────────────────────
  // group_bitpos는 여기서 eval_one_mask에 전달됨
  if (chunk_size <= 0) chunk_size = 131072;

  for (int ii = 0; ii < m; ++ii) {
    const uint64_t bit_i = 1ULL << (uint64_t)ii;
    const uint64_t n_sub = (m == 1) ? 1ULL : (1ULL << (uint64_t)(m - 1));
    double feat_acc = 0.0;

    for (uint64_t sub = 0; sub < n_sub; ) {
      const uint64_t sub_end = std::min(sub + (uint64_t)chunk_size, n_sub);
      for (uint64_t s_idx = sub; s_idx < sub_end; ++s_idx) {
        const uint64_t lo   = s_idx & (bit_i - 1ULL);
        const uint64_t hi   = (s_idx & ~(bit_i - 1ULL)) << 1;
        const uint64_t S_wo = lo | hi;
        const uint64_t S_w  = S_wo | bit_i;
        const double v0 = eval_one_mask(trees, xvec, ok_vec, group_bitpos, S_wo, Vbuf, eff_maxj);
        const double v1 = eval_one_mask(trees, xvec, ok_vec, group_bitpos, S_w,  Vbuf, eff_maxj);
        feat_acc += w[popcount64(S_wo)] * (v1 - v0);
      }
      sub = sub_end;
      Rcpp::checkUserInterrupt();
    }

    const int f_out = fg[ii] - 1;
    if (f_out >= 0 && f_out < p) out[f_out] = feat_acc;
  }

  return out;
}


// =============================================================================
// §11  SECTION F — inner MC 순열 샘플링 (inner_mc_perm_cpp)
//
// [Fix-C02] contrib_abs.reserve 부족 방어: reserve(n_perm*4+8)
//   - adaptive=TRUE 시 n_perm 이상 실행 가능 (max_mc까지)
//   - R 측에서 max_mc를 C++에 전달하지 않으므로 n_perm*4로 보수적 예약
// [II-1] known_flag: vector<uint8_t>
// [II-3] feature→tree 역인덱스 기반 skip
// [Fix-4] od 루프: base_sum 1회 계산, perm_cache std::copy reset
// =============================================================================

// [[Rcpp::export]]
NumericVector inner_mc_perm_cpp(
    SEXP          trees_xptr,
    NumericVector x_row,
    IntegerVector feats_g,
    IntegerVector outer_known_cols,
    int           maxj,
    int           p,
    int           n_perm,
    bool          antithetic,
    double        target_se  = -1.0,  // [Perf-3] adaptive SE 조기 종료. -1 = 비활성
    int           min_perm   = 32     // [Perf-3] adaptive 최소 order 수
) {
  const int m = feats_g.size();
  NumericVector acc(p, 0.0);
  if (m == 0) return acc;
  if (n_perm <= 0) n_perm = 1;

  const TreesCache&          cache_ref = get_cache_ref(trees_xptr);
  const std::vector<TreeDP>& trees     = cache_ref.trees;
  const FeatIndex&           fidx      = cache_ref.fidx;   // [II-3]
  const int B = (int)trees.size();

  // [Fix-C03] eff_maxj 통일
  const int eff_maxj = std::max(maxj, cache_ref.maxj);
  if (eff_maxj <= 0) return acc;

  // [Perf-4] cache에서 직접 참조
  std::vector<double>  Vbuf(cache_ref.max_nn + 1, 0.0);
  std::vector<double>  xvec = as<std::vector<double>>(x_row);
  std::vector<int>     fg   = as<std::vector<int>>(feats_g);

  // [II-1] outer known flag (uint8_t, flip-back)
  std::vector<uint8_t> known_flag(eff_maxj, 0);
  std::vector<int>     base_bits;
  base_bits.reserve(outer_known_cols.size());
  for (int i = 0; i < (int)outer_known_cols.size(); ++i) {
    const int c = outer_known_cols[i] - 1;
    if (c >= 0 && c < eff_maxj && !known_flag[c]) {
      known_flag[c] = 1;
      base_bits.push_back(c);
    }
  }

  // [II-3] base 상태에서 각 tree root value 캐시
  std::vector<double> tree_root_cache(B, 0.0);
  for (int tt = 0; tt < B; ++tt)
    tree_root_cache[tt] = eval_one_tree_flag(trees[tt], xvec, known_flag, Vbuf, eff_maxj);

  // [Fix-4] base_sum: od 루프 밖에서 1회 계산
  double base_sum = 0.0;
  for (int tt = 0; tt < B; ++tt) base_sum += tree_root_cache[tt];

  std::vector<int>    perm(m);
  std::iota(perm.begin(), perm.end(), 0);

  // [Fix-4] perm_cache: 루프 밖 1회 alloc, od마다 std::copy로 reset
  std::vector<double> perm_cache(B);

  const int n_orders = antithetic ? 2 : 1;

  // [Perf-3] adaptive SE: draw마다 |기여| 합 추적
  const bool do_adaptive = (target_se > 0.0) && (min_perm >= 1);
  // [Fix-C02] reserve 부족 방어: n_perm*4+8 (max_mc는 알 수 없으므로 보수적으로)
  std::vector<double> contrib_abs;
  if (do_adaptive) contrib_abs.reserve(n_perm * 4 + 8);
  int total_orders = 0;

  GetRNGstate();
  for (int draw = 0; draw < n_perm; ++draw) {
    // Fisher-Yates shuffle
    for (int i = m - 1; i > 0; --i) {
      int j = (int)(unif_rand() * (i + 1));
      if (j > i) j = i;
      std::swap(perm[i], perm[j]);
    }

    for (int od = 0; od < n_orders; ++od) {
      // [Fix-4] perm_cache reset (재할당 없이 std::copy)
      std::copy(tree_root_cache.begin(), tree_root_cache.end(), perm_cache.begin());

      std::vector<int> added_bits;
      added_bits.reserve(m);

      // [Fix-4] v_curr: base_sum 재사용
      double v_curr = base_sum;
      double abs_sum = 0.0;  // [Perf-3] 이번 order의 |기여| 합

      for (int step = 0; step < m; ++step) {
        const int ii = (od == 0) ? perm[step] : perm[m - 1 - step];
        const int f1 = fg[ii];   // 1-based
        const int f0 = f1 - 1;  // 0-based

        if (f0 >= 0 && f0 < eff_maxj && !known_flag[f0]) {
          known_flag[f0] = 1;
          added_bits.push_back(f0);

          // [II-3] f1을 사용하는 tree만 재계산, Δv를 캐시 차분으로 계산
          const std::vector<int>& affected = fidx.trees_for(f1);
          double v_new = v_curr;
          for (int tt : affected) {
            const double new_val = eval_one_tree_flag(trees[tt], xvec, known_flag, Vbuf, eff_maxj);
            v_new        += (new_val - perm_cache[tt]);
            perm_cache[tt] = new_val;
          }

          const double delta = v_new - v_curr;
          if (f0 < p) acc[f0] += delta;
          abs_sum += std::fabs(delta);
          v_curr = v_new;
        }
      }

      // flip-back: known_flag 복원
      for (int b : added_bits) known_flag[b] = 0;

      total_orders++;
      if (do_adaptive) contrib_abs.push_back(abs_sum);
    }

    // [Perf-3] adaptive 조기 종료: min_perm order 이후 매 draw마다 SE 체크
    // total_orders 기준 — antithetic 시 draw당 2 orders (R fallback n_stat 단위와 일치)
    if (do_adaptive && total_orders >= min_perm) {
      const int n_s = (int)contrib_abs.size();
      if (n_s >= 2) {
        double mean_s = 0.0;
        for (double v : contrib_abs) mean_s += v;
        mean_s /= n_s;
        double var_s = 0.0;
        for (double v : contrib_abs) var_s += (v - mean_s) * (v - mean_s);
        var_s /= (n_s - 1);
        const double se = std::sqrt(var_s / n_s);
        if (se <= target_se) break;
      }
    }

    if ((draw & 63) == 0) Rcpp::checkUserInterrupt();
  }
  PutRNGstate();

  // known_flag cleanup (base outer bits)
  for (int b : base_bits) known_flag[b] = 0;

  // [Perf-3] acc를 실제 처리된 order 수로 정규화
  // total_orders >= 1 이면 항상 나눔
  if (total_orders >= 1) {
    const double inv = 1.0 / total_orders;
    for (int i = 0; i < p; ++i) acc[i] *= inv;
  }

  return acc;
}
