#include "iip.h"
/*-------------------------------------------------------------------*/
/*   ヘッセ行列の固有値からコーナー点を検出するサンプルプログラム    */
/*                       by H.Nagahashi                              */
/*      第一引数：PGM画像オブジェクト                                */
/*      第二引数：固有値（大）配列オブジェクト                       */
/*      第三引数：固有値（小）配列オブジェクト                       */
/*-------------------------------------------------------------------*/
inline bool is_edge(double val0, double val1, double thre1, double thre2)
{
    return (val0 < thre1) && (val0 >= thre2) && (val1 >= thre2);
}

void extractHessianCorner(ARY2I& img, ARY2D& ary0, ARY2D& ary1)
{
  int height    = img.height,  width = img.width;
  int2D pic     = img.array;
  double2D eig0 = ary0.array,   eig1 = ary1.array;

  //局所最適解を探索するための２次元配列
  int2D tmp = allocInt2D(height, width);

  //大きい固有値の最大、最少を検出
  double mn=10000.0, mx = -10000.0;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      if(eig0[i][j] > mx) mx = eig0[i][j];
      if(eig0[i][j] < mn) mn = eig0[i][j];
    }
  }
  // 大きい固有値の最大値の60％(thresh1)以上を持つ点をエッジ成分として却下
  // 大きい固有値の最大値の5％(thresh2)を共通閾値として決定
  // 大きい固有値と小さい固有値が同時にthresh2以上である画素を角点として検出
  //
#if 0
  double thresh1 = mx * 0.6,  thresh2 = mx * 0.05;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      tmp[i][j] = 0;
      if(eig0[i][j] >= thresh1) continue;
      if(eig0[i][j] >= thresh2 && eig1[i][j] >= thresh2){
        tmp[i][j] = 255;
      }
    }
  }
  //条件を満たす点のうち，周囲にも条件を満たす点が存在するものエッジ成分とする
  //基本条件は上と同様

#elseif 0
  double thresh1 = mx * 0.6,  thresh2 = mx * 0.05;
  int window = 3;
  bool exist;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      tmp[i][j] = 0;
      exist = false;
      if (!(is_edge(eig0[i][j], eig1[i][j], thresh1, thresh2)))
        continue;
      for (int ii = -window; ii < window; ++ii){
        for (int jj = -window; jj < window; ++jj){
          int i_ = std::min(std::max(0, ii + i), height);
          int j_ = std::min(std::max(0, jj + j), width);
          if (is_edge(eig0[i_][j_], eig1[i_][j_], thresh1, thresh2)){
            exist = true;
            break;
          }
        }
        if (exist) break;
      }
      if (exist)
        tmp[i][j] = 255;
    }
  }
#else
  double thresh1 = mx * 0.6,  thresh2 = mx * 0.05;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      tmp[i][j] = 0;
      if(eig0[i][j] >= thresh1) continue;
      if(eig0[i][j] + eig1[i][j] >= thresh2 * 3){
        tmp[i][j] = 255;
      }
    }
  }
#endif
  // コーナー点候補領域（複数の点の集合）の点に対し、各点の大きい固有値と
  // 小さい固有値の差が最も少ない点を局所最適点として求める。（他の最適化もあり）
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      if(!tmp[i][j]) continue;
      int i0=i, j0=j, flag = 1;
      while(flag) {
         flag = 0;
         double min_dif = eig0[i0][j0] - eig1[i0][j0];
         for(int k = -5; k <= 5; k++){  //前後５,左右５の領域を探索
           for(int m = -5; m <= 5; m++){
             if(k == 0 && m == 0) continue;
	     if(i0+k < 0 || i0+k >= height) continue;
	     if(j0+m < 0 || j0+m >= width) continue;
             if(!tmp[i0+k][j0+m])continue;

	     if( (eig0[i0+k][j0+m] - eig1[i0+k][j0+m]) < min_dif){
	          min_dif = eig0[i0+k][j0+m] - eig1[i0+k][j0+m];
	          i0 = i0+k; j0 = j0+m;
	          flag = 1;
	      }
	      else tmp[i0+k][j0+m] = 0;
           }
	  }
      }
      pic[i0][j0] = 255;
    }
  }
  //作業領域の解放
  freeInt2D(tmp, height, width);
}
