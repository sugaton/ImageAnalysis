/*-------------------------------------------------------------------*/
/*            ヘッセ行列とその固有値計算サンプルプログラム           */
/*                            by H.Nagahashi                         */
/*   引数：第一引数：入力画像ARY2I型オブジェクト                     */
/*         第二引数：固有値（大）格納用実数型配列オブジェクト        */
/*         第三引数：固有値（小）格納用実数型配列オブジェクト        */
/*         第四引数：ガウス核オブジェクト                            */
/*   ガウス核のスケールパラメータは1.0～3.5程度を選択。              */
/*-------------------------------------------------------------------*/
#include "iip.h"
/*-------------------------------------------------------------------*/
/*                Calculate EigenValue of Hessian                    */
/*-------------------------------------------------------------------*/
void calcHessianEigenVal(ARY2I& img, ARY2D& ary0, ARY2D& ary1, GSKNL0& gk)
{
  int width     = img.width,   height = img.height;
  int2D pic     = img.array;
  double2D eig0 = ary0.array,    eig1   = ary1.array;

  //作業用実数型2次元配列の確保
  double2D gx2  = allocDouble2D(height,width);
  double2D gx1  = allocDouble2D(height,width);
  double2D gx0  = allocDouble2D(height,width);

  int sX = gk.step[0], sY = gk.step[1]; //各方向のマスク中間位置
  ARY2D wx0   = gk.weight[0];  //ｘ方向０次、一次、二次導関数用線形フィルタ配列
  double2D wx = wx0.array;
  ARY2D wy0   = gk.weight[1];  //ｙ方向０次、一次、二次導関数用線形フィルタ配列
  double2D wy = wy0.array;

  double sum1,sum2,sum3;
  for(int i=0 ; i < height; i++) {
    for(int j=0; j < width; j++) {
      sum1 = sum2 = sum3 = 0.0;
      for(int k = -sX; k <= sX; k++){
        if((j+k) < 0 || (j+k) >= width) continue;
        sum1 += (double)pic[i][j+k] * wx[0][sX+k];  //ｘの０次導関数畳み込み
        sum2 += (double)pic[i][j+k] * wx[1][sX+k];  //ｘの一次導関数畳み込み
        sum3 += (double)pic[i][j+k] * wx[2][sX+k];  //ｘの二次導関数畳み込み
      }
      gx0[i][j] = sum1;
      gx1[i][j] = sum2;
      gx2[i][j] = sum3;
    }
  }
  //縦方向ガウス核による縦続畳み込み処理と固有値計算
  double a, b, c;
  for(int i = 0; i < height; i++){
    for(int j = 0; j < width; j++){
      a = b = c  = 0.0;
      for(int k = -sY; k <= sY; k++){
	 if((i+k) < 0 || (i+k) >= height) continue;
	 a += gx2[i+k][j] * wy[0][sY+k];
	 b += gx0[i+k][j] * wy[2][sY+k];
         c += gx1[i+k][j] * wy[1][sY+k];
      }
      //点(i,j)のヘッセ行列 H=[{a  c}{c  b}]
      //ヘッセ行列の固有値計算 特性方程式 |H-rI|=0
      //rの２次方程式：r^2-(a+b)r+(ab-c^2)=0を解く
      double dx = (double)(a + b);
      double dy = sqrt((double)((a-b)*(a-b) + c*c) );
      double r0 = (dx + dy )/2.0;           //大きい固有値
      double r1 = (dx - dy )/2.0;           //小さい固有値
      eig0[i][j] = r0;                      //大きい固有値を格納
      eig1[i][j] = r1;                      //小さい固有値を格納
    }
  }
  //作業領域の解放
  freeDouble2D(gx0, height, width);
  freeDouble2D(gx1, height, width);
  freeDouble2D(gx2, height, width);
}
/*-------------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  //-------------課題1---------------
  double sigma, gamma;
  if (argc != 6){
    std::cout << "Usage: " << argv[0] << " input output1 output2 sigma gamma" << std::endl;
    return(1);
  }
  PGM& input = readPGM(argv[1]);
  sscanf(argv[4], "%lf", &sigma);
  sscanf(argv[5], "%lf", &gamma);
  int height = input.body.height, width = input.body.width;
  ARY2D eig0 = createARY2D(height, width);
  ARY2D eig1 = createARY2D(height, width);
  GSKNL0& gk = createGSKNL0(sigma, gamma);
  calcHessianEigenVal(input.body, eig0, eig1, gk);
  //---------------------------------
  //-------------課題2---------------
  PGM& output1 = createPGM(height, width);
  PGM& output2 = createPGM(height, width);
  double _max=0.0; //固有値の最大値
  for (int i = 0; i < height; ++i)
    for (int j = 0; j < width; ++j)
      if (_max > eig0.array[i][j]) _max = eig0.array[i][j];
  // 0~255に正規化
  for (int i = 0; i < height; ++i)
    for (int j = 0; j < width; ++j){
      output1.body.array[i][j] = (int)(255 * eig0.array[i][j] / _max);
      output2.body.array[i][j] = (int)(255 * eig1.array[i][j] / _max);
    }
  writePGM(output1, argv[2]);
  writePGM(output2, argv[3]);
  //---------------------------------
  deletePGM(output1);
  deletePGM(output2);
  deleteARY2D(eig0);
  deleteARY2D(eig1);
  deleteGSKNL0(gk);
  return 0;
}
