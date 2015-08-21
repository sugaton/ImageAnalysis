/*-------------------------------------------------------------------*/
/*              Intensity Conversion Sample Program                  */
/*                         by H.Nagahashi                            */
/*     Run:  intensityConv in_PGM hst1_PGM out_PGM hst2_PGM          */
/*-------------------------------------------------------------------*/

#include    "iip.h"
#include <cmath>

#define __M 255

/*--------------------------------------------------------------------*/
            //// visualizeHist()のサンプルプログラム ////
/*--------------------------------------------------------------------*/
void visualizeHist(ARY1I& hist, PGM& hstImage)
{
  int bins = hist.width;
  if(bins > hstImage.body.width)
    error("Size mismatch in visualizeHist",1);

  int height = hstImage.body.height, width = hstImage.body.width;
  int2D ary  = hstImage.body.array;
  int1D hst  = hist.array;
  for(int i = 0; i < height; i++)
    for(int j = 0; j < width; j++) ary[i][j] = 80; //背景濃度

  int maxv = 0;
  //最大度数の検出
  for(int i = 0; i < bins; i++) if(hst[i] > maxv ) maxv = hst[i];

  for(int b = 0; b < bins; b++){
    int rep = hst[b] * height / maxv;
    for(int i = 0; i< rep; i++) ary[height-1-i][b] = 255;
  }
}

/*--------------------------------------------------------------------*/
//                makeConvTable()のサンプルプログラム
/*--------------------------------------------------------------------*/
void makeConvTable(PGM& pgm, int1D table, int size)
{
  //単純な変換テーブルの例（階調値を1/10にする変換）
  for(int i = 0; i < size; i++) table[i] = (i/10) * 10;
}
// sqrt(M) * sqrt(f(i,j))
void makeConvTable1(PGM& pgm, int1D table, int size )
{
  for(int i = 0; i < size; i++) table[i] = std::sqrt(__M) * std::sqrt(i);
}
// f(i,j)** 2 / M
void makeConvTable2(PGM& pgm, int1D table, int size )
{
  for(int i = 0; i < size; i++) table[i] = i * i / __M;
}
// 0 (i < 50) : f(i,j) (50<= i < 200) : M (200 <= i)
void makeConvTable3(PGM& pgm, int1D table, int size )
{
  for(int i = 0; i < size; i++)
  {
    if(i < 50) table[i] = 0;
    else if(50 <= i && i < 200) table[i] = i;
    else if(200 <= i) table[i] = __M;
  }
}

void getVMAX_MIN(PGM& pgm, int *vmax, int *vmin)
{
  for(int i = 0; i < pgm.body.height; ++i){
    for (int j = 0; j < pgm.body.width; ++j){
      int val = pgm.body.array[i][j];
      if (val < *vmin) *vmin = val;
      else if (val > *vmax) *vmax = val;
    }
  }
}
void mulhist(ARY1I& hist, ARY1I& mulh)
{
  int sum = 0;
  for (int j = 0; j < hist.width; ++j){
    sum += hist.array[j];
    mulh.array[j] = sum;
  }
}
// prob 2.
void makeConvTable4(PGM& pgm, int1D table, int size, ARY1I& hist)
{
  int N = pgm.body.height * pgm.body.width;
  int vmax = 0, vmin = __M;
  getVMAX_MIN(pgm, &vmax, &vmin);
  for(int i = 0; i < size; i++) table[i] = vmin + (double)hist.array[i] * (vmax - vmin) / N;
}
/*--------------------- メインプログラム -----------------------------*/
int main(int argc, char *argv[])
{
  if (argc != 6 || (argc>1 && strcmp(argv[1],"help")==0)) {
    cout << "\nUsage: " << argv[0] << " in_pgm hst1 out_pgm hst2 type"
         << endl;
    cout << "hst1 and hst2 must be PGM type of images" << endl;
    cout << "type : 1, 2, 3, 4" << endl;
    exit(1);
  }
  int type;
  sscanf(argv[5], "%d", &type);
  /*---------------- 画像データの読み込み --------------------*/
  PGM& pgm1  = readPGM(argv[1]);
  int width = pgm1.body.width,  height = pgm1.body.height;
  int depth = pgm1.depth;

  int asize    = depth+1;
  /*----------------- ヒストグラム用配列確保 -----------------*/
  ARY1I& hist  = createARY1I(asize);
  ARY1I& mulh  = createARY1I(asize);

  PGM& pgm2    = createPGM(height, width); //出力画像用配列確保
  PGM& hstpgm1 = createPGM(256, asize);    //度数を256段階で表示
  PGM& hstpgm2 = createPGM(256, asize);    //ビン総数は(depth+1)

  int1D table  = allocInt1D(asize);        //変換表配列の確保
  histogram(pgm1.body, hist);
  // makeConvTable(pgm1, table, asize);       //輝度変換表の作成
  switch (type) {
    case 1:
      makeConvTable1(pgm1, table, asize);
      break;
    case 2:
      makeConvTable2(pgm1, table, asize);
      break;
    case 3:
      makeConvTable3(pgm1, table, asize);
      break;
    case 4:
      mulhist(hist, mulh);
      makeConvTable4(pgm1, table, asize, mulh);
      break;
  }

  /*-------------- 変換表tableによる輝度変換 -------------------*/
  int2D ary1 = pgm1.body.array, ary2 = pgm2.body.array;
  for(int i=0; i< height; i++)
    for(int j = 0; j < width; j++) ary2[i][j] = table[ary1[i][j]];

  /*------ 原画像および変換画像ヒストグラム計算と可視化 ------*/
  visualizeHist(hist, hstpgm1);
  histogram(pgm2.body, hist);
  visualizeHist(hist, hstpgm2);

  /*--------------------- 結果の出力 -------------------------*/
  writePGM(hstpgm1, argv[2]);    //入力画像のヒストグラム
  writePGM(pgm2,    argv[3]);    //平滑化されたPGM画像
  writePGM(hstpgm2, argv[4]);    //平滑化された画像のヒストグラム

  /*-------------------- 作業領域の解放 ----------------------*/
  freeInt1D(table, asize);
  deletePGM(hstpgm2);
  deletePGM(hstpgm1);
  deleteARY1I(hist);
  deletePGM(pgm2);
  deletePGM(pgm1);
  return 0;
}
