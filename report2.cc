#include "iip.h"
/*-------------------------------------------------------------------*/
int main(int argc,char *argv[])
{
const int ImgSize = 81;
double sigma, gamma;
if (argc!=5 ||(argc>1 && strcmp(argv[1],"help")==0)) {
cout << "\nUsage: " << argv[0] << " sigma gamma outImage1 outImage2\n" <<endl;
exit(1);
}
/* 初期スケール係数の読み込み */
// sscanf_s(argv[1], "%lf", &sigma);
sscanf(argv[1], "%lf", &sigma);
// sscanf_s(argv[2], "%lf", &gamma);
sscanf(argv[2], "%lf", &gamma);
GSKNL0& gk = createGSKNL0(sigma, gamma);
/* Gauss フィルタレスポンスを画像化するための PGM 画像の生成 */
PGM& pgm = createPGM(ImgSize, ImgSize); //81x81 の画像領域を生成
ARY2I& img = pgm.body;
setARY2I(img, 0); //配列の初期化
img.array[ImgSize/2][ImgSize/2]= 1; //2 次元インパルス入力信号
ARY2D& ary = createARY2D(ImgSize,ImgSize); //結果格納のための ARY2D 配列
gaussConv0(img, gk, ary, 1, 0); //1 次 (Gx) ガウス導関数による畳み込み
/*---- 画像化：信号レベル０を 128 階調レベルに設定 ---*/
PGM& outImage = createPGM(ImgSize, ImgSize);
ARY2DtoPGM(ary, outImage);
/*----------------------２次--------------------*/
gaussConv0(img, gk, ary, 2, 0); //2 次 (Gxx) ガウス導関数による畳み込み
PGM& outImage2 = createPGM(ImgSize, ImgSize);
ARY2DtoPGM(ary, outImage2);
/*---------------------- 結果の出力と後処理 -----------------------*/
writePGM(outImage, argv[3]);
writePGM(outImage2, argv[4]);
deletePGM(outImage);
deletePGM(outImage2);
deleteARY2D(ary);
deletePGM(pgm);
deleteGSKNL0(gk);
return 0;
}
