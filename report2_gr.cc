#include "iip.h"
/*-------------------------------------------------------------------*/
int main(int argc,char *argv[])
{
const int ImgSize = 81;
double sigma, gamma, theta, lambda;
if (argc!=6 ||(argc>1 && strcmp(argv[1],"help")==0)) {
cout << "\nUsage: " << argv[0] << " sigma gamma theta lambda outImage\n" <<endl;
exit(1);
}
/* 初期スケール係数の読み込み */
sscanf(argv[1], "%lf", &sigma);
sscanf(argv[2], "%lf", &gamma);
sscanf(argv[3], "%lf", &theta);
sscanf(argv[4], "%lf", &lambda);
GRKNL1& gk = createGRKNL1(sigma, gamma, theta, lambda);
PGM& pgm = createPGM(ImgSize, ImgSize); //81x81 の画像領域を生成
ARY2I& img = pgm.body;
setARY2I(img, 0); //配列の初期化
img.array[ImgSize/2][ImgSize/2]= 1; //2 次元インパルス入力信号
ARY2D& ary = createARY2D(ImgSize,ImgSize); //結果格納のための ARY2D 配列
gaborConv1(img, gk, ary);
PGM& outImage = createPGM(ImgSize, ImgSize);
ARY2DtoPGM(ary, outImage);
/*---------------------- 結果の出力と後処理 -----------------------*/
writePGM(outImage, argv[5]);
deletePGM(outImage);
deleteARY2D(ary);
deletePGM(pgm);
deleteGRKNL1(gk);
return 0;
}
