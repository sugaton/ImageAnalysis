#include 	"iip.h"
#include 	"function.h"

const double MAX_HUE = 360.0;
const double MAX_SAT = 1.0;

double2D plus_(ARY2D ary, double val){
  int height = ary.height, width = ary.width;
  double2D d2 = allocDouble2D(height, width);
  for(int i=0; i<height; ++i)
    for(int j=0; j<width; ++j)
      d2[i][j] = ary.array[i][j] + val;
  return d2;
}

struct YCbCr {
  string            ptype;
  ARY2D             Y, Cb, Cr;
  struct YCbCr     *self;
} ;

YCbCr& createYCbCr(int height, int width)
{
  YCbCr* ptr = new YCbCr;
  YCbCr& ycbcr = *ptr;
  ycbcr.ptype = "YCbCr";
  ycbcr.self = ptr;
  ycbcr.Y = createARY2D(height, width);
  ycbcr.Cr = createARY2D(height, width);
  ycbcr.Cb = createARY2D(height, width);
  return ycbcr;
}

void deleteYCbCr(YCbCr& ycbcr)
{
  deleteARY2D(ycbcr.Y); deleteARY2D(ycbcr.Cr); deleteARY2D(ycbcr.Cb);
  delete ycbcr.self;
}

void RGBtoYCbCr(RGB& rgb, YCbCr& ycbcr)
{
  int width = (rgb.red).width,  height = (rgb.red).height;
  if(width != (ycbcr.Y).width || height != (ycbcr.Y).height)
    error("Array size mismatch error in RGBtoYCbCr", 1);

  double2D r = (rgb.red).array,  g = (rgb.green).array, b = (rgb.blue).array;
  double2D y = (ycbcr.Y).array,  cb = (ycbcr.Cb).array,   cr = (ycbcr.Cr).array;

  for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
      double rr = r[i][j], gg = g[i][j], bb = b[i][j];
      y[i][j] = 0.2989 * rr + 0.5866 * gg + 0.1145 * bb;
      cb[i][j] = -0.1687 * rr -0.3313 * gg + 0.5000 * bb;
      cr[i][j] = 0.5000 * rr - 0.4187 * gg - 0.0813 * bb;
    }
  }
}


inline int reScale(double val, double max_val, int max_grade)
{
  return (int)(val * max_grade / max_val);
}

void create2DcolorMap(RGB& rgb, double max1, double2D ary1,
                      double max2, double2D ary2, PPM& map)
/*------------------------------------------------------------*/
/*   この関数では，予めどの色指標(ary1とary2)に対する2次元の  */
/*   色マップであるのか，そして，その2次元指標色マップの縦横  */
/*   のサイズを予め決めておく必要がある．                     */
/*   色指標とは，RGB表色系を他の表色系に変換後に得られる新た  */
/*   な３つの色系のことである．この３つの中から２つを選択．   */
/*   ２次元色指標マップは，PPM型の配列として表現される．      */
/*------------------------------------------------------------*/
{
  int width    = (rgb.red).width, height = (rgb.red).height;
  double2D red = rgb.red.array;
  double2D grn = rgb.green.array;
  double2D blu = rgb.blue.array;
  int2D MR = (map.red).array, MG = (map.green).array, MB = (map.blue).array;

  //2次元指標色マップの縦，横のサイズ
  int Hsize = map.red.width, Vsize = map.red.height;

  //マップ配列の初期化
  for(int i=0; i<Vsize; i++)
    for(int j=0; j<Hsize; j++)
      MR[i][j] = MG[i][j] = MB[i][j] = 0;

  for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
      double index1 = ary1[i][j];
      double index2 = ary2[i][j];
      int n = reScale(index1, max1, Hsize-1);
      int m = reScale(index2, max2, Vsize-1);

      //2次元配列ppmの対応する位置にRGB画像の各画素をマップ
      MR[m][n] = (int)(red[i][j] * 255.0);
      MG[m][n] = (int)(grn[i][j] * 255.0);
      MB[m][n] = (int)(blu[i][j] * 255.0);
    }
  }
}

inline bool flesh_filter_HS(int v, int h, HSV& hsv){
  return ((hsv.hue.array[v][h]) / MAX_HUE) < 0.12 && (hsv.sat.array[v][h]) / MAX_SAT < 0.6;
}
void flesh_filter(HSV& hsv, PGM& pgm){
  int width = hsv.value.width, height = hsv.value.height;
  // double2D cb_ = plus_(hsv.Cb, 0.5);
  for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
      if(flesh_filter_HS(i, j, hsv)){
        pgm.body.array[i][j] = 255;
        // printf("true\n");
      }
      else
        pgm.body.array[i][j] = 0;
        // pgm.body.array[i][j] = (int)(hsv.value.array[i][j]*255);
    }
  }
}

/// 1. に対応 RGB -> HSV , RGB -> YCbCr
void main1_(RGB& rgb, HSV& hsv, YCbCr& ycbcr){
  RGBtoHSV(rgb, hsv);
  RGBtoYCbCr(rgb, ycbcr);
}
/// 2.に対応 HSV->PGM , YCbCr -> PGM
void main2_(HSV& hsv, YCbCr& ycbcr){
  int height = hsv.value.height, width = hsv.value.width;
  PGM pgm1 = createPGM(height, width);
  ARY2DtoPGM(hsv.value, pgm1);
  writePGM(pgm1, "HSV_gray.pgm");
  PGM pgm2 = createPGM(height, width);
  ARY2DtoPGM(ycbcr.Y, pgm2);
  writePGM(pgm2, "YCbCr_gray.pgm");
  deletePGM(pgm1);
  deletePGM(pgm2);
}

/// 3.に対応
void main3_(PPM& ppm){
  // ppm -> hsv & ppm -> YCbCr
  int width = ppm.red.width, height = ppm.red.height;
  RGB& rgb = createRGB(height, width);
  HSV& hsv = createHSV(height, width);
  YCbCr& ycbcr = createYCbCr(height, width);
  PPMtoRGB(ppm, rgb);
  RGBtoHSV(rgb, hsv);
  RGBtoYCbCr(rgb, ycbcr);
  //（H-S）指標色マップ
  PPM& hs_ppm = createPPM(height, width);
  create2DcolorMap(rgb, 360.0, hsv.hue.array, 1.0, hsv.sat.array, hs_ppm);
  writePPM(hs_ppm, "HS_colormap.ppm");
  //（Cb-Cr）指標色マップ
  PPM& cbcr_ppm = createPPM(height, width);
  // YCbCrのCb,Crの値域は[-0.5 ~ 0.5]なので[0 ~ 1.0]にして関数に渡す
  double2D cb_ = plus_(ycbcr.Cb, 0.5);
  double2D cr_ = plus_(ycbcr.Cr, 0.5);
  create2DcolorMap(rgb, 1.0, cb_, 1.0, cr_, cbcr_ppm);
  writePPM(cbcr_ppm, "CbCr_colormap.ppm");
  // delete object created in this function
  freeDouble2D(cb_, height, width);
  freeDouble2D(cr_, height, width);
  deletePPM(hs_ppm);
  deletePPM(cbcr_ppm);
  deleteRGB(rgb);
  deleteHSV(hsv);
  deleteYCbCr(ycbcr);
}
/// 4.に対応
void main4_(PPM& ppm){
  int width = ppm.red.width, height = ppm.red.height;
  RGB& rgb = createRGB(height, width);
  HSV& hsv = createHSV(height, width);
  PGM& pgm = createPGM(height, width);
  PPMtoRGB(ppm, rgb);
  RGBtoHSV(rgb, hsv);
  flesh_filter(hsv, pgm);
  writePGM(pgm, "filtered_.pgm");
  // free
  deleteRGB(rgb);
  deleteHSV(hsv);
  deletePGM(pgm);
}

int main(int argc, char const *argv[]) {
  // load
  char* filename =  "Majorca.ppm";
  PPM& ppm = readPPM(filename);
  int width = ppm.red.width, height = ppm.red.height;
  printf("load Majorca.ppm : (%d, %d)\n", height, width);
  // allocate
  RGB& rgb = createRGB(height, width);
  HSV& hsv = createHSV(height, width);
  YCbCr& ycbcr  = createYCbCr(height, width);
  PPMtoRGB(ppm, rgb);
  // 課題 1. 2.
  main1_(rgb, hsv, ycbcr);
  main2_(hsv, ycbcr);
  // 課題 3. 4.
  // PPM& hand_ppm = readPPM("hand1.ppm");
  // printf("load hand1.ppm : (%d, %d)\n", hand_ppm.red.height, hand_ppm.red.width);
  PPM& hand_ppm = readPPM("hand2.ppm");
  printf("load hand2.ppm : (%d, %d)\n", hand_ppm.red.height, hand_ppm.red.width);
  main3_(hand_ppm);
  main4_(hand_ppm);
  // free
  deleteRGB(rgb);
  deleteHSV(hsv);
  deleteYCbCr(ycbcr);
  deletePPM(ppm);
  deletePPM(hand_ppm);
  return 0;
}
