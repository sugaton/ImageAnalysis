/*---------------------------------------------------------------------*/
/*                     Non Local Means Filter version 1                */
/*                            by H.Nagahashi                           */
/*---------------------------------------------------------------------*/

#include "iip.h"

/*--------------------------------------------------------------*/
/*           �K�E�X�֐��Ɋ�Â��d�݌W���̌���                   */
/*--------------------------------------------------------------*/
void makeGaussianWeights(double1D gauss, int halfTempSize, double sigma)
{
  //���S��f����̋����ɉ������K�E�X�֐��ɂ��d�݌W���̐���
  int s = 0;
  double sigma2 = 2.0 * sigma * sigma;
  for(int i = -halfTempSize; i <= halfTempSize; i++)
    for(int j = -halfTempSize; j <= halfTempSize; j++){
        gauss[s++] = exp(-(double)(i*i+j*j) / sigma2 );
    }
}

/*--------------------------------------------------------------*/
/*  templateSize x templateSize�̈�̃x�N�g�����C�����C���֐�   */
/*--------------------------------------------------------------*/
inline int tempRasterize(int2D pic, int i0, int j0, int htSize, int tSize,
                         int height, int width, double1D z)
{

  if(i0 < htSize || i0 >= height-htSize) return (0);
  if(j0 < htSize || j0 >= width -htSize) return (0);

  int s = 0;
  for(int k = -htSize ; k <= htSize; k++){
    int *ptr = &pic[i0+k][0] + (j0 - htSize);
    for(int m = 0; m < tSize; m++) {
      z[s++] = (double) *(ptr++);
    }
  }
  return (1);
}

/*--------------------------------------------------------------*/
/*               Non Local Means Filter�̖{��                   */
/*--------------------------------------------------------------*/
void NonLocalMeansFilter(ARY2I& ary0, ARY2I& ary1, double sigma_a,
	                     double sigma_h, int tempSize, int windSize)
{
  int  width = ary0.width, height = ary0.height;
  int2D pic0 = ary0.array, pic1   = ary1.array;

  int squareTempSize = tempSize * tempSize;
  int halfTempSize   = tempSize / 2;
  int halfWindSize   = windSize / 2;

  //tempSize x tempSize �̂Q�����̈���P�����z��Ƃ��ď���
  double1D z0    = allocDouble1D(squareTempSize);
  double1D z1    = allocDouble1D(squareTempSize);

  //�K�E�X�֐����P�����z��
  double1D gauss = allocDouble1D(squareTempSize);
  makeGaussianWeights(gauss, halfTempSize, sigma_a);

  //���_���������̈�ŕ�����
  for(int i0 = halfTempSize; i0 <= height - halfTempSize; i0++){
    for(int j0 = halfTempSize; j0 <= width - halfTempSize; j0++){
      if(!tempRasterize(pic0, i0, j0, halfTempSize, tempSize,
               height, width, z0)) continue; //���ړ_�ߖT�̈�̂P������

      //��������̒l�Əd�݂̐ώZ�p�ϐ�
      double new_val = 0.0, total_weights = 0.0;
      //���ړ_�̈���̊e�_�ɑ΂��ē��l�ɋߖT�̈�̂P������
      for(int i1 = i0 - halfWindSize; i1 <= i0 + halfWindSize; i1++){
        for(int j1 = j0 - halfWindSize; j1 <= j0 + halfWindSize; j1++){
          if(!tempRasterize(pic0, i1,j1, halfTempSize, tempSize,
                            height, width, z1)) continue;

          //�d�ݕt���ގ��x�v�Z
          double wdist =0.0, weight = 0.0;
          for(int s = 0; s < squareTempSize; s++)
            wdist += (z0[s]-z1[s]) * (z0[s]-z1[s]) * gauss[s];

          weight = exp(-wdist/(sigma_h * sigma_h));
          total_weights += weight;
          new_val += pic0[i1][j1] * weight;
        }
      }
      new_val /= total_weights;
      pic1[i0][j0] = (int) new_val;
    }
  }
  // �摜�[�̏���
  for(int i = 0; i < height; i++){
	  for(int j = 0; j < halfTempSize; j++) pic1[i][j] = 0;
	  for(int j = width-halfTempSize; j<width; j++) pic1[i][j] = 0;
  }
  for(int i = 0; i < halfTempSize; i++)
	  for(int j = 0; j < width; j++) pic1[i][j] = 0;
  for(int i = height-halfTempSize; i<height; i++)
	  for(int j = 0; j < width; j++) pic1[i][j] = 0;

  freeDouble1D(gauss, squareTempSize);
  freeDouble1D(z1,    squareTempSize);
  freeDouble1D(z0,    squareTempSize);

}

void BiLateralFilter(ARY2I& ary0, ARY2I& ary1, double sigma_c, double sigma_s)
{
  int  width = ary0.width, height = ary0.height;
  int2D pic0 = ary0.array, pic1   = ary1.array;
  // gauss�j�����O�v�Z
  int Smax = 255;
  double1D gauss_s = allocDouble1D(Smax);  // gauss_s[f(x,y)-f(x',y')]
  double2D gauss_c = allocDouble2D(width, height); // gauss_c[x-x'][y-y']
  for (int i = 0; i < Smax; ++i)
    gauss_s[i] = std::exp(-1 * (i * i) / (2 * sigma_s * sigma_s) );
  //�@�e2�_�Ԃ̋����ɂ��d�݂����O�v�Z
  for (int i = 0; i < width; ++i)
    for (int j = 0; j < height; ++j)
      gauss_c[i][j] = std::exp(-1 * ((i * i) + (j * j)) / (2 * sigma_c * sigma_c) );
  // �𒼂ɂ��ׂĂ�2�_�̑g�ݍ��킹�ɂ��ă��[�v��
  for(int y = 0; y < height; ++y){
    for(int x = 0; x < width; ++x){
      // �Ώۓ_(x, y), �C�ӂ̉�f(x_, y_)
      double weight;
      int val = pic0[y][x];
      double regu = 0.0, sum = 0.0;
      for(int y_ = 0; y_ < height; ++y_)
        for(int x_ = 0; x_ < width; ++x_){
          if (y == y_ && x_ == x) continue;
          weight = gauss_c[std::abs(x - x_)][std::abs(y - y_)] * gauss_s[std::abs(pic0[y_][x_] - val)];
          sum += (double)pic0[y_][x_] * weight;
          regu += weight;
        }
      pic1[y][x] = (int)(sum / regu);
    }
  }
  freeDouble1D(gauss_s, Smax);
  freeDouble2D(gauss_c, width, height);

}
/*--------------------------------------------------------------*/
/*              NL-Means�t�B���^�������C���֐�                  */
/*--------------------------------------------------------------*/
int main(int argc, char *argv[] )
{
  if (argc!=6 ||(argc>1 && strcmp(argv[1],"help")==0)) {
    cout << "\nUsage: " << argv[0] <<" inImage sigma_1 sigma_2 outImage type\n";
    cout << "type: \t1: Non-Local Means Filter\n\t2: Bi-lateral Filter\n";
    exit(1);
  }
  int type; // 1: NL-Means    2: BiLateralFilter
  int templateSize;
  int windowSize = 21;  //�����Csigma_a�̒l��3�ȏ�ł���΁C�v�T�C�Y�g��

  //�摜�̓ǂݍ���
  PGM& src = readPGM(argv[1]);
  int width = src.body.width, height = src.body.height;

  //�K�E�X�֐��p�����[�^�̓ǂݍ���
  double sigma_1, sigma_2;
  sscanf(argv[2], "%lf",&sigma_1); // �K�E�X�p�����[�^
  sscanf(argv[3], "%lf",&sigma_2);
  sscanf(argv[5], "%d",&type);

  if(sigma_1 > 3.0)
    error("Sigma_a is too big! Please choose new one less than 3.0.",1);

  if(sigma_1 < 1.0)
	  templateSize = 5;
  else
	  templateSize = ((int)(sigma_1 * 2.5 + 0.5)) * 2 +1;

  // �o�͉摜�p�z��̊m��
  PGM& des = createPGM(height, width);
  switch (type) {
      case 1:
        NonLocalMeansFilter(src.body, des.body, sigma_1, sigma_2,
                            templateSize, windowSize);
      case 2:
        BiLateralFilter(src.body, des.body, sigma_1, sigma_2);
  }

  writePGM(des, argv[4]);

  deletePGM(des);
  deletePGM(src);
  return 0;
}
